import { compileSurface } from '../surface.js';
import { createExample } from '../examples.js';
import {
  KERNEL_NAME,
  KERNEL_VERSION,
  NORMALIZED_DOMAIN,
  NurbsSurfaceKernel
} from './nurbs-kernel.js';
import { buildKernelSurface } from './adapter.js';
import { norm3, sub3, maxComponentError } from '../math-utils.js';

export const ARCHIVE_SCHEMA_VERSION = 2;
export const DECAL_BAKED_SEAM = 'baked-seam';
export const DECAL_ARC_LENGTH = 'arc-length';
export const MIGRATION_TARGET_DOMAIN = Object.freeze([2, 5]);

function arcLength3D(points) {
  let total = 0;
  const cumulative = [0];
  for (let i = 1; i < points.length; i++) {
    total += norm3(sub3(points[i], points[i - 1]));
    cumulative.push(total);
  }
  return { total, cumulative };
}

function makeLoopParameterSet(surface) {
  const [u0, u1] = surface.domainU;
  const [v0, v1] = surface.domainV;
  const uMid = (u0 + u1) / 2;
  const seamU = u0 + 0.35 * (u1 - u0);
  const params = [];
  const N = 48;
  for (let k = 0; k < N; k++) {
    const theta = (2 * Math.PI * k) / N;
    const cu = Math.cos(theta);
    const su = Math.sin(theta);
    params.push({
      u: uMid + 0.28 * (u1 - u0) * cu,
      v: (v0 + v1) / 2 + 0.32 * (v1 - v0) * su,
      seam: theta === 0 ? { u: seamU, v: (v0 + v1) / 2 } : null
    });
  }
  return { params, seamU };
}

function pointFromEval(surface, u, v) {
  const r = surface.evaluate(u, v);
  if (!r.ok) throw new Error(`TRIM_FIXTURE_EVAL_FAILED: ${r.code}`);
  return r.position;
}

function buildTrimLoop(surface) {
  const { params, seamU } = makeLoopParameterSet(surface);
  const points3D = params.map(p => pointFromEval(surface, p.u, p.v));
  const { total, cumulative } = arcLength3D(points3D);
  return {
    params: params.map(p => ({ u: p.u, v: p.v })),
    points3D,
    arcLength: { total, cumulative },
    seamU
  };
}

function buildDecals(surface, loop, mode) {
  const count = 8;
  const decals = [];
  for (let i = 0; i < count; i++) {
    const index = Math.round((i * (loop.points3D.length - 1)) / (count - 1));
    const targetFraction = i / (count - 1);
    const base = {
      id: `decal-${i}`,
      index,
      fraction: targetFraction,
      placement3D: loop.points3D[index].slice()
    };
    if (mode === DECAL_BAKED_SEAM) {
      decals.push({
        ...base,
        mode: DECAL_BAKED_SEAM,
        seamProvenance: {
          seamU: loop.seamU,
          u: loop.params[index].u,
          v: loop.params[index].v,
          uv: [
            (loop.params[index].u - loop.seamU) / (surface.domainU[1] - surface.domainU[0]),
            loop.params[index].v
          ]
        }
      });
    } else {
      const located = pointAtArcFraction(loop, targetFraction);
      decals.push({
        ...base,
        index: located.index,
        placement3D: located.point,
        mode: DECAL_ARC_LENGTH,
        arcLengthFraction: targetFraction
      });
    }
  }
  return decals;
}

export function createTrimmedArchive(mode, options = {}) {
  const base = options.surfaceSpec ?? createExample(options.exampleId ?? 'cylinder');
  const spec = structuredClone(base);
  const surface = compileSurface(spec);
  const loop = buildTrimLoop(surface);
  const decals = buildDecals(surface, loop, mode);
  const archive = {
    format: 'rational-spline-studio-trimmed',
    schemaVersion: 1,
    kernel: { name: 'legacy-basis-sum', version: '1.0.0' },
    id: `trim-${mode}`,
    surfaceSpec: spec,
    trimLoops: [
      {
        id: 'loop-0',
        domain: 'legacy',
        params: loop.params,
        points3D: loop.points3D,
        arcLength: loop.arcLength
      }
    ],
    decals,
    seamSource: mode === DECAL_BAKED_SEAM
      ? { seamU: loop.seamU, bakedInto: 'decal-uv' }
      : null
  };
  if (options.omitSeamSource) archive.seamSource = null;
  return archive;
}

function kernelFromArchive(archive, targetDomain = NORMALIZED_DOMAIN) {
  const surface = compileSurface(archive.surfaceSpec);
  const bundle = buildKernelSurface(archive.surfaceSpec, { reparametrize: true, targetDomain });
  return { surface, kernel: bundle.kernel, transforms: bundle.transforms };
}

function invertLandmark(kernel, target, hintU, hintV, scale) {
  const guess = {
    u: hintU === null ? null : hintU,
    v: hintV === null ? null : hintV
  };
  const result = kernel.invertPoint(target, guess.u, guess.v, { scale, tolerance: 1e-10, iterations: 64 });
  return result;
}

function migrateTrimLoop(archive, kernel, transforms, surface, trace) {
  const loops = [];
  for (const loop of archive.trimLoops) {
    const params = [];
    const inversionResiduals = [];
    for (let i = 0; i < loop.points3D.length; i++) {
      const target = loop.points3D[i];
      const oldParam = loop.params[i];
      const hintU = transforms.U.map(oldParam.u);
      const hintV = transforms.V.map(oldParam.v);
      const inv = invertLandmark(kernel, target, hintU, hintV, surface.scale);
      if (!inv.ok) {
        trace.push({ code: 'LANDMARK_INVERSION_FAILED', loop: loop.id, index: i, detail: inv.code });
        return { ok: false, code: 'LANDMARK_INVERSION_FAILED', trace };
      }
      inversionResiduals.push(inv.residual);
      params.push({ u: inv.u, v: inv.v, residual: inv.residual, iterations: inv.iterations });
    }
    const points3D = params.map(p => {
      const r = kernel.evaluate(p.u, p.v);
      return r.position;
    });
    const { total, cumulative } = arcLength3D(points3D);
    loops.push({
      id: loop.id,
      domain: 'kernel-v2-normalized',
      params: params.map(p => ({ u: p.u, v: p.v })),
      points3D,
      arcLength: { total, cumulative },
      inversionResiduals
    });
  }
  return { ok: true, loops, trace };
}

function migrateBakedSeamDecal(decal, transforms) {
  const sp = decal.seamProvenance;
  return {
    id: decal.id,
    mode: DECAL_BAKED_SEAM,
    index: decal.index,
    fraction: decal.fraction,
    placement3D: decal.placement3D.slice(),
    seamProvenance: {
      seamU: transforms.U.map(sp.seamU),
      u: transforms.U.map(sp.u),
      v: transforms.V.map(sp.v),
      uv: sp.uv.slice()
    }
  };
}

function lerp3(a, b, t) {
  return [
    a[0] + (b[0] - a[0]) * t,
    a[1] + (b[1] - a[1]) * t,
    a[2] + (b[2] - a[2]) * t
  ];
}

function pointAtArcFraction(loop, fraction) {
  const targetLength = fraction * loop.arcLength.total;
  const cumulative = loop.arcLength.cumulative;
  let chosen = 0;
  for (let i = 1; i < cumulative.length; i++) {
    if (cumulative[i] <= targetLength) chosen = i;
  }
  const segStart = cumulative[chosen];
  const segEnd = chosen + 1 < cumulative.length ? cumulative[chosen + 1] : segStart;
  const segLen = segEnd - segStart;
  const local = segLen === 0 ? 0 : (targetLength - segStart) / segLen;
  const nextIndex = Math.min(chosen + 1, loop.points3D.length - 1);
  return {
    index: chosen,
    point: lerp3(loop.points3D[chosen], loop.points3D[nextIndex], local),
    param: loop.params[chosen]
  };
}

function migrateArcLengthDecal(decal, migratedLoop, kernel, trace) {
  const located = pointAtArcFraction(migratedLoop, decal.arcLengthFraction);
  const evalResult = kernel.evaluate(located.param.u, located.param.v);
  if (!evalResult.ok) {
    trace.push({ code: 'ARC_DECAL_PLACEMENT_FAILED', decal: decal.id, detail: evalResult.code });
    return null;
  }
  return {
    id: decal.id,
    mode: DECAL_ARC_LENGTH,
    index: located.index,
    fraction: decal.arcLengthFraction,
    arcLengthFraction: decal.arcLengthFraction,
    placement3D: located.point
  };
}

export function migrateArchive(archive, options = {}) {
  const trace = [];
  const targetDomain = options.targetDomain ?? MIGRATION_TARGET_DOMAIN;
  const needsSeam = archive.decals.some(d => d.mode === DECAL_BAKED_SEAM);
  if (needsSeam && !archive.seamSource) {
    return {
      ok: false,
      status: 'SUSPENDED',
      code: 'MISSING_SEAM_SOURCE',
      archiveId: archive.id,
      message: '旧档包含接缝烘入贴花但缺失 seamSource；单独挂起，等待补齐接缝来源后再迁移。',
      trace: [{ code: 'MISSING_SEAM_SOURCE', message: '无法在无接缝来源的情况下重定位烘入贴花。' }]
    };
  }

  const { surface, kernel, transforms } = kernelFromArchive(archive, targetDomain);
  const loopResult = migrateTrimLoop(archive, kernel, transforms, surface, trace);
  if (!loopResult.ok) return { ok: false, status: 'FAILED', ...loopResult };

  const migratedLoops = loopResult.loops;
  const decals = [];
  for (const decal of archive.decals) {
    if (decal.mode === DECAL_BAKED_SEAM) {
      decals.push(migrateBakedSeamDecal(decal, transforms));
    } else {
      const migrated = migrateArcLengthDecal(decal, migratedLoops[0], kernel, trace);
      if (!migrated) return { ok: false, status: 'FAILED', code: 'ARC_DECAL_PLACEMENT_FAILED', trace };
      decals.push(migrated);
    }
  }

  const migrated = {
    format: archive.format,
    schemaVersion: ARCHIVE_SCHEMA_VERSION,
    kernel: { name: KERNEL_NAME, version: KERNEL_VERSION, parameterDomain: targetDomain.slice() },
    id: archive.id,
    surfaceSpec: reparametrizedSpec(archive.surfaceSpec, transforms),
    sourceSpec: archive.surfaceSpec,
    transforms,
    targetDomain: targetDomain.slice(),
    trimLoops: migratedLoops,
    decals,
    seamSource: archive.seamSource
      ? { seamU: transforms.U.map(archive.seamSource.seamU), bakedInto: archive.seamSource.bakedInto }
      : null
  };
  return { ok: true, status: 'MIGRATED', archive: migrated, trace };
}

function reparametrizedSpec(spec, transforms) {
  return {
    ...structuredClone(spec),
    knotsU: spec.knotsU.map(transforms.U.map),
    knotsV: spec.knotsV.map(transforms.V.map)
  };
}

export function resumeSuspendedArchive(archive, seamSource) {
  const restored = structuredClone(archive);
  restored.seamSource = seamSource;
  return migrateArchive(restored);
}

function decalIntersections(decals, loops) {
  return decals.map(decal => {
    const loop = loops[0];
    const idx = Math.min(decal.index, loop.points3D.length - 1);
    const point = decal.placement3D ? decal.placement3D.slice() : loop.points3D[idx].slice();
    return {
      id: decal.id,
      mode: decal.mode,
      index: idx,
      point,
      param: loop.params[idx]
    };
  });
}

function trimGrid(loop, samples = 24) {
  const vertices = [];
  const n = loop.points3D.length;
  for (let k = 0; k < samples; k++) {
    const idx = Math.round((k * (n - 1)) / (samples - 1));
    vertices.push(loop.points3D[idx].slice());
  }
  const edges = vertices.slice(1).map((_, k) => [k, k + 1]);
  return { vertices, edges };
}

export function checkRoundTrip(archive, migrated, options = {}) {
  const positionTol = options.positionTol ?? 1e-8;
  const gridTol = options.gridTol ?? 1e-8;
  const oldSurface = compileSurface(archive.surfaceSpec);
  const newBundle = buildKernelSurface(migrated.surfaceSpec, { reparametrize: false });
  const kernel = newBundle.kernel;

  const oldIntersections = decalIntersections(archive.decals, archive.trimLoops);
  const newIntersections = decalIntersections(migrated.decals, migrated.trimLoops);
  const decalErrors = [];
  let decalsReversible = true;
  for (let i = 0; i < oldIntersections.length; i++) {
    const a = oldIntersections[i];
    const b = newIntersections[i];
    const error = maxComponentError(a.point, b.point);
    const modeMatch = a.mode === b.mode;
    if (error > positionTol || !modeMatch) decalsReversible = false;
    decalErrors.push({ id: a.id, mode: a.mode, positionError: error, modePreserved: modeMatch });
  }

  const oldGrid = trimGrid(archive.trimLoops[0]);
  const newGrid = trimGrid(migrated.trimLoops[0]);
  let gridMaxError = 0;
  const gridErrors = [];
  for (let i = 0; i < oldGrid.vertices.length; i++) {
    const error = maxComponentError(oldGrid.vertices[i], newGrid.vertices[i]);
    gridErrors.push(error);
    gridMaxError = Math.max(gridMaxError, error);
  }

  const inverseLoopErrors = [];
  let inverseReversible = true;
  for (let i = 0; i < migrated.trimLoops[0].params.length; i++) {
    const newParam = migrated.trimLoops[0].params[i];
    const target = migrated.trimLoops[0].points3D[i];
    const oldU = migrated.transforms.U.inverse(newParam.u);
    const oldV = migrated.transforms.V.inverse(newParam.v);
    const back = oldSurface.evaluate(oldU, oldV);
    if (!back.ok) {
      inverseReversible = false;
      inverseLoopErrors.push({ index: i, error: Infinity, code: back.code });
      continue;
    }
    const error = maxComponentError(back.position, target);
    if (error > positionTol) inverseReversible = false;
    inverseLoopErrors.push({ index: i, error });
  }

  return {
    ok: decalsReversible && inverseReversible && gridMaxError <= gridTol,
    reversible: decalsReversible && inverseReversible,
    decalsReversible,
    inverseReversible,
    gridMaxError,
    gridWithinTolerance: gridMaxError <= gridTol,
    decalErrors,
    inverseLoopErrors,
    tolerances: { positionTol, gridTol },
    kernel: { name: KERNEL_NAME, version: KERNEL_VERSION }
  };
}

export function compareMigrationStrategies(archive, options = {}) {
  const targetDomain = options.targetDomain ?? MIGRATION_TARGET_DOMAIN;
  const adapted = buildKernelSurface(archive.surfaceSpec, { reparametrize: false });
  const reparam = buildKernelSurface(archive.surfaceSpec, { reparametrize: true, targetDomain });
  const loop = archive.trimLoops[0];
  const midIndex = Math.floor(loop.params.length / 2);
  const probe = loop.params[midIndex];
  const target = loop.points3D[midIndex];

  const adaptedEval = adapted.kernel.evaluate(probe.u, probe.v);
  const reParam = { u: reparam.transforms.U.map(probe.u), v: reparam.transforms.V.map(probe.v) };
  const reparamEval = reparam.kernel.evaluate(reParam.u, reParam.v);

  const adaptedInverse = adapted.kernel.invertPoint(target, probe.u, probe.v, { scale: 1 });
  const reparamInverse = reparam.kernel.invertPoint(target, reParam.u, reParam.v, { scale: 1 });

  const sameGeometry = adaptedEval.ok && reparamEval.ok
    && maxComponentError(adaptedEval.position, reparamEval.position) < 1e-10;

  return {
    targetDomain: targetDomain.slice(),
    strategyA: {
      name: 'keep-legacy-domain-adapt-new-kernel',
      domainU: adapted.kernel.domainU,
      domainV: adapted.kernel.domainV,
      positionError: adaptedEval.ok ? maxComponentError(adaptedEval.position, target) : null,
      inversion: adaptedInverse
    },
    strategyB: {
      name: 'reparametrize-source-and-trim-loops',
      domainU: reparam.kernel.domainU,
      domainV: reparam.kernel.domainV,
      positionError: reparamEval.ok ? maxComponentError(reparamEval.position, target) : null,
      inversion: reparamInverse
    },
    geometryIdenticalAcrossDomains: sameGeometry,
    domainsDiffer: JSON.stringify(adapted.kernel.domainU) !== JSON.stringify(reparam.kernel.domainU)
      || JSON.stringify(adapted.kernel.domainV) !== JSON.stringify(reparam.kernel.domainV),
    selected: 'strategyB-reparametrize',
    reason: '裁剪曲线与贴花必须随节点域统一换域，否则同形几何的两套贴花会被按同一参数位置错误处理。'
  };
}

export { kernelFromArchive, NurbsSurfaceKernel };
