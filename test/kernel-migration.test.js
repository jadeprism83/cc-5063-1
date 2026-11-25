import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createExample } from '../src/examples.js';
import { compileSurface } from '../src/surface.js';
import { maxComponentError } from '../src/math-utils.js';
import {
  KERNEL_NAME,
  KERNEL_VERSION,
  NORMALIZED_DOMAIN,
  NurbsCurveKernel,
  kernelBasisValues,
  kernelBasisDerivatives,
  knotSpanIndex
} from '../src/kernel/nurbs-kernel.js';
import { buildKernelSurface } from '../src/kernel/adapter.js';
import {
  createTrimmedArchive,
  migrateArchive,
  resumeSuspendedArchive,
  checkRoundTrip,
  compareMigrationStrategies,
  MIGRATION_TARGET_DOMAIN,
  ARCHIVE_SCHEMA_VERSION,
  DECAL_BAKED_SEAM,
  DECAL_ARC_LENGTH
} from '../src/kernel/trimming.js';

const POSITION_TOL = 1e-8;

function maxPos(a, b) {
  return maxComponentError(a, b);
}

function derivativeNormError(actual, expected) {
  const num = Math.hypot(...actual.map((x, i) => x - expected[i]));
  const den = Math.max(1, Math.hypot(...expected));
  return num / den;
}

test('build-lock and package kernelLock pin the new kernel and its versions', () => {
  const lock = JSON.parse(readFileSync(new URL('../build-lock.json', import.meta.url)));
  const pkg = JSON.parse(readFileSync(new URL('../package.json', import.meta.url)));
  assert.equal(lock.engine.name, KERNEL_NAME);
  assert.equal(lock.engine.version, KERNEL_VERSION);
  assert.equal(lock.evaluationPath.primary, KERNEL_NAME);
  assert.equal(lock.evaluationPath.meshGeneration, KERNEL_NAME);
  assert.equal(pkg.kernelLock.engine, KERNEL_NAME);
  assert.equal(pkg.kernelLock.engineVersion, KERNEL_VERSION);
  assert.equal(pkg.version, KERNEL_VERSION);
  assert.deepEqual(lock.migration.defaultTargetDomain, [...MIGRATION_TARGET_DOMAIN]);
});

test('surface evaluation runs through the new kernel path, not the legacy engine', () => {
  for (const id of ['plane', 'cylinder', 'nonuniform', 'repeated-heavy']) {
    const s = compileSurface(createExample(id));
    assert.equal(s.engine, KERNEL_NAME);
    assert.equal(s.kernelVersion, KERNEL_VERSION);
    const r = s.evaluate(0.37, 0.61);
    assert.equal(r.ok, true, `${id}: ${r.code}`);
    assert.equal(r.engine, KERNEL_NAME);
    assert.equal(r.branch.method, 'de-boor-kernel-v2');
    assert.ok(r.kernelParameters);
  }
});

test('new kernel and legacy basis-sum agree on position and derivatives across all examples', () => {
  const fractions = [0, 0.13, 0.5, 0.5001, 0.87, 1];
  for (const id of ['plane', 'cylinder', 'nonuniform', 'repeated-heavy', 'zero-weight-endpoint']) {
    const s = compileSurface(createExample(id));
    const [u0, u1] = s.domainU;
    const [v0, v1] = s.domainV;
    let worstPosition = 0;
    let worstDerivative = 0;
    for (const f of fractions) {
      for (const g of [0, 0.37, 1]) {
        const u = u0 + f * (u1 - u0);
        const v = v0 + g * (v1 - v0);
        const neu = s.evaluate(u, v);
        const legacy = s.evaluateLegacy(u, v);
        if (!neu.ok || !legacy.ok) {
          assert.equal(neu.ok, legacy.ok, `${id} status mismatch at ${u},${v}`);
          assert.equal(neu.code, legacy.code);
          continue;
        }
        worstPosition = Math.max(worstPosition, maxPos(neu.position, legacy.position));
        worstDerivative = Math.max(
          worstDerivative,
          derivativeNormError(neu.derivatives.u, legacy.derivatives.u),
          derivativeNormError(neu.derivatives.v, legacy.derivatives.v)
        );
      }
    }
    assert.ok(worstPosition < 1e-11, `${id} position error ${worstPosition}`);
    assert.ok(worstDerivative < 1e-10, `${id} derivative error ${worstDerivative}`);
  }
});

test('repeated knots resolve on both one-sided branches and stay a partition of unity', () => {
  const knots = [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1];
  const degree = 3;
  for (const t of [0, 0.25, 0.5, 0.75, 1]) {
    const sides = t === 0 ? ['R'] : t === 1 ? ['L'] : ['L', 'R'];
    for (const side of sides) {
      const values = kernelBasisValues(knots, degree, t, side);
      const sum = values.reduce((a, x) => a + x, 0);
      assert.ok(Math.abs(sum - 1) < 1e-12, `t=${t} side=${side} sum=${sum}`);
    }
  }
  // A double interior knot in degree 3 is C1: basis value AND derivative stay
  // continuous across sides, while the resolved knot span switches one-sidedly.
  const leftV = kernelBasisValues(knots, degree, 0.5, 'L');
  const rightV = kernelBasisValues(knots, degree, 0.5, 'R');
  assert.ok(leftV.every((v, i) => Math.abs(v - rightV[i]) < 1e-12));
  const leftD = kernelBasisDerivatives(knots, degree, 0.5, 'L');
  const rightD = kernelBasisDerivatives(knots, degree, 0.5, 'R');
  assert.ok(leftD.every((v, i) => Math.abs(v - rightD[i]) < 1e-12));
  assert.notEqual(knotSpanIndex(knots, degree, 0.5, 'L'), knotSpanIndex(knots, degree, 0.5, 'R'));
});

test('extreme-weight surface still matches the legacy engine within tolerance', () => {
  const s = compileSurface(createExample('repeated-heavy'));
  let worst = 0;
  for (const f of [0, 0.25, 0.4999, 0.5, 0.5001, 0.75, 1]) {
    const u = f;
    const neu = s.evaluate(u, 0.5);
    const legacy = s.evaluateLegacy(u, 0.5);
    assert.equal(neu.ok, legacy.ok);
    if (neu.ok) worst = Math.max(worst, maxPos(neu.position, legacy.position));
  }
  assert.ok(worst < 1e-9, `extreme weight position error ${worst}`);
});

test('zero-weight degenerate input returns an explicit status, never a fabricated point', () => {
  const s = compileSurface(createExample('zero-weight-endpoint'));
  const r = s.evaluate(s.domainU[0], s.domainV[0]);
  assert.equal(r.ok, false);
  assert.equal(r.code, 'DEGENERATE_DENOMINATOR');
  assert.equal(r.engine, KERNEL_NAME);
  assert.ok(r.limitCandidates.alongU || r.limitCandidates.alongV);
  assert.equal(r.position, undefined);
  assert.equal(r.normal, undefined);
});

test('curve kernel rejects malformed degree and knot length with explicit codes', () => {
  assert.throws(
    () => new NurbsCurveKernel({ degree: -1, knots: [0, 1], homogeneousPoints: [[1, 1]] }),
    /KERNEL_INVALID_DEGREE/
  );
  assert.throws(
    () => new NurbsCurveKernel({ degree: 2, knots: [0, 1], homogeneousPoints: [[1, 1]] }),
    /KERNEL_KNOT_LENGTH_MISMATCH/
  );
});

test('reparametrization maps knots and jacobian consistently while preserving geometry', () => {
  const spec = createExample('nonuniform');
  const legacy = buildKernelSurface(spec, { reparametrize: false });
  const target = [2, 5];
  const moved = buildKernelSurface(spec, { reparametrize: true, targetDomain: target });
  assert.deepEqual(moved.kernel.domainU, target);
  assert.deepEqual(moved.kernel.domainV, target);
  assert.notDeepEqual(legacy.kernel.domainU, moved.kernel.domainU);
  assert.equal(moved.transforms.U.jacobian, target[1] - target[0]);

  const u = 0.43;
  const v = 0.6;
  const a = legacy.kernel.evaluate(u, v);
  const ku = moved.transforms.U.map(u);
  const kv = moved.transforms.V.map(v);
  const b = moved.kernel.evaluate(ku, kv);
  assert.ok(a.ok && b.ok);
  assert.ok(maxPos(a.position, b.position) < 1e-12);
  // inverse map round-trips the parameter exactly
  assert.ok(Math.abs(moved.transforms.U.inverse(ku) - u) < 1e-14);
});

test('spatial landmark inversion recovers parameters on the reparametrized kernel', () => {
  const spec = createExample('cylinder');
  const moved = buildKernelSurface(spec, { reparametrize: true, targetDomain: MIGRATION_TARGET_DOMAIN });
  const s = compileSurface(spec);
  for (const [u, v] of [[0.2, 0.3], [0.55, 0.8], [0.9, 0.1]]) {
    const target = s.evaluate(u, v).position;
    const ku = moved.transforms.U.map(u);
    const kv = moved.transforms.V.map(v);
    const inv = moved.kernel.invertPoint(target, ku, kv, { scale: s.scale, tolerance: 1e-10 });
    assert.equal(inv.ok, true, `inversion failed at ${u},${v}: ${inv.code}`);
    assert.ok(Math.abs(inv.u - ku) < 1e-7 && Math.abs(inv.v - kv) < 1e-7);
    const back = moved.kernel.evaluate(inv.u, inv.v);
    assert.ok(maxPos(back.position, target) < POSITION_TOL);
  }
});

test('both decal modes migrate to the new domain and remain reversibly round-trippable', () => {
  for (const mode of [DECAL_BAKED_SEAM, DECAL_ARC_LENGTH]) {
    const archive = createTrimmedArchive(mode);
    const migration = migrateArchive(archive);
    assert.equal(migration.ok, true, `${mode}: ${migration.code}`);
    const migrated = migration.archive;
    assert.equal(migrated.schemaVersion, ARCHIVE_SCHEMA_VERSION);
    assert.equal(migrated.kernel.name, KERNEL_NAME);
    assert.deepEqual(migrated.kernel.parameterDomain, [...MIGRATION_TARGET_DOMAIN]);
    assert.ok(migrated.decals.every(d => d.mode === mode));

    const rt = checkRoundTrip(archive, migrated);
    assert.ok(rt.decalErrors.every(e => e.modePreserved), `${mode} decal mode not preserved`);
    assert.ok(Math.max(...rt.decalErrors.map(e => e.positionError)) <= POSITION_TOL);
    assert.ok(rt.gridMaxError <= rt.tolerances.gridTol);
    assert.ok(Math.max(...rt.inverseLoopErrors.map(e => e.error)) <= POSITION_TOL);
    assert.equal(rt.ok, true);
    assert.equal(rt.reversible, true);
  }
});

test('same-shaped geometry does not collapse the two decal sets into identical handling', () => {
  const baked = migrateArchive(createTrimmedArchive(DECAL_BAKED_SEAM)).archive;
  const arc = migrateArchive(createTrimmedArchive(DECAL_ARC_LENGTH)).archive;
  assert.equal(baked.decals[0].mode, DECAL_BAKED_SEAM);
  assert.equal(arc.decals[0].mode, DECAL_ARC_LENGTH);
  assert.ok(baked.decals[0].seamProvenance);
  assert.equal(arc.decals[0].seamProvenance, undefined);
  assert.ok(arc.decals[0].arcLengthFraction !== undefined);
  // identical 3D placement can coincide, but provenance/param semantics must stay distinct
  assert.notDeepEqual(
    baked.decals.map(d => d.seamProvenance.uv),
    arc.decals.map(d => [d.arcLengthFraction, d.index])
  );
});

test('archives missing seam source are individually suspended and can be resumed', () => {
  const orphan = createTrimmedArchive(DECAL_BAKED_SEAM, { omitSeamSource: true });
  const suspended = migrateArchive(orphan);
  assert.equal(suspended.ok, false);
  assert.equal(suspended.status, 'SUSPENDED');
  assert.equal(suspended.code, 'MISSING_SEAM_SOURCE');
  assert.equal(suspended.archiveId, orphan.id);

  const resumed = resumeSuspendedArchive(orphan, { seamU: 0.35, bakedInto: 'decal-uv' });
  assert.equal(resumed.ok, true);
  const rt = checkRoundTrip(orphan, resumed.archive);
  assert.equal(rt.ok, true);
  assert.ok(rt.gridMaxError <= rt.tolerances.gridTol);
});

test('strategy comparison shows genuinely different domains with identical geometry', () => {
  const archive = createTrimmedArchive(DECAL_ARC_LENGTH, { exampleId: 'nonuniform' });
  const result = compareMigrationStrategies(archive);
  assert.deepEqual(result.strategyA.domainU, NORMALIZED_DOMAIN);
  assert.deepEqual(result.strategyB.domainU, [...MIGRATION_TARGET_DOMAIN]);
  assert.equal(result.domainsDiffer, true);
  assert.equal(result.geometryIdenticalAcrossDomains, true);
  assert.equal(result.selected, 'strategyB-reparametrize');
  assert.ok(result.strategyB.positionError < POSITION_TOL);
});

test('display mesh is generated by the new kernel and reports its engine on every vertex', () => {
  const s = compileSurface(createExample('cylinder'));
  const mesh = s.tessellate(12, 8);
  assert.equal(mesh.engine, KERNEL_NAME);
  assert.equal(mesh.kernelVersion, KERNEL_VERSION);
  assert.ok(mesh.indices.length > 0);
  assert.equal(mesh.failures.length, 0);
  assert.ok(mesh.evaluations.every(r => r.engine === KERNEL_NAME));
});

test('mutating one old control line still yields a reversible migration on real intersections', () => {
  // Perturb one control line (an entire U-column) programmatically — no hand-edited data files.
  const spec = createExample('cylinder');
  const col = 1;
  for (let i = 0; i < spec.points.length; i++) {
    spec.points[i][col].y += 0.35;
  }
  // Build the archive from the mutated spec so its trim loop/decals are real geometry.
  const archive = createTrimmedArchive(DECAL_ARC_LENGTH, { surfaceSpec: spec });
  const migration = migrateArchive(archive);
  assert.equal(migration.ok, true, migration.code);
  const rt = checkRoundTrip(archive, migration.archive);
  assert.equal(rt.ok, true);
  assert.ok(Math.max(...rt.decalErrors.map(e => e.positionError)) <= POSITION_TOL);
  assert.ok(rt.gridMaxError <= rt.tolerances.gridTol);
});
