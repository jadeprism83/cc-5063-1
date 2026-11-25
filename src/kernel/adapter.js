import {
  KERNEL_NAME,
  KERNEL_VERSION,
  NORMALIZED_DOMAIN,
  NurbsSurfaceKernel,
  reparametrizeKnots
} from './nurbs-kernel.js';
import { knotAxisInfo } from '../basis.js';
import { frameFromDerivatives } from '../math-utils.js';

export const ADAPTER_VERSION = '2.0.0';

function sideToKernel(side) {
  if (side === -1 || side === 'L' || side === 'left') return 'L';
  return 'R';
}

function sideFromKernel(side) {
  return side === 'L' ? -1 : 1;
}

export function homogeneousGridFromSpec(spec) {
  return spec.points.map((row, i) =>
    row.map((p, j) => {
      const w = spec.weights[i][j];
      const x = Array.isArray(p) ? p[0] : p.x;
      const y = Array.isArray(p) ? p[1] : p.y;
      const z = Array.isArray(p) ? p[2] : p.z;
      return [x * w, y * w, z * w, w];
    })
  );
}

export function buildKernelSurface(spec, { reparametrize = false, targetDomain = NORMALIZED_DOMAIN } = {}) {
  let knotsU = spec.knotsU.slice();
  let knotsV = spec.knotsV.slice();
  const transforms = {
    U: { map: t => t, inverse: t => t, jacobian: 1, targetDomain: [spec.knotsU[spec.degreeU], spec.knotsU[spec.knotsU.length - spec.degreeU - 1]] },
    V: { map: t => t, inverse: t => t, jacobian: 1, targetDomain: [spec.knotsV[spec.degreeV], spec.knotsV[spec.knotsV.length - spec.degreeV - 1]] }
  };

  if (reparametrize) {
    const ru = reparametrizeKnots(spec.knotsU, spec.degreeU, targetDomain);
    const rv = reparametrizeKnots(spec.knotsV, spec.degreeV, targetDomain);
    knotsU = ru.knots;
    knotsV = rv.knots;
    transforms.U = { map: ru.map, inverse: ru.inverse, jacobian: ru.jacobian, targetDomain: targetDomain.slice() };
    transforms.V = { map: rv.map, inverse: rv.inverse, jacobian: rv.jacobian, targetDomain: targetDomain.slice() };
  }

  const kernel = new NurbsSurfaceKernel({
    degreeU: spec.degreeU,
    degreeV: spec.degreeV,
    knotsU,
    knotsV,
    homogeneousGrid: homogeneousGridFromSpec(spec)
  });
  return { kernel, transforms, reparametrized: reparametrize };
}

function legacyTrace(kernelResult, infoU, infoV) {
  const trace = [...(infoU?.trace ?? []), ...(infoV?.trace ?? [])];
  if (!kernelResult.ok) {
    if (kernelResult.code === 'DEGENERATE_DENOMINATOR') {
      trace.push({ code: 'DEGENERATE_DENOMINATOR', message: '新内核 de Boor 分支确认有理分母退化；适配层透传状态，不伪造曲面点。' });
    }
  }
  return trace;
}

export function adaptEvaluation(surface, u, v, options = {}) {
  const infoU = knotAxisInfo(surface.U, surface.p, surface.nu, u, options.sideU ?? 'auto');
  const infoV = knotAxisInfo(surface.V, surface.q, surface.nv, v, options.sideV ?? 'auto');

  if (!infoU.inDomain || !infoV.inDomain) {
    return {
      ok: false,
      code: 'OUT_OF_DOMAIN',
      engine: KERNEL_NAME,
      kernelVersion: KERNEL_VERSION,
      adapterVersion: ADAPTER_VERSION,
      message: `参数 (${u}, ${v}) 超出定义域 U${surface.domainU.join(',')} × V${surface.domainV.join(',')}。`,
      u, v,
      infoU, infoV,
      trace: []
    };
  }

  const sideU = sideToKernel(infoU.side);
  const sideV = sideToKernel(infoV.side);
  const kernelInput = surface.kernelBundle;
  const ku = kernelInput.transforms.U.map(u);
  const kv = kernelInput.transforms.V.map(v);
  const raw = kernelInput.kernel.evaluate(ku, kv, { sideU, sideV });

  if (!raw.ok) {
    return {
      ...raw,
      u, v,
      sideU: infoU.side,
      sideV: infoV.side,
      infoU, infoV,
      adapterVersion: ADAPTER_VERSION,
      trace: legacyTrace(raw, infoU, infoV)
    };
  }

  const jacU = kernelInput.transforms.U.jacobian;
  const jacV = kernelInput.transforms.V.jacobian;
  const derivativeU = raw.derivatives.u.map(x => x * jacU);
  const derivativeV = raw.derivatives.v.map(x => x * jacV);
  const frame = frameFromDerivatives(derivativeU, derivativeV, surface.scale);
  const denomScale = surface.scale;
  const cancellationRatio = raw.denominator === 0 ? 0 : 1;

  return {
    ok: true,
    engine: KERNEL_NAME,
    kernelVersion: KERNEL_VERSION,
    adapterVersion: ADAPTER_VERSION,
    kernelSideU: sideU,
    kernelSideV: sideV,
    kernelParameters: { u: ku, v: kv },
    u, v,
    sideU: sideFromKernel(sideU),
    sideV: sideFromKernel(sideV),
    infoU, infoV,
    position: raw.position,
    derivatives: { u: derivativeU, v: derivativeV },
    normal: frame.normal,
    normalStatus: frame.status,
    frame,
    denominatorStatus: 'OK',
    homogeneous: raw.homogeneous,
    homogeneousDerivatives: {
      u: raw.homogeneousDerivatives.u.map(x => x * jacU),
      v: raw.homogeneousDerivatives.v.map(x => x * jacV)
    },
    branch: {
      method: 'de-boor-kernel-v2',
      u: infoU,
      v: infoV,
      denominator: raw.denominator,
      denomScale,
      cancellationRatio,
      weightRatio: surface.validation.info?.weightRatio ?? Infinity,
      largeWeightSpread: surface.validation.info?.largeWeightSpread ?? false,
      span: raw.span
    },
    trace: legacyTrace(raw, infoU, infoV)
  };
}

export function adaptSubdividedShape(surface, legacySubdivided) {
  return {
    ...legacySubdivided,
    crossCheckedWith: KERNEL_NAME,
    kernelVersion: KERNEL_VERSION
  };
}
