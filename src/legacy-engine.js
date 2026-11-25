import { basisValues, basisDerivatives, knotAxisInfo, findDomainSpan } from './basis.js';
import { frameFromDerivatives } from './math-utils.js';

function legacyEvaluate(surface, u, v, options = {}) {
  const infoU = knotAxisInfo(surface.U, surface.p, surface.nu, u, options.sideU ?? 'auto');
  const infoV = knotAxisInfo(surface.V, surface.q, surface.nv, v, options.sideV ?? 'auto');
  const trace = [];
  if (!infoU.inDomain || !infoV.inDomain) {
    return {
      ok: false,
      code: 'OUT_OF_DOMAIN',
      message: `参数 (${u}, ${v}) 超出定义域 U${surface.domainU.join(',')} × V${surface.domainV.join(',')}。`,
      u, v, infoU, infoV, trace
    };
  }
  trace.push(...infoU.trace, ...infoV.trace);

  const bu = basisValues(surface.U, surface.p, u, infoU.side);
  const bv = basisValues(surface.V, surface.q, v, infoV.side);
  const du = basisDerivatives(surface.U, surface.p, u, infoU.side);
  const dv = basisDerivatives(surface.V, surface.q, v, infoV.side);

  const H = [0, 0, 0, 0];
  const Hu = [0, 0, 0, 0];
  const Hv = [0, 0, 0, 0];
  let denomScale = 0;
  let denomUScale = 0;
  let denomVScale = 0;
  let xyzScale = 0;
  const activeZeroWeights = [];
  const activeIndicesU = [];
  const activeIndicesV = [];

  for (let i = 0; i < surface.nu; i++) {
    if (Math.abs(bu[i]) > 1e-14) activeIndicesU.push(i);
    for (let j = 0; j < surface.nv; j++) {
      if (i === 0 && Math.abs(bv[j]) > 1e-14) activeIndicesV.push(j);
      const a = bu[i] * bv[j];
      const au = du[i] * bv[j];
      const av = bu[i] * dv[j];
      const w = surface.W[i][j];
      const x = surface.P[i][j][0] * w;
      const y = surface.P[i][j][1] * w;
      const z = surface.P[i][j][2] * w;
      const h = [x, y, z, w];
      for (let k = 0; k < 4; k++) {
        H[k] += a * h[k];
        Hu[k] += au * h[k];
        Hv[k] += av * h[k];
      }
      denomScale += Math.abs(a * w);
      denomUScale += Math.abs(au * w);
      denomVScale += Math.abs(av * w);
      xyzScale += Math.abs(a) * (Math.abs(x) + Math.abs(y) + Math.abs(z));
      if (a !== 0 && w === 0) activeZeroWeights.push([i, j]);
    }
  }

  const denominator = H[3];
  const cancellationRatio = denomScale > 0 ? Math.abs(denominator) / denomScale : 0;
  const denominatorFinite = Number.isFinite(denominator) && H.every(Number.isFinite) && Hu.every(Number.isFinite) && Hv.every(Number.isFinite);

  const branch = {
    method: 'homogeneous',
    u: infoU,
    v: infoV,
    activeIndicesU,
    activeIndicesV,
    activeZeroWeights,
    denominator,
    denomScale,
    cancellationRatio,
    weightRatio: surface.validation.info?.weightRatio ?? Infinity,
    largeWeightSpread: surface.validation.info?.largeWeightSpread ?? false
  };

  if (!denominatorFinite || denominator === 0 || denomScale === 0) {
    trace.push({ code: 'DEGENERATE_DENOMINATOR', message: '齐次有理分母为零、非有限或无有效权重贡献；不进行除法，也不生成法线。' });
    const limitCandidates = legacyLimitCandidates(H, Hu, Hv, denomUScale, denomVScale, xyzScale, trace);
    return {
      ok: false,
      code: 'DEGENERATE_DENOMINATOR',
      message: '有理分母退化；已保留端点/单侧与齐次导数分支。',
      engine: 'legacy-basis-sum',
      u, v,
      sideU: infoU.side,
      sideV: infoV.side,
      homogeneous: H,
      homogeneousDerivatives: { u: Hu, v: Hv },
      branch,
      limitCandidates,
      trace
    };
  }

  const invW = 1 / denominator;
  const position = [H[0] * invW, H[1] * invW, H[2] * invW];
  const derivativeU = [
    (Hu[0] - position[0] * Hu[3]) * invW,
    (Hu[1] - position[1] * Hu[3]) * invW,
    (Hu[2] - position[2] * Hu[3]) * invW
  ];
  const derivativeV = [
    (Hv[0] - position[0] * Hv[3]) * invW,
    (Hv[1] - position[1] * Hv[3]) * invW,
    (Hv[2] - position[2] * Hv[3]) * invW
  ];

  let denominatorStatus = 'OK';
  if (cancellationRatio < 1e-12) {
    denominatorStatus = 'ILL_CONDITIONED';
    trace.push({ code: 'ILL_CONDITIONED_DENOMINATOR', message: `分母抵消比 |W|/Σ|Nw|=${cancellationRatio.toExponential(3)} 过小；必须比较局部细分结果。` });
  }
  if (branch.largeWeightSpread) {
    trace.push({ code: 'LARGE_WEIGHT_SPREAD_COMPARE_SUBDIVISION', message: '权重差异超过 1e6；已具备齐次求值与节点插入/Bézier 局部细分对照分支。' });
  }
  if (activeZeroWeights.length) {
    trace.push({ code: 'ACTIVE_ZERO_WEIGHT', message: `当前基函数贡献内存在零权重控制点 ${JSON.stringify(activeZeroWeights)}；齐次贡献为零。` });
  }

  const frame = frameFromDerivatives(derivativeU, derivativeV, surface.scale);
  trace.push(...frame.trace.map(message => ({ code: frame.status, message })));

  return {
    ok: true,
    engine: 'legacy-basis-sum',
    u, v,
    sideU: infoU.side,
    sideV: infoV.side,
    position,
    derivatives: { u: derivativeU, v: derivativeV },
    normal: frame.normal,
    normalStatus: frame.status,
    frame,
    denominatorStatus,
    homogeneous: H,
    homogeneousDerivatives: { u: Hu, v: Hv },
    branch,
    trace
  };
}

function legacyLimitCandidates(H, Hu, Hv, denomUScale, denomVScale, xyzScale, trace) {
  const candidates = {};
  const hxyzNorm = Math.hypot(H[0], H[1], H[2]);
  const homogeneousPointIsZero = hxyzNorm <= 1e-11 * Math.max(xyzScale, 1e-300);

  if (Number.isFinite(Hu[3]) && Math.abs(Hu[3]) > 1e-12 * Math.max(denomUScale, 1e-300)) {
    candidates.alongU = {
      position: [Hu[0] / Hu[3], Hu[1] / Hu[3], Hu[2] / Hu[3]],
      rule: 'H_u^xyz / H_u^w',
      unique: homogeneousPointIsZero
    };
    trace.push({ code: 'ENDPOINT_LIMIT_U_BRANCH', message: '分母为零时保留沿 U 向内的一阶齐次极限候选 H_u^xyz/H_u^w。' });
  }
  if (Number.isFinite(Hv[3]) && Math.abs(Hv[3]) > 1e-12 * Math.max(denomVScale, 1e-300)) {
    candidates.alongV = {
      position: [Hv[0] / Hv[3], Hv[1] / Hv[3], Hv[2] / Hv[3]],
      rule: 'H_v^xyz / H_v^w',
      unique: homogeneousPointIsZero
    };
    trace.push({ code: 'ENDPOINT_LIMIT_V_BRANCH', message: '分母为零时保留沿 V 向内的一阶齐次极限候选 H_v^xyz/H_v^w；角点处不擅自合并为唯一值。' });
  }
  return candidates;
}

export { legacyEvaluate };
