import { validateSpec, normalizeSpec } from './validation.js';
import { basisValues, basisDerivatives, knotAxisInfo, findDomainSpan } from './basis.js';
import { frameFromDerivatives } from './math-utils.js';
import { evaluateSubdivided } from './subdivision.js';

export function compileSurface(input) {
  const validation = validateSpec(input);
  if (!validation.valid) {
    const error = new Error(validation.errors.map(e => `${e.code}: ${e.message}`).join('\n'));
    error.validation = validation;
    throw error;
  }
  return new RationalSurface(normalizeSpec(input));
}

export class RationalSurface {
  constructor(spec) {
    this.spec = spec;
    this.id = spec.id;
    this.name = spec.name;
    this.p = spec.degreeU;
    this.q = spec.degreeV;
    this.U = spec.knotsU;
    this.V = spec.knotsV;
    this.P = spec.points;
    this.W = spec.weights;
    this.nu = this.P.length;
    this.nv = this.P[0].length;
    this.domainU = [this.U[this.p], this.U[this.nu]];
    this.domainV = [this.V[this.q], this.V[this.nv]];
    this.scale = spec.info?.scale ?? 1;
    this.validation = validateSpec(spec);
  }

  evaluate(u, v, options = {}) {
    const infoU = knotAxisInfo(this.U, this.p, this.nu, u, options.sideU ?? 'auto');
    const infoV = knotAxisInfo(this.V, this.q, this.nv, v, options.sideV ?? 'auto');
    const trace = [];
    if (!infoU.inDomain || !infoV.inDomain) {
      return {
        ok: false,
        code: 'OUT_OF_DOMAIN',
        message: `参数 (${u}, ${v}) 超出定义域 U${this.domainU.join(',')} × V${this.domainV.join(',')}。`,
        u, v, infoU, infoV, trace
      };
    }
    trace.push(...infoU.trace, ...infoV.trace);

    const bu = basisValues(this.U, this.p, u, infoU.side);
    const bv = basisValues(this.V, this.q, v, infoV.side);
    const du = basisDerivatives(this.U, this.p, u, infoU.side);
    const dv = basisDerivatives(this.V, this.q, v, infoV.side);

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

    for (let i = 0; i < this.nu; i++) {
      if (Math.abs(bu[i]) > 1e-14) activeIndicesU.push(i);
      for (let j = 0; j < this.nv; j++) {
        if (i === 0 && Math.abs(bv[j]) > 1e-14) activeIndicesV.push(j);
        const a = bu[i] * bv[j];
        const au = du[i] * bv[j];
        const av = bu[i] * dv[j];
        const w = this.W[i][j];
        const x = this.P[i][j][0] * w;
        const y = this.P[i][j][1] * w;
        const z = this.P[i][j][2] * w;
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
      weightRatio: this.validation.info?.weightRatio ?? Infinity,
      largeWeightSpread: this.validation.info?.largeWeightSpread ?? false
    };

    if (!denominatorFinite || denominator === 0 || denomScale === 0) {
      trace.push({ code: 'DEGENERATE_DENOMINATOR', message: '齐次有理分母为零、非有限或无有效权重贡献；不进行除法，也不生成法线。' });
      const limitCandidates = this.#limitCandidates(H, Hu, Hv, denomUScale, denomVScale, xyzScale, trace);
      return {
        ok: false,
        code: 'DEGENERATE_DENOMINATOR',
        message: '有理分母退化；已保留端点/单侧与齐次导数分支。',
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

    const frame = frameFromDerivatives(derivativeU, derivativeV, this.scale);
    trace.push(...frame.trace.map(message => ({ code: frame.status, message })));

    return {
      ok: true,
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

  #limitCandidates(H, Hu, Hv, denomUScale, denomVScale, xyzScale, trace) {
    const candidates = {};
    const hxyzNorm = Math.hypot(H[0], H[1], H[2]);
    const homogeneousPointIsZero = hxyzNorm <= 1e-11 * Math.max(xyzScale, 1e-300);

    if (Number.isFinite(Hu[3]) && Math.abs(Hu[3]) > 1e-12 * Math.max(denomUScale, 1e-300)) {
      candidates.alongU = {
        position: [Hu[0] / Hu[3], Hu[1] / Hu[3], Hu[2] / Hu[3]],
        rule: "H_u^xyz / H_u^w",
        unique: homogeneousPointIsZero
      };
      trace.push({ code: 'ENDPOINT_LIMIT_U_BRANCH', message: '分母为零时保留沿 U 向内的一阶齐次极限候选 H_u^xyz/H_u^w。' });
    }
    if (Number.isFinite(Hv[3]) && Math.abs(Hv[3]) > 1e-12 * Math.max(denomVScale, 1e-300)) {
      candidates.alongV = {
        position: [Hv[0] / Hv[3], Hv[1] / Hv[3], Hv[2] / Hv[3]],
        rule: "H_v^xyz / H_v^w",
        unique: homogeneousPointIsZero
      };
      trace.push({ code: 'ENDPOINT_LIMIT_V_BRANCH', message: '分母为零时保留沿 V 向内的一阶齐次极限候选 H_v^xyz/H_v^w；角点处不擅自合并为唯一值。' });
    }
    return candidates;
  }

  evaluateSubdivided(u, v, options = {}) {
    return evaluateSubdivided(this, u, v, options);
  }

  spanAt(u, v, options = {}) {
    const iu = knotAxisInfo(this.U, this.p, this.nu, u, options.sideU ?? 'auto');
    const iv = knotAxisInfo(this.V, this.q, this.nv, v, options.sideV ?? 'auto');
    if (!iu.inDomain || !iv.inDomain) return null;
    const ku = findDomainSpan(this.U, this.p, this.nu, u, iu.side);
    const kv = findDomainSpan(this.V, this.q, this.nv, v, iv.side);
    return { ku, kv, infoU: iu, infoV: iv };
  }

  tessellate(uSegments = 64, vSegments = 32) {
    uSegments = Math.max(1, Math.floor(uSegments));
    vSegments = Math.max(1, Math.floor(vSegments));
    const [u0, u1] = this.domainU;
    const [v0, v1] = this.domainV;
    const cols = uSegments + 1;
    const rows = vSegments + 1;
    const positions = [];
    const normals = [];
    const uvs = [];
    const valid = new Array(cols * rows).fill(false);
    const evaluations = new Array(cols * rows).fill(null);
    const failures = [];

    for (let j = 0; j <= vSegments; j++) {
      const v = v0 + ((v1 - v0) * j) / vSegments;
      for (let i = 0; i <= uSegments; i++) {
        const u = u0 + ((u1 - u0) * i) / uSegments;
        const idx = j * cols + i;
        const r = this.evaluate(u, v);
        evaluations[idx] = r;
        uvs.push(i / uSegments, j / vSegments);
        if (r.ok && r.normal) {
          valid[idx] = true;
          positions.push(...r.position);
          normals.push(...r.normal);
        } else {
          positions.push(NaN, NaN, NaN);
          normals.push(0, 0, 0);
          failures.push({ i, j, u, v, code: r.code ?? r.normalStatus, message: r.message });
        }
      }
    }

    const indices = [];
    for (let j = 0; j < vSegments; j++) {
      for (let i = 0; i < uSegments; i++) {
        const a = j * cols + i;
        const b = a + 1;
        const c = a + cols;
        const d = c + 1;
        if (valid[a] && valid[b] && valid[c]) indices.push(a, b, c);
        if (valid[b] && valid[d] && valid[c]) indices.push(b, d, c);
      }
    }

    return { positions, normals, uvs, indices, valid, evaluations, failures, cols, rows };
  }

  controlNet() {
    const lineIndices = [];
    for (let i = 0; i < this.nu; i++) {
      for (let j = 0; j < this.nv; j++) {
        const idx = i * this.nv + j;
        if (i + 1 < this.nu) lineIndices.push(idx, (i + 1) * this.nv + j);
        if (j + 1 < this.nv) lineIndices.push(idx, i * this.nv + j + 1);
      }
    }
    const points = this.P.flatMap((row, i) => row.flatMap((p, j) => [...p, this.W[i][j]]));
    return { points, lineIndices };
  }
}
