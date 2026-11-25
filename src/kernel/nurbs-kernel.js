export const KERNEL_NAME = 'nurbs-kernel-v2';
export const KERNEL_VERSION = '2.0.0';
export const NORMALIZED_DOMAIN = Object.freeze([0, 1]);

const SIDE_CODE = { L: -1, left: -1, R: 1, right: 1 };

export function normalizeSide(side) {
  if (side === -1 || side === 'L' || side === 'left') return 'L';
  return 'R';
}

export function domainOf(knots, degree, controlCount = knots.length - degree - 1) {
  return [knots[degree], knots[controlCount]];
}

export function normalizeParameter(t, knots, degree) {
  const [lo, hi] = domainOf(knots, degree);
  if (hi === lo) return 0;
  return (t - lo) / (hi - lo);
}

export function knotSpanIndex(knots, degree, t, side = 'R') {
  const n = knots.length - degree - 2;
  if (t >= knots[n + 1]) return n;
  const code = SIDE_CODE[side] ?? 1;
  if (code < 0 && t > knots[degree]) {
    for (let k = n; k >= degree; k--) {
      if (knots[k] < t && t <= knots[k + 1]) return k;
    }
  }
  for (let k = degree; k <= n; k++) {
    if (knots[k] <= t && t < knots[k + 1]) return k;
  }
  return degree;
}

function sideSign(side) {
  return SIDE_CODE[side] ?? 1;
}

function basisMatrix(knots, degree, t, side) {
  const sign = sideSign(side);
  const count = knots.length - degree - 1;
  const N = Array.from({ length: count + 1 }, () => new Array(degree + 1).fill(0));
  for (let i = 0; i < count; i++) {
    if (knots[i] === knots[i + 1]) {
      N[i][0] = 0;
      continue;
    }
    N[i][0] = sign < 0
      ? (knots[i] < t && t <= knots[i + 1] ? 1 : 0)
      : (knots[i] <= t && t < knots[i + 1] ? 1 : 0);
  }
  for (let q = 1; q <= degree; q++) {
    for (let i = 0; i < count; i++) {
      const d1 = knots[i + q] - knots[i];
      const d2 = knots[i + q + 1] - knots[i + 1];
      let value = 0;
      if (d1 !== 0) value += ((t - knots[i]) / d1) * N[i][q - 1];
      if (d2 !== 0) value += ((knots[i + q + 1] - t) / d2) * N[i + 1][q - 1];
      N[i][q] = value;
    }
  }
  return N;
}

export function kernelBasisValues(knots, degree, t, side = 'R') {
  const s = normalizeSide(side);
  const matrix = basisMatrix(knots, degree, t, s);
  return matrix.map(row => row[degree]);
}

export function kernelBasisDerivatives(knots, degree, t, side = 'R') {
  const s = normalizeSide(side);
  const count = knots.length - degree - 1;
  if (degree === 0) return new Array(count).fill(0);
  const lower = basisMatrix(knots, degree - 1, t, s).map(row => row[degree - 1]);
  const result = new Array(count).fill(0);
  for (let i = 0; i < count; i++) {
    const d1 = knots[i + degree] - knots[i];
    if (d1 !== 0) result[i] += (degree * lower[i]) / d1;
    const d2 = knots[i + degree + 1] - knots[i + 1];
    if (d2 !== 0) result[i] -= (degree * (lower[i + 1] ?? 0)) / d2;
  }
  return result;
}

function weightedSum(values, points, dim) {
  const out = new Array(dim).fill(0);
  for (let i = 0; i < values.length; i++) {
    const a = values[i];
    if (a === 0) continue;
    const p = points[i];
    for (let d = 0; d < dim; d++) out[d] += a * p[d];
  }
  return out;
}

export class NurbsCurveKernel {
  constructor({ degree, knots, homogeneousPoints }) {
    if (!Number.isInteger(degree) || degree < 0) throw new Error('KERNEL_INVALID_DEGREE');
    if (!Array.isArray(knots) || knots.length !== homogeneousPoints.length + degree + 1) {
      throw new Error('KERNEL_KNOT_LENGTH_MISMATCH');
    }
    this.degree = degree;
    this.knots = knots.slice();
    this.hpoints = homogeneousPoints.map(p => p.slice());
    this.dim = this.hpoints[0].length;
    this.domain = domainOf(this.knots, this.degree);
  }

  inDomain(t) {
    const [lo, hi] = this.domain;
    return t >= lo && t <= hi;
  }

  evaluateHomogeneous(t, side = 'R') {
    if (!this.inDomain(t)) {
      return { ok: false, code: 'OUT_OF_DOMAIN', domain: this.domain, t };
    }
    const s = normalizeSide(side);
    const values = kernelBasisValues(this.knots, this.degree, t, s);
    const h = weightedSum(values, this.hpoints, this.dim);
    return { ok: true, h, span: knotSpanIndex(this.knots, this.degree, t, s), side: s };
  }

  evaluateDerivativeHomogeneous(t, side = 'R') {
    if (!this.inDomain(t)) {
      return { ok: false, code: 'OUT_OF_DOMAIN', domain: this.domain, t };
    }
    const s = normalizeSide(side);
    const derivs = kernelBasisDerivatives(this.knots, this.degree, t, s);
    const dh = weightedSum(derivs, this.hpoints, this.dim);
    return { ok: true, dh, side: s };
  }

  evaluate(t, side = 'R') {
    const point = this.evaluateHomogeneous(t, side);
    if (!point.ok) return point;
    const h = point.h;
    const w = h[h.length - 1];
    if (!Number.isFinite(w) || w === 0) {
      return {
        ok: false,
        code: 'DEGENERATE_DENOMINATOR',
        homogeneous: h,
        side: point.side,
        t,
        message: '新内核：齐次权重为零，拒绝除法并返回退化状态。'
      };
    }
    const position = h.slice(0, -1).map(x => x / w);
    const deriv = this.evaluateDerivativeHomogeneous(t, side);
    if (!deriv.ok) return { ok: false, ...deriv };
    const dh = deriv.dh;
    const dw = dh[dh.length - 1];
    const derivative = dh.slice(0, -1).map((x, i) => (x - position[i] * dw) / w);
    return {
      ok: true,
      position,
      derivative,
      homogeneous: h,
      homogeneousDerivative: dh,
      denominator: w,
      side: point.side,
      t,
      engine: KERNEL_NAME
    };
  }

  invertPoint(target, hint = null, options = {}) {
    const [lo, hi] = this.domain;
    let t = hint === null ? (lo + hi) / 2 : Math.min(hi, Math.max(lo, hint));
    const iterations = options.iterations ?? 48;
    const tolerance = options.tolerance ?? 1e-11;
    const spatialScale = options.scale ?? 1;
    let residual = null;
    for (let n = 0; n < iterations; n++) {
      const r = this.evaluate(t);
      if (!r.ok) return { ok: false, code: r.code, t };
      const f = r.position.map((x, i) => x - target[i]);
      residual = Math.hypot(...f);
      if (residual <= tolerance * Math.max(1, spatialScale)) {
        return { ok: true, parameter: t, residual, iterations: n + 1 };
      }
      let jtj = 0;
      let jtr = 0;
      for (let i = 0; i < target.length; i++) {
        jtj += r.derivative[i] * r.derivative[i];
        jtr += r.derivative[i] * f[i];
      }
      if (jtj <= 1e-300) break;
      const step = jtr / jtj;
      let next = t - step;
      let accepted = false;
      for (let back = 0; back < 24; back++) {
        next = Math.min(hi, Math.max(lo, next));
        const probe = this.evaluate(next);
        if (probe.ok) {
          const probeResidual = Math.hypot(...probe.position.map((x, i) => x - target[i]));
          if (probeResidual <= residual) {
            accepted = true;
            break;
          }
        }
        next = t - step * Math.pow(0.5, back + 1);
      }
      if (!accepted) break;
      t = next;
    }
    return { ok: false, code: 'INVERSION_DID_NOT_CONVERGE', parameter: t, residual, iterations };
  }
}

function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

export class NurbsSurfaceKernel {
  constructor({ degreeU, degreeV, knotsU, knotsV, homogeneousGrid }) {
    this.p = degreeU;
    this.q = degreeV;
    this.U = knotsU.slice();
    this.V = knotsV.slice();
    this.grid = homogeneousGrid.map(row => row.map(p => p.slice()));
    this.nu = this.grid.length;
    this.nv = this.grid[0].length;
    this.dim = this.grid[0][0].length;
    this.domainU = domainOf(this.U, this.p);
    this.domainV = domainOf(this.V, this.q);
  }

  #uCurves() {
    return Array.from({ length: this.nv }, (_, j) =>
      new NurbsCurveKernel({
        degree: this.p,
        knots: this.U,
        homogeneousPoints: this.grid.map(row => row[j])
      })
    );
  }

  evaluateHomogeneous(u, v, sideU = 'R', sideV = 'R') {
    const su = normalizeSide(sideU);
    const sv = normalizeSide(sideV);
    if (u < this.domainU[0] || u > this.domainU[1] || v < this.domainV[0] || v > this.domainV[1]) {
      return { ok: false, code: 'OUT_OF_DOMAIN', domainU: this.domainU, domainV: this.domainV, u, v };
    }
    const iso = this.#uCurves().map(curve => curve.evaluateHomogeneous(u, su).h);
    const vCurve = new NurbsCurveKernel({ degree: this.q, knots: this.V, homogeneousPoints: iso });
    const result = vCurve.evaluateHomogeneous(v, sv);
    if (!result.ok) return result;
    return {
      ok: true,
      h: result.h,
      sideU: su,
      sideV: sv,
      spanU: knotSpanIndex(this.U, this.p, u, su),
      spanV: knotSpanIndex(this.V, this.q, v, sv)
    };
  }

  evaluate(u, v, options = {}) {
    const sideU = normalizeSide(options.sideU ?? 'R');
    const sideV = normalizeSide(options.sideV ?? 'R');
    if (u < this.domainU[0] || u > this.domainU[1] || v < this.domainV[0] || v > this.domainV[1]) {
      return {
        ok: false,
        code: 'OUT_OF_DOMAIN',
        engine: KERNEL_NAME,
        domainU: this.domainU,
        domainV: this.domainV,
        u, v
      };
    }

    const uCurves = this.#uCurves();
    const isoH = [];
    const isoHu = [];
    for (const curve of uCurves) {
      const h = curve.evaluateHomogeneous(u, sideU);
      const dh = curve.evaluateDerivativeHomogeneous(u, sideU);
      if (!h.ok || !dh.ok) return { ok: false, code: h.code ?? dh.code, engine: KERNEL_NAME, u, v };
      isoH.push(h.h);
      isoHu.push(dh.dh);
    }
    const vCurveH = new NurbsCurveKernel({ degree: this.q, knots: this.V, homogeneousPoints: isoH });
    const vCurveHu = new NurbsCurveKernel({ degree: this.q, knots: this.V, homogeneousPoints: isoHu });
    const Hpoint = vCurveH.evaluateHomogeneous(v, sideV);
    const HuPoint = vCurveHu.evaluateHomogeneous(v, sideV);
    const HvPoint = vCurveH.evaluateDerivativeHomogeneous(v, sideV);
    if (!Hpoint.ok || !HuPoint.ok || !HvPoint.ok) {
      return { ok: false, code: 'KERNEL_EVALUATION_FAILED', engine: KERNEL_NAME, u, v };
    }
    const H = Hpoint.h;
    const Hu = HuPoint.h;
    const Hv = HvPoint.dh;
    const w = H[this.dim - 1];

    if (![...H, ...Hu, ...Hv].every(Number.isFinite) || w === 0) {
      const limitCandidates = {};
      if (Number.isFinite(Hu[this.dim - 1]) && Hu[this.dim - 1] !== 0) {
        limitCandidates.alongU = { position: Hu.slice(0, -1).map(x => x / Hu[this.dim - 1]), rule: 'H_u^xyz / H_u^w' };
      }
      if (Number.isFinite(Hv[this.dim - 1]) && Hv[this.dim - 1] !== 0) {
        limitCandidates.alongV = { position: Hv.slice(0, -1).map(x => x / Hv[this.dim - 1]), rule: 'H_v^xyz / H_v^w' };
      }
      return {
        ok: false,
        code: 'DEGENERATE_DENOMINATOR',
        engine: KERNEL_NAME,
        kernelVersion: KERNEL_VERSION,
        u, v,
        sideU, sideV,
        homogeneous: H,
        homogeneousDerivatives: { u: Hu, v: Hv },
        limitCandidates,
        message: '新内核：有理分母退化，返回明确状态与齐次极限候选。'
      };
    }

    const invW = 1 / w;
    const position = H.slice(0, -1).map(x => x * invW);
    const derivativeU = new Array(this.dim - 1);
    const derivativeV = new Array(this.dim - 1);
    for (let i = 0; i < this.dim - 1; i++) {
      derivativeU[i] = (Hu[i] - position[i] * Hu[this.dim - 1]) * invW;
      derivativeV[i] = (Hv[i] - position[i] * Hv[this.dim - 1]) * invW;
    }
    return {
      ok: true,
      engine: KERNEL_NAME,
      kernelVersion: KERNEL_VERSION,
      u, v,
      sideU, sideV,
      position,
      derivatives: { u: derivativeU, v: derivativeV },
      homogeneous: H,
      homogeneousDerivatives: { u: Hu, v: Hv },
      denominator: w,
      span: {
        u: knotSpanIndex(this.U, this.p, u, sideU),
        v: knotSpanIndex(this.V, this.q, v, sideV)
      }
    };
  }

  invertPoint(target, hintU = null, hintV = null, options = {}) {
    const candidates = [];
    const midU = (this.domainU[0] + this.domainU[1]) / 2;
    const midV = (this.domainV[0] + this.domainV[1]) / 2;
    if (hintU !== null && hintV !== null) candidates.push([hintU, hintV]);
    candidates.push([midU, midV]);
    for (const fu of [0.25, 0.75]) {
      for (const fv of [0.25, 0.75]) {
        candidates.push([
          this.domainU[0] + fu * (this.domainU[1] - this.domainU[0]),
          this.domainV[0] + fv * (this.domainV[1] - this.domainV[0])
        ]);
      }
    }

    let lastFailure = null;
    for (const [startU, startV] of candidates) {
      const attempt = this.#invertFrom(target, startU, startV, options);
      if (attempt.ok) return attempt;
      lastFailure = attempt;
    }
    return lastFailure ?? { ok: false, code: 'INVERSION_DID_NOT_CONVERGE' };
  }

  #invertFrom(target, hintU, hintV, options) {
    let u = Math.min(this.domainU[1], Math.max(this.domainU[0], hintU));
    let v = Math.min(this.domainV[1], Math.max(this.domainV[0], hintV));
    const iterations = options.iterations ?? 64;
    const tolerance = options.tolerance ?? 1e-10;
    const spatialScale = options.scale ?? 1;
    const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));
    let residual = null;
    for (let n = 0; n < iterations; n++) {
      const r = this.evaluate(u, v, options);
      if (!r.ok) return { ok: false, code: r.code, u, v };
      const f = r.position.map((x, i) => x - target[i]);
      residual = Math.hypot(...f);
      if (residual <= tolerance * Math.max(1, spatialScale)) {
        return { ok: true, u, v, residual, iterations: n + 1 };
      }
      const J = [r.derivatives.u, r.derivatives.v];
      const jtj = [
        [dot(J[0], J[0]), dot(J[0], J[1])],
        [dot(J[0], J[1]), dot(J[1], J[1])]
      ];
      const jtr = [dot(J[0], f), dot(J[1], f)];
      const det = jtj[0][0] * jtj[1][1] - jtj[0][1] * jtj[0][1];
      if (Math.abs(det) <= 1e-300) break;
      const du = (jtj[1][1] * jtr[0] - jtj[0][1] * jtr[1]) / det;
      const dv = (jtj[0][0] * jtr[1] - jtj[0][1] * jtr[0]) / det;
      let accepted = false;
      let nextU = u;
      let nextV = v;
      for (let back = 0; back < 24; back++) {
        const shrink = Math.pow(0.5, back);
        nextU = clamp(u - du * shrink, this.domainU[0], this.domainU[1]);
        nextV = clamp(v - dv * shrink, this.domainV[0], this.domainV[1]);
        const probe = this.evaluate(nextU, nextV, options);
        if (probe.ok) {
          const probeResidual = Math.hypot(...probe.position.map((x, i) => x - target[i]));
          if (probeResidual <= residual) {
            accepted = true;
            break;
          }
        }
      }
      if (!accepted) break;
      u = nextU;
      v = nextV;
    }
    return { ok: false, code: 'INVERSION_DID_NOT_CONVERGE', u, v, residual, iterations };
  }
}

export function reparametrizeKnots(knots, degree, targetDomain = NORMALIZED_DOMAIN) {
  const [lo, hi] = domainOf(knots, degree);
  const [t0, t1] = targetDomain;
  if (hi === lo) throw new Error('KERNEL_EMPTY_DOMAIN');
  const scale = (t1 - t0) / (hi - lo);
  return {
    knots: knots.map(k => t0 + (k - lo) * scale),
    map: t => t0 + (t - lo) * scale,
    inverse: s => lo + (s - t0) / scale,
    jacobian: scale
  };
}
