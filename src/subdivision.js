import { bernstein, bernsteinDerivative, deCasteljau4 } from './bernstein.js';
import { frameFromDerivatives } from './math-utils.js';

function homogeneousGrid(surface) {
  return surface.P.map((row, i) =>
    row.map((p, j) => [p[0] * surface.W[i][j], p[1] * surface.W[i][j], p[2] * surface.W[i][j], surface.W[i][j]])
  );
}

function mix4(a, b, alpha) {
  return [
    alpha * a[0] + (1 - alpha) * b[0],
    alpha * a[1] + (1 - alpha) * b[1],
    alpha * a[2] + (1 - alpha) * b[2],
    alpha * a[3] + (1 - alpha) * b[3]
  ];
}

function insertionSpan(knots, t) {
  for (let k = 0; k < knots.length - 1; k++) {
    if (knots[k] <= t && t < knots[k + 1]) return k;
  }
  return -1;
}

export function insertKnotU(grid, knots, degree, t) {
  const oldCount = grid.length;
  const k = insertionSpan(knots, t);
  if (k < 0) throw new Error(`U 向节点插入失败：t=${t} 不在可插入区间。`);
  const nextKnots = knots.slice();
  nextKnots.splice(k + 1, 0, t);

  const next = new Array(oldCount + 1);
  const fixedEnd = Math.min(k - degree, oldCount - 1);
  for (let i = 0; i <= fixedEnd; i++) next[i] = grid[i].map(row => row.slice());

  for (let i = Math.max(0, k - degree + 1); i <= k; i++) {
    const d = knots[i + degree] - knots[i];
    const alpha = d === 0 ? 0 : (t - knots[i]) / d;
    next[i] = grid[0].map((_, j) => {
      const current = i < oldCount ? grid[i][j] : null;
      const previous = grid[i - 1][j];
      if (!current) return previous.slice();
      return mix4(current, previous, alpha);
    });
  }

  for (let i = k + 1; i <= oldCount; i++) {
    next[i] = grid[i - 1].map(row => row.slice());
  }
  return { grid: next, knots: nextKnots, span: k };
}

export function insertKnotV(grid, knots, degree, t) {
  const oldRows = grid.length;
  const oldCount = grid[0].length;
  const k = insertionSpan(knots, t);
  if (k < 0) throw new Error(`V 向节点插入失败：t=${t} 不在可插入区间。`);
  const nextKnots = knots.slice();
  nextKnots.splice(k + 1, 0, t);

  const next = Array.from({ length: oldRows }, () => new Array(oldCount + 1));
  for (let i = 0; i < oldRows; i++) {
    for (let j = 0; j <= Math.min(k - degree, oldCount - 1); j++) next[i][j] = grid[i][j].slice();

    for (let j = Math.max(0, k - degree + 1); j <= k; j++) {
      const d = knots[j + degree] - knots[j];
      const alpha = d === 0 ? 0 : (t - knots[j]) / d;
      const previous = grid[i][j - 1];
      const current = j < oldCount ? grid[i][j] : null;
      next[i][j] = current ? mix4(current, previous, alpha) : previous.slice();
    }

    for (let j = k + 1; j <= oldCount; j++) next[i][j] = grid[i][j - 1].slice();
  }
  return { grid: next, knots: nextKnots, span: k };
}

function multiplicity(knots, t) {
  return knots.reduce((n, k) => n + (k === t ? 1 : 0), 0);
}

function raiseAxis(grid, knots, degree, value, axis) {
  let g = grid;
  let ks = knots;
  const initial = multiplicity(ks, value);
  const inserts = Math.max(0, degree - initial);
  for (let n = 0; n < inserts; n++) {
    const fn = axis === 'u' ? insertKnotU : insertKnotV;
    const result = fn(g, ks, degree, value);
    g = result.grid;
    ks = result.knots;
  }
  return { grid: g, knots: ks, initialMultiplicity: initial, insertions: inserts };
}

function findTransitionSpan(knots, a, b) {
  for (let k = 0; k < knots.length - 1; k++) {
    if (knots[k] === a && knots[k + 1] === b) return k;
  }
  return -1;
}

function sumPatch(patch, bu, bv, dbu, dbv, lenU, lenV) {
  const H = [0, 0, 0, 0];
  const Hu = [0, 0, 0, 0];
  const Hv = [0, 0, 0, 0];
  for (let i = 0; i < bu.length; i++) {
    for (let j = 0; j < bv.length; j++) {
      const h = patch[i][j];
      const a = bu[i] * bv[j];
      const au = (dbu[i] / lenU) * bv[j];
      const av = bu[i] * (dbv[j] / lenV);
      for (let k = 0; k < 4; k++) {
        H[k] += a * h[k];
        Hu[k] += au * h[k];
        Hv[k] += av * h[k];
      }
    }
  }
  return { H, Hu, Hv };
}

export function evaluateSubdivided(surface, u, v, options = {}) {
  const located = surface.spanAt(u, v, options);
  if (!located || located.ku < 0 || located.kv < 0) {
    return { ok: false, code: 'OUT_OF_DOMAIN', method: 'local-subdivision' };
  }

  const { infoU, infoV } = located;
  const a = surface.U[located.ku];
  const b = surface.U[located.ku + 1];
  const c = surface.V[located.kv];
  const d = surface.V[located.kv + 1];
  if (!(b > a && d > c)) {
    return { ok: false, code: 'INVALID_EXTRACTION_SPAN', method: 'local-subdivision', a, b, c, d };
  }

  let grid = homogeneousGrid(surface);
  let U = surface.U.slice();
  let V = surface.V.slice();

  const raiseA = raiseAxis(grid, U, surface.p, a, 'u');
  grid = raiseA.grid; U = raiseA.knots;
  const raiseB = raiseAxis(grid, U, surface.p, b, 'u');
  grid = raiseB.grid; U = raiseB.knots;
  const raiseC = raiseAxis(grid, V, surface.q, c, 'v');
  grid = raiseC.grid; V = raiseC.knots;
  const raiseD = raiseAxis(grid, V, surface.q, d, 'v');
  grid = raiseD.grid; V = raiseD.knots;

  const ku = findTransitionSpan(U, a, b);
  const kv = findTransitionSpan(V, c, d);
  if (ku < surface.p || kv < surface.q) {
    return { ok: false, code: 'BEZIER_EXTRACTION_FAILED', method: 'local-subdivision', ku, kv };
  }

  const startU = ku - surface.p;
  const startV = kv - surface.q;
  const patch = Array.from({ length: surface.p + 1 }, (_, i) =>
    Array.from({ length: surface.q + 1 }, (_, j) => grid[startU + i][startV + j].slice())
  );

  const s = (u - a) / (b - a);
  const t = (v - c) / (d - c);
  const bu = bernstein(surface.p, s);
  const bv = bernstein(surface.q, t);
  const dbu = bernsteinDerivative(surface.p, s);
  const dbv = bernsteinDerivative(surface.q, t);
  const { H, Hu, Hv } = sumPatch(patch, bu, bv, dbu, dbv, b - a, d - c);
  const deCastH = deCasteljau4(patch, s, t);

  if (H[3] === 0 || !H.every(Number.isFinite) || deCastH[3] === 0) {
    return {
      ok: false,
      code: 'DEGENERATE_DENOMINATOR',
      method: 'local-subdivision',
      homogeneous: H,
      deCasteljauHomogeneous: deCastH,
      extraction: { a, b, c, d, s, t, startU, startV }
    };
  }

  const invW = 1 / H[3];
  const bernsteinPosition = [H[0] * invW, H[1] * invW, H[2] * invW];
  const deCastPosition = [deCastH[0] / deCastH[3], deCastH[1] / deCastH[3], deCastH[2] / deCastH[3]];
  const position = deCastPosition;
  const derivativeU = [
    (Hu[0] - bernsteinPosition[0] * Hu[3]) * invW,
    (Hu[1] - bernsteinPosition[1] * Hu[3]) * invW,
    (Hu[2] - bernsteinPosition[2] * Hu[3]) * invW
  ];
  const derivativeV = [
    (Hv[0] - bernsteinPosition[0] * Hv[3]) * invW,
    (Hv[1] - bernsteinPosition[1] * Hv[3]) * invW,
    (Hv[2] - bernsteinPosition[2] * Hv[3]) * invW
  ];
  const frame = frameFromDerivatives(derivativeU, derivativeV, surface.scale);

  return {
    ok: true,
    method: 'local-subdivision-de-casteljau',
    u, v,
    sideU: infoU.side,
    sideV: infoV.side,
    position,
    bernsteinPosition,
    deCasteljauPosition: deCastPosition,
    deCasteljauInternalError: Math.max(
      Math.abs(deCastPosition[0] - bernsteinPosition[0]),
      Math.abs(deCastPosition[1] - bernsteinPosition[1]),
      Math.abs(deCastPosition[2] - bernsteinPosition[2])
    ),
    derivatives: { u: derivativeU, v: derivativeV },
    normal: frame.normal,
    normalStatus: frame.status,
    frame,
    extraction: {
      a, b, c, d, s, t, startU, startV,
      insertions: {
        a: raiseA.insertions,
        b: raiseB.insertions,
        c: raiseC.insertions,
        d: raiseD.insertions
      }
    }
  };
}
