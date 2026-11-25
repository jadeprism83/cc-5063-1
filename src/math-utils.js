export function add3(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

export function sub3(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function scale3(a, s) {
  return [a[0] * s, a[1] * s, a[2] * s];
}

export function dot3(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function cross3(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}

export function norm3(a) {
  return Math.hypot(a[0], a[1], a[2]);
}

export function maxComponentError(a, b) {
  return Math.max(
    Math.abs(a[0] - b[0]),
    Math.abs(a[1] - b[1]),
    Math.abs(a[2] - b[2])
  );
}

export function vectorError(actual, expected, scale = 1) {
  const absolute = norm3(sub3(actual, expected));
  const denominator = Math.max(1e-14, norm3(expected), scale * 1e-12);
  return { absolute, relative: absolute / denominator };
}

export function frameFromDerivatives(du, dv, modelScale = 1) {
  const lu = norm3(du);
  const lv = norm3(dv);
  const tol = 1e-11 * Math.max(modelScale, 1e-12);
  const trace = [];

  if (![...du, ...dv].every(Number.isFinite)) {
    return { normal: null, status: 'NON_FINITE_DERIVATIVE', trace: ['导数含非有限分量，不生成法线。'] };
  }
  if (lu <= tol) {
    return { normal: null, status: 'DEGENERATE_U_DERIVATIVE', lengthU: lu, lengthV: lv, trace: [`|S_u|=${lu.toExponential(4)} <= ${tol.toExponential(4)}；U 向导数退化，不输出任意单位法线。`] };
  }
  if (lv <= tol) {
    return { normal: null, status: 'DEGENERATE_V_DERIVATIVE', lengthU: lu, lengthV: lv, trace: [`|S_v|=${lv.toExponential(4)} <= ${tol.toExponential(4)}；V 向导数退化，不输出任意单位法线。`] };
  }

  const cross = cross3(du, dv);
  const lc = norm3(cross);
  const areaScale = lu * lv;
  if (lc <= 1e-11 * areaScale || lc <= tol * tol) {
    return {
      normal: null,
      status: 'LINEARLY_DEPENDENT_TANGENTS',
      cross,
      crossLength: lc,
      lengthU: lu,
      lengthV: lv,
      trace: [`|S_u×S_v|/(|S_u||S_v|)=${(lc / areaScale).toExponential(4)}；切向量共线，法线不唯一。`]
    };
  }

  return {
    normal: scale3(cross, 1 / lc),
    status: 'OK',
    cross,
    crossLength: lc,
    lengthU: lu,
    lengthV: lv,
    trace: []
  };
}

export function percentile(sorted, p) {
  if (!sorted.length) return 0;
  const x = p * (sorted.length - 1);
  const i = Math.floor(x);
  const j = Math.min(i + 1, sorted.length - 1);
  return sorted[i] + (sorted[j] - sorted[i]) * (x - i);
}

export function distribution(values) {
  const finite = values.filter(Number.isFinite);
  if (!finite.length) {
    return { count: 0, min: null, max: null, mean: null, rms: null, p50: null, p95: null, histogram: [] };
  }
  const sorted = finite.slice().sort((a, b) => a - b);
  const min = sorted[0];
  const max = sorted[sorted.length - 1];
  const sum = finite.reduce((a, b) => a + b, 0);
  const squares = finite.reduce((a, b) => a + b * b, 0);
  const bins = new Array(12).fill(0);
  const span = max - min;
  for (const x of finite) {
    let b = span === 0 ? 0 : Math.min(bins.length - 1, Math.floor(((x - min) / span) * bins.length));
    bins[b]++;
  }
  return {
    count: finite.length,
    min,
    max,
    mean: sum / finite.length,
    rms: Math.sqrt(squares / finite.length),
    p50: percentile(sorted, 0.5),
    p95: percentile(sorted, 0.95),
    histogram: bins,
    edges: Array.from({ length: bins.length + 1 }, (_, i) => min + (span * i) / bins.length)
  };
}
