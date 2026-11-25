import { distribution, norm3, sub3, dot3 } from './math-utils.js';
import { basisValues } from './basis.js';

function domainLength(surface, axis) {
  const d = axis === 'u' ? surface.domainU : surface.domainV;
  return d[1] - d[0];
}

function finiteDifferenceAxis(surface, axis, u, v, located) {
  const isU = axis === 'u';
  const info = isU ? located.infoU : located.infoV;
  const spanIndex = isU ? located.ku : located.kv;
  const knots = isU ? surface.U : surface.V;
  const domain = domainLength(surface, axis);
  const atBoundaryRule = info.atKnot || info.endpoint;

  let maxH = 1e-3 * domain;
  if (spanIndex >= 0) {
    const a = knots[spanIndex];
    const b = knots[spanIndex + 1];
    if (atBoundaryRule) {
      maxH = Math.min(maxH, 0.1 * (b - a));
    } else {
      const value = isU ? u : v;
      maxH = Math.min(maxH, 0.2 * Math.max(1e-14, Math.min(value - a, b - value)));
    }
  }

  function evaluateAt(delta) {
    const args = isU
      ? [u + delta, v, { sideU: delta === 0 ? info.side : undefined, sideV: located.infoV.side }]
      : [u, v + delta, { sideU: located.infoU.side, sideV: delta === 0 ? info.side : undefined }];
    return surface.evaluate(...args);
  }

  function derivativeAt(h) {
    const p0 = evaluateAt(0);
    if (!p0.ok) return null;
    const key = isU ? 'u' : 'v';

    if (atBoundaryRule && info.side < 0) {
      const p1 = evaluateAt(-h);
      const p2 = evaluateAt(-2 * h);
      if (!p1.ok || !p2.ok) return null;
      return [0, 1, 2].map(k => (3 * p0.position[k] - 4 * p1.position[k] + p2.position[k]) / (2 * h));
    }
    if (atBoundaryRule && info.side > 0) {
      const p1 = evaluateAt(h);
      const p2 = evaluateAt(2 * h);
      if (!p1.ok || !p2.ok) return null;
      return [0, 1, 2].map(k => (-3 * p0.position[k] + 4 * p1.position[k] - p2.position[k]) / (2 * h));
    }

    const pm = evaluateAt(-h);
    const pp = evaluateAt(h);
    if (!pm.ok || !pp.ok) return null;
    return [0, 1, 2].map(k => (pp.position[k] - pm.position[k]) / (2 * h));
  }

  const estimates = [];
  let h = maxH;
  for (let n = 0; n < 9; n++) {
    const coarse = derivativeAt(h);
    const half = derivativeAt(h / 2);
    if (coarse && half && [...coarse, ...half].every(Number.isFinite)) {
      const richardson = coarse.map((x, k) => (4 * half[k] - x) / 3);
      const score = norm3(sub3(richardson, half));
      estimates.push({ h, vector: richardson, raw: half, score });
    }
    h /= Math.sqrt(10);
  }
  if (!estimates.length) return { ok: false, code: `FD_SAMPLE_FAILED_${axis.toUpperCase()}` };

  estimates.sort((a, b) => a.score - b.score);
  const best = estimates[0];
  return {
    ok: true,
    vector: best.vector,
    h: best.h,
    score: best.score,
    side: info.side,
    candidates: estimates.length
  };
}

export function finiteDifferenceDerivatives(surface, u, v, options = {}) {
  const located = surface.spanAt(u, v, options);
  if (!located) return { ok: false, code: 'OUT_OF_DOMAIN' };

  const du = finiteDifferenceAxis(surface, 'u', u, v, located);
  if (!du.ok) return du;
  const dv = finiteDifferenceAxis(surface, 'v', u, v, located);
  if (!dv.ok) return dv;

  return {
    ok: true,
    u: du.vector,
    v: dv.vector,
    hu: du.h,
    hv: dv.h,
    scoreU: du.score,
    scoreV: dv.score,
    sideU: located.infoU.side,
    sideV: located.infoV.side
  };
}

function relativeVectorError(a, b, scale) {
  const abs = norm3(sub3(a, b));
  return { absolute: abs, relative: abs / Math.max(norm3(b), 1e-12, scale * 1e-12) };
}

export function compareAtPoint(surface, u, v, reference = null, options = {}) {
  const actual = surface.evaluate(u, v, options);
  const subdivided = surface.evaluateSubdivided(u, v, options);
  const fd = finiteDifferenceDerivatives(surface, u, v, options);
  const located = surface.spanAt(u, v, options);

  const result = {
    u, v,
    sideU: actual.sideU ?? located?.infoU?.side,
    sideV: actual.sideV ?? located?.infoV?.side,
    actual,
    subdivided,
    finiteDifference: fd,
    ok: actual.ok,
    errors: {}
  };

  if (located) {
    const bu = basisValues(surface.U, surface.p, u, located.infoU.side);
    const bv = basisValues(surface.V, surface.q, v, located.infoV.side);
    result.basisPartitionErrorU = Math.abs(bu.reduce((a, b) => a + b, 0) - 1);
    result.basisPartitionErrorV = Math.abs(bv.reduce((a, b) => a + b, 0) - 1);
  }

  if (actual.ok) {
    if (subdivided.ok) {
      result.errors.subdivisionPosition = relativeVectorError(subdivided.position, actual.position, surface.scale);
      result.errors.subdivisionDerivativeU = relativeVectorError(subdivided.derivatives.u, actual.derivatives.u, surface.scale);
      result.errors.subdivisionDerivativeV = relativeVectorError(subdivided.derivatives.v, actual.derivatives.v, surface.scale);
    }
    if (fd.ok) {
      result.errors.finiteDifferenceU = relativeVectorError(actual.derivatives.u, fd.u, surface.scale);
      result.errors.finiteDifferenceV = relativeVectorError(actual.derivatives.v, fd.v, surface.scale);
    }
    if (reference?.position) {
      const expectedPosition = reference.position(u, v);
      result.errors.analyticPosition = relativeVectorError(actual.position, expectedPosition, surface.scale);
    }
    if (reference?.derivatives) {
      const expectedDerivatives = reference.derivatives(u, v);
      result.errors.analyticDerivativeU = relativeVectorError(actual.derivatives.u, expectedDerivatives.u, surface.scale);
      result.errors.analyticDerivativeV = relativeVectorError(actual.derivatives.v, expectedDerivatives.v, surface.scale);
    }
    if (reference?.normal && actual.normal) {
      const n = reference.normal(u, v);
      result.normalDot = dot3(actual.normal, n);
      result.errors.analyticNormal = { absolute: Math.abs(1 - Math.abs(result.normalDot)), relative: Math.abs(1 - Math.abs(result.normalDot)) };
    }
  }

  return result;
}

export function runErrorDistribution(surface, reference = null, count = 25) {
  const [u0, u1] = surface.domainU;
  const [v0, v1] = surface.domainV;
  const samples = [];
  const metrics = {
    analyticPosition: [],
    analyticDerivative: [],
    finiteDifferenceDerivative: [],
    subdivisionPosition: [],
    subdivisionDerivative: [],
    basisPartition: []
  };
  const heatmap = [];
  let failures = 0;

  for (let j = 0; j < count; j++) {
    const v = v0 + ((j + 0.5) / count) * (v1 - v0);
    const row = [];
    for (let i = 0; i < count; i++) {
      const u = u0 + ((i + 0.5) / count) * (u1 - u0);
      const r = compareAtPoint(surface, u, v, reference);
      samples.push(r);
      if (!r.ok) failures++;
      const derivativeError = Math.max(
        r.errors.finiteDifferenceU?.relative ?? NaN,
        r.errors.finiteDifferenceV?.relative ?? NaN
      );
      row.push(Number.isFinite(derivativeError) ? derivativeError : null);

      if (r.errors.analyticPosition) metrics.analyticPosition.push(r.errors.analyticPosition.absolute);
      if (r.errors.analyticDerivativeU) {
        metrics.analyticDerivative.push(Math.max(r.errors.analyticDerivativeU.absolute, r.errors.analyticDerivativeV.absolute));
      }
      if (Number.isFinite(derivativeError)) metrics.finiteDifferenceDerivative.push(derivativeError);
      if (r.errors.subdivisionPosition) metrics.subdivisionPosition.push(r.errors.subdivisionPosition.absolute);
      if (r.errors.subdivisionDerivativeU) {
        metrics.subdivisionDerivative.push(Math.max(
          r.errors.subdivisionDerivativeU.relative,
          r.errors.subdivisionDerivativeV.relative
        ));
      }
      if (Number.isFinite(r.basisPartitionErrorU)) metrics.basisPartition.push(Math.max(r.basisPartitionErrorU, r.basisPartitionErrorV));
    }
    heatmap.push(row);
  }

  const stats = Object.fromEntries(Object.entries(metrics).map(([key, values]) => [key, distribution(values)]));
  return { count, samples, heatmap, stats, failures };
}

function distinctInternalKnots(knots, p, count) {
  const lo = knots[p];
  const hi = knots[count];
  return [...new Set(knots.filter(k => k > lo && k < hi))];
}

export function buildKeyLocations(surface) {
  const [u0, u1] = surface.domainU;
  const [v0, v1] = surface.domainV;
  const us = [u0, u0 + 0.25 * (u1 - u0), (u0 + u1) / 2, u0 + 0.75 * (u1 - u0), u1];
  const vs = [v0, (v0 + v1) / 2, v1];
  const locations = [];

  for (const u of us) {
    for (const v of vs) locations.push({ u, v, label: `常规关键点 u=${u.toFixed(4)}, v=${v.toFixed(4)}` });
  }

  for (const knot of distinctInternalKnots(surface.U, surface.p, surface.nu)) {
    locations.push({ u: knot, v: (v0 + v1) / 2, sideU: 'left', label: `U 重复节点 ${knot} 左侧` });
    locations.push({ u: knot, v: (v0 + v1) / 2, sideU: 'right', label: `U 重复节点 ${knot} 右侧` });
    const locatedR = surface.spanAt(knot, (v0 + v1) / 2, { sideU: 'right' });
    const locatedL = surface.spanAt(knot, (v0 + v1) / 2, { sideU: 'left' });
    if (locatedR?.ku >= 0) {
      const k = locatedR.ku;
      const delta = 0.01 * (surface.U[k + 1] - surface.U[k]);
      locations.push({ u: knot + delta, v: (v0 + v1) / 2, label: `U 节点右邻 ${knot}+δ` });
    }
    if (locatedL?.ku >= 0) {
      const k = locatedL.ku;
      const delta = 0.01 * (surface.V ? surface.U[k + 1] - surface.U[k] : 1);
      locations.push({ u: knot - delta, v: (v0 + v1) / 2, label: `U 节点左邻 ${knot}-δ` });
    }
  }

  for (const knot of distinctInternalKnots(surface.V, surface.q, surface.nv)) {
    locations.push({ u: (u0 + u1) / 2, v: knot, sideV: 'left', label: `V 重复节点 ${knot} 左侧` });
    locations.push({ u: (u0 + u1) / 2, v: knot, sideV: 'right', label: `V 重复节点 ${knot} 右侧` });
  }

  return locations;
}

export function runKeyChecks(surface, reference = null) {
  const records = [];
  for (const loc of buildKeyLocations(surface)) {
    const comparison = compareAtPoint(surface, loc.u, loc.v, reference, loc);
    const checks = [];
    const add = (name, pass, detail = '') => checks.push({ name, pass, detail });

    add('参数域与分支可追踪', Boolean(comparison.actual?.branch || comparison.actual?.infoU), comparison.actual?.trace?.map(t => t.message).join(' / ') ?? '');
    add('基函数单位分解', Math.max(comparison.basisPartitionErrorU ?? 1, comparison.basisPartitionErrorV ?? 1) < 1e-10,
      `U=${comparison.basisPartitionErrorU?.toExponential?.(2)}, V=${comparison.basisPartitionErrorV?.toExponential?.(2)}`);

    if (!comparison.ok) {
      add('退化分母未被伪装成曲面点', comparison.actual?.code === 'DEGENERATE_DENOMINATOR', comparison.actual?.message ?? '');
      records.push({ ...loc, comparison, checks, status: checks.every(c => c.pass) ? 'branch' : 'fail' });
      continue;
    }

    const fdU = comparison.errors.finiteDifferenceU;
    const fdV = comparison.errors.finiteDifferenceV;
    add('一阶导数 vs 单侧数值差分', Boolean(fdU && fdV && fdU.relative < 2e-5 && fdV.relative < 2e-5),
      `relU=${fdU?.relative?.toExponential?.(2)}, relV=${fdV?.relative?.toExponential?.(2)}`);

    const sp = comparison.errors.subdivisionPosition;
    add('齐次求值 vs 局部细分位置', Boolean(sp && sp.absolute < 1e-7 * Math.max(1, surface.scale)),
      `abs=${sp?.absolute?.toExponential?.(3)}`);

    if (comparison.errors.subdivisionDerivativeU) {
      const su = comparison.errors.subdivisionDerivativeU.relative;
      const sv = comparison.errors.subdivisionDerivativeV.relative;
      add('齐次/细分一阶导数一致', su < 2e-5 && sv < 2e-5, `relU=${su?.toExponential?.(2)}, relV=${sv?.toExponential?.(2)}`);
    }

    if (comparison.errors.analyticPosition) {
      add('解析位置误差', comparison.errors.analyticPosition.absolute < 1e-9 * Math.max(1, surface.scale),
        comparison.errors.analyticPosition.absolute.toExponential(3));
    }
    if (comparison.errors.analyticDerivativeU) {
      const u = comparison.errors.analyticDerivativeU.absolute;
      const v = comparison.errors.analyticDerivativeV.absolute;
      add('解析一阶导数误差', u < 1e-8 * Math.max(1, surface.scale) && v < 1e-8 * Math.max(1, surface.scale),
        `U=${u.toExponential(3)}, V=${v.toExponential(3)}`);
    }
    if (comparison.normalDot !== undefined) {
      add('解析法线误差', Math.abs(1 - Math.abs(comparison.normalDot)) < 1e-7,
        `dot=${comparison.normalDot.toFixed(9)}`);
    }
    if (comparison.actual.normalStatus !== 'OK') {
      add('退化法线返回空分支', comparison.actual.normal === null, comparison.actual.normalStatus);
    }

    const failed = checks.some(c => !c.pass);
    records.push({ ...loc, comparison, checks, status: failed ? 'fail' : 'pass' });
  }
  return records;
}
