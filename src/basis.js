export function basisValues(knots, degree, t, side = 1) {
  const L = knots.length;
  const memo = new Map();

  function N(i, k) {
    if (i < 0 || i > L - k - 2) return 0;
    const key = i * 64 + k;
    if (memo.has(key)) return memo.get(key);

    let value = 0;
    if (k === 0) {
      if (knots[i] !== knots[i + 1]) {
        if (side < 0) {
          value = knots[i] < t && t <= knots[i + 1] ? 1 : 0;
        } else {
          value = knots[i] <= t && t < knots[i + 1] ? 1 : 0;
        }
      }
    } else {
      const d1 = knots[i + k] - knots[i];
      if (d1 !== 0) value += ((t - knots[i]) / d1) * N(i, k - 1);

      const d2 = knots[i + k + 1] - knots[i + 1];
      if (d2 !== 0) value += ((knots[i + k + 1] - t) / d2) * N(i + 1, k - 1);
    }

    memo.set(key, value);
    return value;
  }

  const count = L - degree - 1;
  return Array.from({ length: count }, (_, i) => N(i, degree));
}

export function basisDerivatives(knots, degree, t, side = 1) {
  const count = knots.length - degree - 1;
  if (degree === 0) return new Array(count).fill(0);

  const lower = basisValues(knots, degree - 1, t, side);
  const result = new Array(count).fill(0);

  for (let i = 0; i < count; i++) {
    const d1 = knots[i + degree] - knots[i];
    if (d1 !== 0) result[i] += (degree * lower[i]) / d1;

    const d2 = knots[i + degree + 1] - knots[i + 1];
    if (d2 !== 0) result[i] -= (degree * lower[i + 1]) / d2;
  }
  return result;
}

export function knotAxisInfo(knots, degree, count, t, requestedSide = 'auto') {
  const lo = knots[degree];
  const hi = knots[count];
  if (!(t >= lo && t <= hi)) {
    return { inDomain: false, lo, hi, t };
  }

  const indices = [];
  knots.forEach((k, i) => {
    if (k === t) indices.push(i);
  });
  const atKnot = indices.length > 0;
  const multiplicity = indices.length;
  const continuity = atKnot ? degree - multiplicity : Infinity;

  let side;
  let endpoint = null;
  const trace = [];

  if (t === lo) {
    side = 1;
    endpoint = 'low';
    trace.push({ code: 'ENDPOINT_LOW_INWARD_RIGHT', message: `参数域下端点 K[${degree}]=${lo}：强制使用右向内单侧规则。` });
  } else if (t === hi) {
    side = -1;
    endpoint = 'high';
    trace.push({ code: 'ENDPOINT_HIGH_INWARD_LEFT', message: `参数域上端点 K[${count}]=${hi}：强制使用左向内单侧规则。` });
  } else if (atKnot) {
    if (requestedSide === -1 || requestedSide === 'left') side = -1;
    else if (requestedSide === 1 || requestedSide === 'right') side = 1;
    else side = 1;

    if (multiplicity > degree) {
      trace.push({
        code: 'KNOT_DISCONTINUOUS',
        message: `内部节点 ${t} 重复度 ${multiplicity}=p+1；位置无唯一双侧值，当前选择${side < 0 ? '左' : '右'}侧分支。`
      });
    } else if (continuity === 0) {
      trace.push({
        code: 'KNOT_DERIVATIVE_ONESIDED',
        message: `内部节点 ${t} 为 C0；一阶导数不唯一，返回明确的${side < 0 ? '左' : '右'}单侧导数。`
      });
    } else {
      trace.push({
        code: 'KNOT_CONTINUOUS',
        message: `节点 ${t} 重复度 ${multiplicity}，C${continuity}；一阶导数双侧连续，仍按${side < 0 ? '左' : '右'}侧规则求值。`
      });
    }
  } else {
    side = 1;
  }

  return {
    inDomain: true,
    lo,
    hi,
    t,
    side,
    endpoint,
    atKnot,
    multiplicity,
    firstIndex: indices[0] ?? null,
    lastIndex: indices[indices.length - 1] ?? null,
    continuity,
    positionUnique: !atKnot || multiplicity <= degree || endpoint !== null,
    derivativeUnique: !atKnot || continuity >= 1,
    trace
  };
}

export function findDomainSpan(knots, degree, count, t, side) {
  const from = degree;
  const to = count - 1;
  if (side < 0) {
    for (let k = from; k <= to; k++) {
      if (knots[k] < knots[k + 1] && knots[k] < t && t <= knots[k + 1]) return k;
    }
  } else {
    for (let k = from; k <= to; k++) {
      if (knots[k] < knots[k + 1] && knots[k] <= t && t < knots[k + 1]) return k;
    }
  }
  if (t === knots[count]) return count - 1;
  return -1;
}

export function partitionUnityError(knots, degree, count, t, side) {
  const b = basisValues(knots, degree, t, side);
  let sum = 0;
  for (let i = 0; i < count; i++) sum += b[i] ?? 0;
  return Math.abs(sum - 1);
}
