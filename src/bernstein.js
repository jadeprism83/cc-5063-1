export function bernstein(degree, t) {
  let values = [1];
  for (let k = 1; k <= degree; k++) {
    const next = new Array(k + 1).fill(0);
    for (let i = 0; i <= k; i++) {
      if (i < k) next[i] += (1 - t) * values[i];
      if (i > 0) next[i] += t * values[i - 1];
    }
    values = next;
  }
  return values;
}

export function bernsteinDerivative(degree, t) {
  if (degree === 0) return [0];
  const lower = bernstein(degree - 1, t);
  return Array.from({ length: degree + 1 }, (_, i) => {
    const previous = i > 0 ? lower[i - 1] : 0;
    const current = i < lower.length ? lower[i] : 0;
    return degree * (previous - current);
  });
}

function deCasteljauLine(points, t) {
  let current = points.map(p => p.slice());
  while (current.length > 1) {
    const next = [];
    for (let i = 0; i < current.length - 1; i++) {
      next.push([
        (1 - t) * current[i][0] + t * current[i + 1][0],
        (1 - t) * current[i][1] + t * current[i + 1][1],
        (1 - t) * current[i][2] + t * current[i + 1][2],
        (1 - t) * current[i][3] + t * current[i + 1][3]
      ]);
    }
    current = next;
  }
  return current[0];
}

export function deCasteljau4(patch, s, t) {
  const qPlusOne = patch[0].length;
  const reducedU = [];
  for (let j = 0; j < qPlusOne; j++) {
    const line = patch.map(row => row[j]);
    reducedU.push(deCasteljauLine(line, s));
  }
  return deCasteljauLine(reducedU, t);
}
