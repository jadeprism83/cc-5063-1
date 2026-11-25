const EPS = 1e-12;
const FD_REL = 1e-6;

export function isFiniteNumber(x) {
  return typeof x === 'number' && Number.isFinite(x);
}

export function cloneSpec(spec) {
  return JSON.parse(JSON.stringify(spec));
}

function pointAt(p, i, j) {
  const raw = p.points[i][j];
  if (Array.isArray(raw)) return { x: raw[0], y: raw[1], z: raw[2] };
  return { x: raw.x, y: raw.y, z: raw.z };
}

export function validateSpec(input) {
  const errors = [];
  const warnings = [];
  const pushError = (code, message, path) => errors.push({ code, message, path });
  const pushWarning = (code, message, path) => warnings.push({ code, message, path });

  if (!input || typeof input !== 'object') {
    pushError('SPEC_NOT_OBJECT', '工程必须是 JSON 对象。');
    return { valid: false, errors, warnings, info: null };
  }

  const p = input.degreeU;
  const q = input.degreeV;
  if (!Number.isInteger(p) || p < 0 || p > 31) {
    pushError('INVALID_DEGREE_U', 'U 次数必须是 0..31 的整数。', 'degreeU');
  }
  if (!Number.isInteger(q) || q < 0 || q > 31) {
    pushError('INVALID_DEGREE_V', 'V 次数必须是 0..31 的整数。', 'degreeV');
  }

  if (!Array.isArray(input.points) || input.points.length === 0) {
    pushError('INVALID_POINTS_U', '控制点必须是非空二维数组。', 'points');
    return { valid: false, errors, warnings, info: null };
  }
  const nu = input.points.length;
  let nv = null;
  for (let i = 0; i < nu; i++) {
    const row = input.points[i];
    if (!Array.isArray(row) || row.length === 0) {
      pushError('INVALID_POINTS_V', `第 ${i} 行控制点必须是非空数组。`, `points[${i}]`);
      continue;
    }
    if (nv === null) nv = row.length;
    if (row.length !== nv) {
      pushError('IRREGULAR_GRID', `控制点网格不是矩形：第 0 行有 ${nv} 个，第 ${i} 行有 ${row.length} 个。`, `points[${i}]`);
    }
    for (let j = 0; j < row.length; j++) {
      const pt = row[j];
      const x = Array.isArray(pt) ? pt[0] : pt?.x;
      const y = Array.isArray(pt) ? pt[1] : pt?.y;
      const z = Array.isArray(pt) ? pt[2] : pt?.z;
      if (![x, y, z].every(isFiniteNumber)) {
        pushError('NON_FINITE_POINT', `控制点 [${i},${j}] 含非有限坐标。`, `points[${i}][${j}]`);
      }
    }
  }

  if (!Array.isArray(input.weights) || input.weights.length !== nu) {
    pushError('INVALID_WEIGHT_GRID', '权重必须是与控制点同维的二维数组。', 'weights');
  } else {
    if (nv !== null) {
      for (let i = 0; i < nu; i++) {
        if (!Array.isArray(input.weights[i]) || input.weights[i].length !== nv) {
          pushError('INVALID_WEIGHT_ROW', `权重行 ${i} 长度必须为 ${nv}。`, `weights[${i}]`);
        } else {
          for (let j = 0; j < nv; j++) {
            const w = input.weights[i][j];
            if (!isFiniteNumber(w)) {
              pushError('NON_FINITE_WEIGHT', `权重 [${i},${j}] 必须是有限数。`, `weights[${i}][${j}]`);
            }
          }
        }
      }
    }
  }

  for (const axis of ['U', 'V']) {
    const degree = axis === 'U' ? p : q;
    const count = axis === 'U' ? nu : nv;
    const knots = axis === 'U' ? input.knotsU : input.knotsV;
    const key = `knots${axis}`;
    if (!Number.isInteger(degree) || degree < 0 || count === null) continue;
    if (degree + 1 > count) {
      pushError('DEGREE_TOO_HIGH', `${axis} 向次数 ${degree} 需要至少 ${degree + 1} 个控制点，当前只有 ${count} 个。`, `degree${axis}`);
    }
    const expectedLength = count + degree + 1;
    if (!Array.isArray(knots)) {
      pushError('KNOTS_NOT_ARRAY', `${axis} 向节点向量必须是数组。`, key);
      continue;
    }
    if (knots.length !== expectedLength) {
      pushError('KNOT_LENGTH', `${axis} 向节点长度应为 ${count}+${degree}+1=${expectedLength}，实际为 ${knots.length}。`, key);
    }
    if (!knots.every(isFiniteNumber)) {
      pushError('KNOT_NON_FINITE', `${axis} 向节点必须全部是有限数。`, key);
      continue;
    }
    for (let i = 1; i < knots.length; i++) {
      if (knots[i] < knots[i - 1]) {
        pushError('KNOT_ORDER', `${axis} 向节点在索引 ${i - 1}..${i} 处递减：${knots[i - 1]} > ${knots[i]}。`, `${key}[${i}]`);
        break;
      }
    }
    if (knots.length === expectedLength) {
      const lo = knots[degree];
      const hi = knots[count];
      if (!(hi > lo)) {
        pushError('EMPTY_DOMAIN', `${axis} 向参数域 [K[${degree}], K[${count}]] 必须非空，实际为 [${lo}, ${hi}]。`, key);
      }
      const mult = knotMultiplicities(knots);
      for (const [value, r, first] of mult) {
        if (r > degree + 1) {
          pushError('KNOT_MULTIPLICITY', `${axis} 向节点 ${formatNumber(value)} 的重复度为 ${r}，超过次数+1=${degree + 1}；基函数/连续性无效。`, `${key}[${first}]`);
        } else if (r > degree && degree > 0) {
          pushWarning('DISCONTINUOUS_KNOT', `${axis} 向节点 ${formatNumber(value)} 重复度 ${r}=p+1；曲面位置在内部节点处可能不连续，求值必须指定单侧。`, key);
        } else if (r > 1) {
          pushWarning('REDUCED_CONTINUITY', `${axis} 向节点 ${formatNumber(value)} 重复度 ${r}，连续性为 C^{${degree - r}}；一阶导数需按单侧解释。`, key);
        }
      }
    }
  }

  let zeroWeights = [];
  let negativeWeights = [];
  let minPositive = Infinity;
  let maxAbs = 0;
  if (Array.isArray(input.weights) && nv !== null) {
    for (let i = 0; i < nu; i++) {
      if (!Array.isArray(input.weights[i])) continue;
      for (let j = 0; j < nv; j++) {
        const w = input.weights[i][j];
        if (!isFiniteNumber(w)) continue;
        maxAbs = Math.max(maxAbs, Math.abs(w));
        if (w === 0) zeroWeights.push([i, j]);
        if (w < 0) negativeWeights.push([i, j]);
        if (w > 0) minPositive = Math.min(minPositive, w);
      }
    }
  }
  if (zeroWeights.length) {
    pushWarning('ZERO_WEIGHT', `存在 ${zeroWeights.length} 个零权重控制点；若其齐次基函数贡献非零，分母分支必须可追踪。`, 'weights');
  }
  if (negativeWeights.length) {
    pushWarning('NEGATIVE_WEIGHT', `存在 ${negativeWeights.length} 个负权重；有理分母可能抵消并产生极点。`, 'weights');
  }
  const weightRatio = Number.isFinite(minPositive) && minPositive > 0 ? maxAbs / minPositive : Infinity;
  if (weightRatio > 1e6) {
    pushWarning('LARGE_WEIGHT_SPREAD', `权重最大绝对值/最小正权重=${weightRatio.toExponential(3)}；应同时检查齐次求值与局部细分。`, 'weights');
  }

  const valid = errors.length === 0;
  return {
    valid,
    errors,
    warnings,
    info: valid ? {
      nu,
      nv,
      degreeU: p,
      degreeV: q,
      domainU: [input.knotsU[p], input.knotsU[nu]],
      domainV: [input.knotsV[q], input.knotsV[nv]],
      zeroWeights,
      negativeWeights,
      weightRatio,
      largeWeightSpread: weightRatio > 1e6,
      scale: computeScale(input, nu, nv)
    } : null
  };
}

function computeScale(input, nu, nv) {
  let maxNorm = 1;
  for (let i = 0; i < nu; i++) {
    for (let j = 0; j < nv; j++) {
      const pt = pointAt(input, i, j);
      if ([pt.x, pt.y, pt.z].every(isFiniteNumber)) {
        maxNorm = Math.max(maxNorm, Math.hypot(pt.x, pt.y, pt.z));
      }
    }
  }
  return maxNorm;
}

export function knotMultiplicities(knots) {
  const out = [];
  for (let i = 0; i < knots.length;) {
    let j = i + 1;
    while (j < knots.length && knots[j] === knots[i]) j++;
    out.push([knots[i], j - i, i]);
    i = j;
  }
  return out;
}

export function formatNumber(x, digits = 6) {
  if (!Number.isFinite(x)) return String(x);
  if (Math.abs(x) < 1e-12 && x !== 0) return x.toExponential(2);
  const rounded = Number(x.toFixed(digits));
  return String(rounded);
}

export function normalizeSpec(input) {
  const v = validateSpec(input);
  if (!v.valid) {
    const err = new Error(v.errors.map(e => `${e.code}: ${e.message}`).join('\n'));
    err.validation = v;
    throw err;
  }
  const nu = input.points.length;
  const nv = input.points[0].length;
  const points = Array.from({ length: nu }, (_, i) =>
    Array.from({ length: nv }, (_, j) => {
      const pt = pointAt(input, i, j);
      return [pt.x, pt.y, pt.z];
    })
  );
  const weights = Array.from({ length: nu }, (_, i) =>
    Array.from({ length: nv }, (_, j) => input.weights[i][j])
  );
  return {
    id: input.id || 'surface',
    name: input.name || '未命名有理样条曲面',
    degreeU: input.degreeU,
    degreeV: input.degreeV,
    knotsU: input.knotsU.slice(),
    knotsV: input.knotsV.slice(),
    points,
    weights,
    info: v.info
  };
}
