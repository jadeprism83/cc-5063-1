export const EXAMPLE_IDS = ['plane', 'cylinder', 'nonuniform', 'repeated-heavy', 'zero-weight-endpoint'];

export function createExample(id) {
  if (id === 'plane') return planeExample();
  if (id === 'cylinder') return cylinderExample();
  if (id === 'nonuniform') return nonuniformExample();
  if (id === 'repeated-heavy') return repeatedHeavyExample();
  if (id === 'zero-weight-endpoint') return zeroWeightEndpointExample();
  return planeExample();
}

function grid(nu, nv, fn) {
  return Array.from({ length: nu }, (_, i) =>
    Array.from({ length: nv }, (_, j) => fn(i, j))
  );
}

function constantWeights(nu, nv, value = 1) {
  return grid(nu, nv, () => value);
}

function planeExample() {
  const nu = 3;
  const nv = 3;
  return {
    id: 'plane',
    name: '平面：二次 Bézier 网格',
    degreeU: 2,
    degreeV: 2,
    knotsU: [0, 0, 0, 1, 1, 1],
    knotsV: [0, 0, 0, 1, 1, 1],
    points: grid(nu, nv, (i, j) => ({ x: -2 + 2 * i, y: -1 + j, z: 0 })),
    weights: constantWeights(nu, nv, 1),
    tessellation: { uSegments: 48, vSegments: 32 }
  };
}

function cylinderExample() {
  const nu = 3;
  const nv = 2;
  const R = 1.2;
  const H = 2;
  const xy = [
    [R, 0],
    [R, R],
    [0, R]
  ];
  const uWeights = [1, Math.SQRT1_2, 1];
  return {
    id: 'cylinder',
    name: '圆柱片：有理二次 90° 弧',
    degreeU: 2,
    degreeV: 1,
    knotsU: [0, 0, 0, 1, 1, 1],
    knotsV: [0, 0, 1, 1],
    points: grid(nu, nv, (i, j) => ({ x: xy[i][0], y: xy[i][1], z: -H / 2 + H * j })),
    weights: grid(nu, nv, (i) => uWeights[i]),
    tessellation: { uSegments: 72, vSegments: 16 }
  };
}

function nonuniformExample() {
  const nu = 7;
  const nv = 4;
  return {
    id: 'nonuniform',
    name: '非均匀节点：三次 U / 二次 V',
    degreeU: 3,
    degreeV: 2,
    knotsU: [0, 0, 0, 0, 0.18, 0.43, 0.72, 1, 1, 1, 1],
    knotsV: [0, 0, 0, 0.38, 1, 1, 1],
    points: grid(nu, nv, (i, j) => ({
      x: -3 + i,
      y: -1.5 + j,
      z: 0.42 * Math.sin(i * 1.13) * Math.cos(j * 0.91) + 0.12 * (i - 3) * (j - 1.5)
    })),
    weights: grid(nu, nv, (i, j) => 1 + 0.08 * Math.sin(i + j * 1.7)),
    tessellation: { uSegments: 96, vSegments: 40 }
  };
}

function repeatedHeavyExample() {
  const nu = 6;
  const nv = 2;
  const y = [0, 0.85, 1.15, 1.05, 0.75, 0];
  const w = [1, 1e-4, 1e3, 1e-3, 10, 1];
  return {
    id: 'repeated-heavy',
    name: '重复节点 × 极端权重（齐次/细分对照）',
    degreeU: 3,
    degreeV: 1,
    knotsU: [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1],
    knotsV: [0, 0, 1, 1],
    points: grid(nu, nv, (i, j) => ({
      x: -2 + 0.8 * i,
      y: y[i],
      z: j === 0 ? -0.55 : 0.55
    })),
    weights: grid(nu, nv, (i) => w[i]),
    tessellation: { uSegments: 120, vSegments: 12 }
  };
}

function zeroWeightEndpointExample() {
  const spec = planeExample();
  return {
    ...spec,
    id: 'zero-weight-endpoint',
    name: '零权重端点：可追踪退化分母分支',
    weights: spec.weights.map((row, i) => row.map((w, j) => (i === 0 && j === 0 ? 0 : w)))
  };
}

export const REFERENCES = {
  plane: {
    position: (u, v) => [-2 + 4 * u, -1 + 2 * v, 0],
    derivatives: () => ({ u: [4, 0, 0], v: [0, 2, 0] }),
    normal: () => [0, 0, 1]
  },
  cylinder: (() => {
    const R = 1.2;
    const H = 2;
    const h = Math.SQRT2;
    function theta(u) {
      const denominator = h + (1 - h) * u;
      return 2 * Math.atan(u / denominator);
    }
    return {
      position(u, v) {
        const t = theta(u);
        return [R * Math.cos(t), R * Math.sin(t), -H / 2 + H * v];
      },
      derivatives(u, v) {
        const denominator = h + (1 - h) * u;
        const s = u / denominator;
        const sprime = h / (denominator * denominator);
        const tp = (2 * sprime) / (1 + s * s);
        const t = theta(u);
        return {
          u: [-R * Math.sin(t) * tp, R * Math.cos(t) * tp, 0],
          v: [0, 0, H]
        };
      },
      normal(u) {
        const t = theta(u);
        return [Math.cos(t), Math.sin(t), 0];
      }
    };
  })()
};
