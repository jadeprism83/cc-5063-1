import { validateSpec, normalizeSpec } from './validation.js';
import { knotAxisInfo, findDomainSpan } from './basis.js';
import { evaluateSubdivided } from './subdivision.js';
import { legacyEvaluate } from './legacy-engine.js';
import { buildKernelSurface, adaptEvaluation } from './kernel/adapter.js';
import { KERNEL_NAME, KERNEL_VERSION } from './kernel/nurbs-kernel.js';

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
    this.engine = KERNEL_NAME;
    this.kernelVersion = KERNEL_VERSION;
    this.kernelBundle = buildKernelSurface(spec, { reparametrize: false });
  }

  evaluate(u, v, options = {}) {
    return adaptEvaluation(this, u, v, options);
  }

  evaluateLegacy(u, v, options = {}) {
    return legacyEvaluate(this, u, v, options);
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

  isoparametricCurve(fixed, value, segments = 48) {
    const points = [];
    const failures = [];
    if (fixed === 'v') {
      const [u0, u1] = this.domainU;
      for (let i = 0; i <= segments; i++) {
        const u = u0 + ((u1 - u0) * i) / segments;
        const r = this.evaluate(u, value);
        if (r.ok) points.push(r.position);
        else failures.push({ u, v: value, code: r.code });
      }
    } else {
      const [v0, v1] = this.domainV;
      for (let i = 0; i <= segments; i++) {
        const v = v0 + ((v1 - v0) * i) / segments;
        const r = this.evaluate(value, v);
        if (r.ok) points.push(r.position);
        else failures.push({ u: value, v, code: r.code });
      }
    }
    return { fixed, value, points, failures, engine: this.engine };
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
          failures.push({ i, j, u, v, code: r.code ?? r.normalStatus, message: r.message, engine: r.engine });
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

    return { positions, normals, uvs, indices, valid, evaluations, failures, cols, rows, engine: this.engine, kernelVersion: this.kernelVersion };
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
