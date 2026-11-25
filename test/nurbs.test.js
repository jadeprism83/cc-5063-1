import test from 'node:test';
import assert from 'node:assert/strict';
import { compileSurface } from '../src/surface.js';
import { validateSpec } from '../src/validation.js';
import { createExample, REFERENCES } from '../src/examples.js';
import { compareAtPoint, runKeyChecks, finiteDifferenceDerivatives } from '../src/verification.js';
import { exportOBJ, exportSTL } from '../src/exporters.js';
import { basisValues } from '../src/basis.js';

function surface(id) {
  return compileSurface(createExample(id));
}

test('validation rejects knot length, order, multiplicity and empty domain', () => {
  const base = createExample('plane');
  assert.equal(validateSpec(base).valid, true);

  const badLength = structuredClone(base);
  badLength.knotsU.push(2);
  assert.ok(validateSpec(badLength).errors.some(e => e.code === 'KNOT_LENGTH'));

  const badOrder = structuredClone(base);
  badOrder.knotsU = [0, 0, 1, 0.5, 1, 1];
  assert.ok(validateSpec(badOrder).errors.some(e => e.code === 'KNOT_ORDER'));

  const badMult = structuredClone(base);
  badMult.knotsU = [0, 0, 0, 0, 1, 1];
  assert.ok(validateSpec(badMult).errors.some(e => e.code === 'KNOT_MULTIPLICITY'));

  const emptyDomain = structuredClone(base);
  emptyDomain.knotsU = [1, 1, 1, 1, 2, 2];
  assert.ok(validateSpec(emptyDomain).errors.some(e => e.code === 'EMPTY_DOMAIN'));

  const badWeight = structuredClone(base);
  badWeight.weights[0][0] = Number.NaN;
  assert.ok(validateSpec(badWeight).errors.some(e => e.code === 'NON_FINITE_WEIGHT'));
});

test('basis forms partition of unity on resolved endpoint and one-sided knot branches', () => {
  const knots = [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1];
  for (const t of [0, 0.25, 0.5, 0.75, 1]) {
    const sides = t === 0 ? [1] : t === 1 ? [-1] : [-1, 1];
    for (const side of sides) {
      const b = basisValues(knots, 3, t, side);
      assert.ok(Math.abs(b.reduce((a, x) => a + x, 0) - 1) < 1e-12, `t=${t}, side=${side}`);
    }
  }
});

test('plane evaluates analytic position, derivatives and normal', () => {
  const s = surface('plane');
  const ref = REFERENCES.plane;
  for (const [u, v] of [[0, 0], [1, 1], [0.31, 0.72], [0, 0.5], [0.5, 1]]) {
    const r = s.evaluate(u, v);
    assert.equal(r.ok, true);
    const p = ref.position(u, v);
    assert.ok(Math.max(...r.position.map((x, i) => Math.abs(x - p[i]))) < 1e-12);
    assert.ok(Math.hypot(...r.derivatives.u.map((x, i) => x - ref.derivatives().u[i])) < 1e-12);
    assert.ok(Math.hypot(...r.derivatives.v.map((x, i) => x - ref.derivatives().v[i])) < 1e-12);
    assert.ok(Math.abs(r.normal[2] - 1) < 1e-12);
  }
});

test('rational cylinder is analytically on radius and matches angle parameterization', () => {
  const s = surface('cylinder');
  const ref = REFERENCES.cylinder;
  for (const [u, v] of [[0, 0], [0.17, 0.2], [0.5, 1], [0.83, 0.8], [1, 1]]) {
    const r = s.evaluate(u, v);
    assert.equal(r.ok, true);
    const expected = ref.position(u, v);
    assert.ok(Math.hypot(...r.position.map((x, i) => x - expected[i])) < 1e-10);
    const radius = Math.hypot(r.position[0], r.position[1]);
    assert.ok(Math.abs(radius - 1.2) < 1e-11);
    assert.ok(Math.abs(r.normal[0] ** 2 + r.normal[1] ** 2 - 1) < 1e-12);
    assert.ok(Math.abs(r.normal[2]) < 1e-12);
  }
});

test('analytic derivatives agree with one-sided finite differences at endpoints and knots', () => {
  for (const id of ['plane', 'cylinder', 'nonuniform']) {
    const s = surface(id);
    const [u0, u1] = s.domainU;
    const [v0, v1] = s.domainV;
    const points = [
      [u0, v0], [u1, v1], [(u0 + u1) / 2, (v0 + v1) / 2],
      ...[...new Set(s.U.filter(x => x > u0 && x < u1))].map(u => [u, (v0 + v1) / 2])
    ];
    for (const [u, v] of points) {
      const fd = finiteDifferenceDerivatives(s, u, v);
      assert.equal(fd.ok, true, `${id} fd ${u},${v}: ${fd.code}`);
      const r = s.evaluate(u, v);
      const eu = Math.hypot(...r.derivatives.u.map((x, i) => x - fd.u[i]));
      const ev = Math.hypot(...r.derivatives.v.map((x, i) => x - fd.v[i]));
      assert.ok(eu / Math.max(1, Math.hypot(...fd.u)) < 2e-5, `${id} U derivative ${eu}`);
      assert.ok(ev / Math.max(1, Math.hypot(...fd.v)) < 2e-5, `${id} V derivative ${ev}`);
    }
  }
});

test('homogeneous evaluation matches local knot-insertion subdivision', () => {
  for (const id of ['plane', 'cylinder', 'nonuniform', 'repeated-heavy']) {
    const s = surface(id);
    const [u0, u1] = s.domainU;
    const [v0, v1] = s.domainV;
    for (const f of [0, 0.13, 0.5, 0.5001, 0.87, 1]) {
      for (const g of [0, 0.37, 1]) {
        const u = u0 + f * (u1 - u0);
        const v = v0 + g * (v1 - v0);
        const c = compareAtPoint(s, u, v, REFERENCES[id]);
        assert.equal(c.actual.ok, true, `${id} actual ${f},${g}`);
        assert.equal(c.subdivided.ok, true, `${id} subdivided ${f},${g}: ${c.subdivided.code}`);
        assert.ok(c.errors.subdivisionPosition.absolute < 1e-7 * Math.max(1, s.scale),
          `${id} position mismatch ${c.errors.subdivisionPosition.absolute}`);
      }
    }
  }
});

test('key checks pass including repeated knot left and right branches', () => {
  const s = surface('repeated-heavy');
  const records = runKeyChecks(s);
  const failures = records.filter(r => r.status === 'fail');
  assert.equal(failures.length, 0, failures.map(f => `${f.label}: ${f.checks.filter(c => !c.pass).map(c => c.name + ' ' + c.detail).join('; ')}`).join('\n'));
  const knotRecords = records.filter(r => r.label.includes('重复节点 0.5'));
  assert.ok(knotRecords.length >= 2);
  for (const r of knotRecords) {
    assert.ok(r.comparison.actual.trace.some(t => t.code.startsWith('KNOT_')));
  }
});

test('zero endpoint weight has explicit traceable denominator branch and no normal', () => {
  const s = surface('zero-weight-endpoint');
  const r = s.evaluate(s.domainU[0], s.domainV[0]);
  assert.equal(r.ok, false);
  assert.equal(r.code, 'DEGENERATE_DENOMINATOR');
  assert.ok(r.trace.some(t => t.code === 'DEGENERATE_DENOMINATOR'));
  assert.ok(r.limitCandidates.alongU || r.limitCandidates.alongV);
  assert.equal(r.normal, undefined);

  const inside = s.evaluate(0.02, 0.02);
  assert.equal(inside.ok, true);
  assert.equal(inside.normalStatus, 'OK');
});

test('degenerate tangent returns null normal rather than arbitrary unit vector', () => {
  const spec = createExample('plane');
  spec.points = spec.points.map(row => row.map(p => ({ ...p, z: 0 })));
  for (let i = 1; i < spec.points.length; i++) {
    spec.points[i] = spec.points[0].map(p => ({ ...p }));
  }
  const s = compileSurface(spec);
  const r = s.evaluate(0.5, 0.5);
  assert.equal(r.ok, true);
  assert.equal(r.normal, null);
  assert.notEqual(r.normalStatus, 'OK');
});

test('mesh exporters contain actual triangles and valid binary STL size', () => {
  const s = surface('plane');
  const mesh = s.tessellate(8, 6);
  assert.ok(mesh.indices.length > 0);
  const obj = exportOBJ(mesh, 'plane');
  assert.match(obj, /^v /m);
  assert.match(obj, /^f /m);
  const stl = exportSTL(mesh);
  const triangles = mesh.indices.length / 3;
  assert.equal(stl.byteLength, 84 + triangles * 50);
});
