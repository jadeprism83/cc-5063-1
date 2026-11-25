import { cross3, sub3, norm3, scale3 } from './math-utils.js';

export function exportOBJ(mesh, name = 'rational-spline-surface') {
  const used = [...new Set(mesh.indices)].sort((a, b) => a - b);
  const remap = new Map(used.map((oldIndex, newIndex) => [oldIndex, newIndex + 1]));
  const lines = [`o ${name.replace(/\s+/g, '_')}`, '# Rational B-spline tessellation; control net is not triangulated.'];

  for (const oldIndex of used) {
    const k = oldIndex * 3;
    lines.push(`v ${mesh.positions[k].toFixed(9)} ${mesh.positions[k + 1].toFixed(9)} ${mesh.positions[k + 2].toFixed(9)}`);
  }
  for (const oldIndex of used) {
    const k = oldIndex * 2;
    lines.push(`vt ${mesh.uvs[k].toFixed(6)} ${mesh.uvs[k + 1].toFixed(6)}`);
  }
  for (const oldIndex of used) {
    const k = oldIndex * 3;
    lines.push(`vn ${mesh.normals[k].toFixed(9)} ${mesh.normals[k + 1].toFixed(9)} ${mesh.normals[k + 2].toFixed(9)}`);
  }

  for (let i = 0; i < mesh.indices.length; i += 3) {
    const a = remap.get(mesh.indices[i]);
    const b = remap.get(mesh.indices[i + 1]);
    const c = remap.get(mesh.indices[i + 2]);
    lines.push(`f ${a}/${a}/${a} ${b}/${b}/${b} ${c}/${c}/${c}`);
  }
  return lines.join('\n') + '\n';
}

function triangleNormal(a, b, c) {
  const n = cross3(sub3(b, a), sub3(c, a));
  const len = norm3(n);
  return len > 0 ? scale3(n, 1 / len) : [0, 0, 0];
}

export function exportSTL(mesh) {
  const triangleCount = mesh.indices.length / 3;
  const buffer = new ArrayBuffer(84 + triangleCount * 50);
  const view = new DataView(buffer);
  const bytes = new Uint8Array(buffer);
  const encoder = new TextEncoder();
  const header = encoder.encode('RationalSplineStudio NURBS mesh');
  bytes.set(header.subarray(0, Math.min(80, header.length)), 0);
  view.setUint32(80, triangleCount, true);

  let offset = 84;
  const readPoint = index => {
    const k = index * 3;
    return [mesh.positions[k], mesh.positions[k + 1], mesh.positions[k + 2]];
  };
  for (let t = 0; t < triangleCount; t++) {
    const ia = mesh.indices[t * 3];
    const ib = mesh.indices[t * 3 + 1];
    const ic = mesh.indices[t * 3 + 2];
    const a = readPoint(ia);
    const b = readPoint(ib);
    const c = readPoint(ic);
    const n = triangleNormal(a, b, c);
    for (const value of [...n, ...a, ...b, ...c]) {
      view.setFloat32(offset, value, true);
      offset += 4;
    }
    view.setUint16(offset, 0, true);
    offset += 2;
  }
  return bytes;
}

export function downloadBlob(filename, data, type) {
  const blob = data instanceof Blob ? data : new Blob([data], { type });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
