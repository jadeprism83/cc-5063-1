import { mkdir, writeFile } from 'node:fs/promises';
import { join, resolve } from 'node:path';
import { createExample, EXAMPLE_IDS, REFERENCES } from './src/examples.js';
import { compileSurface } from './src/surface.js';
import { exportOBJ, exportSTL } from './src/exporters.js';
import { runKeyChecks, runErrorDistribution } from './src/verification.js';

const requested = process.argv[2] || 'all';
const outDir = resolve(process.argv[3] || './exports');
const ids = requested === 'all' ? EXAMPLE_IDS : [requested];

await mkdir(outDir, { recursive: true });

for (const id of ids) {
  const spec = createExample(id);
  const surface = compileSurface(spec);
  const mesh = surface.tessellate(spec.tessellation?.uSegments ?? 80, spec.tessellation?.vSegments ?? 40);
  const reference = REFERENCES[id] || null;
  const keyChecks = runKeyChecks(surface, reference);
  const distribution = runErrorDistribution(surface, reference, 21);
  const summary = {
    id,
    name: spec.name,
    domain: { u: surface.domainU, v: surface.domainV },
    controlPoints: { u: surface.nu, v: surface.nv },
    degrees: { u: surface.p, v: surface.q },
    mesh: { uSegments: mesh.cols - 1, vSegments: mesh.rows - 1, triangles: mesh.indices.length / 3, skippedSamples: mesh.failures.length },
    keyChecks: {
      pass: keyChecks.filter(r => r.status === 'pass').length,
      branch: keyChecks.filter(r => r.status === 'branch').length,
      fail: keyChecks.filter(r => r.status === 'fail').length,
      failures: keyChecks.filter(r => r.status === 'fail').map(r => ({ label: r.label, checks: r.checks.filter(c => !c.pass) }))
    },
    stats: Object.fromEntries(Object.entries(distribution.stats).map(([k, v]) => [k, {
      count: v.count, min: v.min, max: v.max, mean: v.mean, rms: v.rms, p95: v.p95
    }]))
  };

  await writeFile(join(outDir, `${id}.json`), JSON.stringify({ format: 'rational-spline-studio-project', version: 1, spec }, null, 2));
  await writeFile(join(outDir, `${id}.obj`), exportOBJ(mesh, id));
  await writeFile(join(outDir, `${id}.stl`), Buffer.from(exportSTL(mesh)));
  await writeFile(join(outDir, `${id}.report.json`), JSON.stringify(summary, null, 2));
  console.log(id, summary.keyChecks, `triangles=${summary.mesh.triangles}`);
}

console.log(`Exported ${ids.length} model(s) to ${outDir}`);
