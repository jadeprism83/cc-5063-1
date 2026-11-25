import { createExample } from './src/examples.js';
import { compileSurface } from './src/surface.js';
import {
  createTrimmedArchive,
  migrateArchive,
  checkRoundTrip,
  compareMigrationStrategies,
  DECAL_BAKED_SEAM,
  DECAL_ARC_LENGTH,
  MIGRATION_TARGET_DOMAIN
} from './src/kernel/trimming.js';
import { maxComponentError } from './src/math-utils.js';

function relDerivative(a, b) {
  return Math.hypot(...a.map((x, i) => x - b[i])) / Math.max(1, Math.hypot(...b));
}

console.log('=== 新内核 vs 旧内核（位置/导数误差执行证据）===');
for (const id of ['plane', 'cylinder', 'nonuniform', 'repeated-heavy']) {
  const s = compileSurface(createExample(id));
  let worstPos = 0;
  let worstDeriv = 0;
  for (const f of [0, 0.25, 0.5, 0.75, 1]) {
    for (const g of [0.1, 0.9]) {
      const u = s.domainU[0] + f * (s.domainU[1] - s.domainU[0]);
      const v = s.domainV[0] + g * (s.domainV[1] - s.domainV[0]);
      const neu = s.evaluate(u, v);
      const legacy = s.evaluateLegacy(u, v);
      if (!neu.ok || !legacy.ok) continue;
      worstPos = Math.max(worstPos, maxComponentError(neu.position, legacy.position));
      worstDeriv = Math.max(
        worstDeriv,
        relDerivative(neu.derivatives.u, legacy.derivatives.u),
        relDerivative(neu.derivatives.v, legacy.derivatives.v)
      );
    }
  }
  console.log(`${id.padEnd(16)} engine=${s.engine}  maxPosErr=${worstPos.toExponential(3)}  maxRelDerivErr=${worstDeriv.toExponential(3)}`);
}

console.log('\n=== 迁移可逆性（真实贴花交点 + 裁剪网格）===');
for (const mode of [DECAL_BAKED_SEAM, DECAL_ARC_LENGTH]) {
  const archive = createTrimmedArchive(mode);
  const { archive: migrated } = migrateArchive(archive);
  const rt = checkRoundTrip(archive, migrated);
  const decal = Math.max(...rt.decalErrors.map(e => e.positionError));
  const inv = Math.max(...rt.inverseLoopErrors.map(e => e.error));
  console.log(`${mode.padEnd(12)} domain=${migrated.kernel.parameterDomain.join('->')}  decalErr=${decal.toExponential(3)}  gridErr=${rt.gridMaxError.toExponential(3)}  inverseErr=${inv.toExponential(3)}  reversible=${rt.reversible}`);
}

console.log('\n=== 策略 A/B 对照（同形几何，不同节点域）===');
const strategies = compareMigrationStrategies(createTrimmedArchive(DECAL_ARC_LENGTH, { exampleId: 'nonuniform' }));
console.log(`strategyA domainU=${strategies.strategyA.domainU}  strategyB domainU=${strategies.strategyB.domainU}  target=${MIGRATION_TARGET_DOMAIN.join(',')}`);
console.log(`domainsDiffer=${strategies.domainsDiffer}  geometryIdentical=${strategies.geometryIdenticalAcrossDomains}  selected=${strategies.selected}`);

console.log('\n=== 退化输入显式状态（不伪造点）===');
const deg = compileSurface(createExample('zero-weight-endpoint'));
const r = deg.evaluate(deg.domainU[0], deg.domainV[0]);
console.log(`ok=${r.ok}  code=${r.code}  engine=${r.engine}  hasPosition=${r.position !== undefined}  limitCandidates=${Object.keys(r.limitCandidates).join('|')}`);
