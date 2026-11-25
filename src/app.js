import { createExample, EXAMPLE_IDS, REFERENCES } from './examples.js';
import { compileSurface } from './surface.js';
import { validateSpec } from './validation.js';
import { CanvasViewer } from './viewer.js';
import { runErrorDistribution, runKeyChecks, compareAtPoint } from './verification.js';
import { exportOBJ, exportSTL, downloadBlob } from './exporters.js';
import { formatNumber } from './validation.js';
import { KERNEL_NAME, KERNEL_VERSION } from './kernel/nurbs-kernel.js';
import { maxComponentError } from './math-utils.js';

const $ = id => document.getElementById(id);
const canvas = $('viewer');
const viewer = new CanvasViewer(canvas);

let spec = createExample('plane');
let surface = null;
let mesh = null;
let validation = null;
let selected = [0, 0];
let referenceDirty = false;
let keyRecords = [];
let inspectorTimer = null;
let isoparametric = { fixed: 'u', value: 0.5, points: [], failures: [], visible: true, engine: null };

const display = {
  showSurface: $('showSurface').checked,
  showWire: $('showWire').checked,
  showNet: $('showNet').checked,
  showNormals: $('showNormals').checked
};

function activeReference() {
  return referenceDirty ? null : REFERENCES[spec.id] || null;
}

function pointArray(p) {
  if (Array.isArray(p)) return p;
  return [p.x, p.y, p.z];
}

function rawControlPoints() {
  const points = [];
  spec.points.forEach((row, i) => row.forEach((p, j) => {
    points.push({ index: i * spec.points[0].length + j, i, j, position: pointArray(p), weight: spec.weights[i][j] });
  }));
  return points;
}

function rawLinePairs() {
  const nv = spec.points[0].length;
  const pairs = [];
  spec.points.forEach((row, i) => row.forEach((_, j) => {
    const idx = i * nv + j;
    if (i + 1 < spec.points.length) pairs.push([idx, (i + 1) * nv + j]);
    if (j + 1 < nv) pairs.push([idx, idx + 1]);
  }));
  return pairs;
}

function currentEvalRequest() {
  if (!surface) return null;
  const [u0, u1] = surface.domainU;
  const [v0, v1] = surface.domainV;
  const fu = Number($('evalU').value) / 1000;
  const fv = Number($('evalV').value) / 1000;
  const sideUValue = $('sideU').value;
  const sideVValue = $('sideV').value;
  return {
    u: u0 + fu * (u1 - u0),
    v: v0 + fv * (v1 - v0),
    options: {
      sideU: sideUValue === 'auto' ? 'auto' : Number(sideUValue),
      sideV: sideVValue === 'auto' ? 'auto' : Number(sideVValue)
    }
  };
}

function buildScene(marker = null) {
  return {
    mesh,
    marker,
    isoparametric,
    showSurface: display.showSurface,
    showWire: display.showWire,
    showNet: display.showNet,
    showNormals: display.showNormals,
    net: {
      controlPoints: rawControlPoints(),
      linePairs: rawLinePairs(),
      selected
    }
  };
}

function renderViewer(marker = null) {
  viewer.render(buildScene(marker));
}

function renderValidation() {
  const box = $('validationBox');
  box.innerHTML = '';
  if (!validation) return;
  for (const e of validation.errors) {
    const div = document.createElement('div');
    div.className = 'msg error';
    div.textContent = `${e.code}: ${e.message}`;
    box.appendChild(div);
  }
  for (const w of validation.warnings) {
    const div = document.createElement('div');
    div.className = 'msg warning';
    div.textContent = `${w.code}: ${w.message}`;
    box.appendChild(div);
  }
  if (validation.valid) {
    const ok = document.createElement('div');
    ok.className = 'msg ok';
    const info = validation.info;
    ok.textContent = `有效：${info.nu}×${info.nv} 控制点；U∈[${formatNumber(info.domainU[0])}, ${formatNumber(info.domainU[1])}]，V∈[${formatNumber(info.domainV[0])}, ${formatNumber(info.domainV[1])}]；权重比 ${Number.isFinite(info.weightRatio) ? info.weightRatio.toExponential(2) : '∞'}`;
    box.prepend(ok);
  }
}

function renderTitle(failureCount = 0) {
  if (!surface) {
    $('viewerTitle').textContent = '当前参数无效：请查看左侧校验信息';
    return;
  }
  const segU = Number($('segU').value) || 64;
  const segV = Number($('segV').value) || 32;
  $('viewerTitle').textContent = `${spec.name} · ${surface.nu}×${surface.nv} CP · p=${surface.p}, q=${surface.q} · 三角网格 ${segU}×${segV}${failureCount ? ` · ${failureCount} 个退化采样已跳过` : ''}`;
}

function renderKernelBadge() {
  const badge = $('kernelBadge');
  if (!surface) {
    badge.textContent = '内核：无有效曲面';
    badge.classList.add('degraded');
    return;
  }
  const meshDegraded = mesh && mesh.failures.length > 0;
  badge.classList.toggle('degraded', meshDegraded);
  badge.textContent = meshDegraded
    ? `内核：${surface.engine} v${surface.kernelVersion} · ${mesh.failures.length} 退化采样`
    : `内核：${surface.engine} v${surface.kernelVersion}`;
}

function recomputeIsoparametric() {
  if (!surface) {
    isoparametric = { ...isoparametric, points: [], failures: [], engine: null };
    renderIsoReport();
    return;
  }
  const fixed = $('isoFixed').value === 'v' ? 'v' : 'u';
  const value = Number($('isoValue').value);
  const curve = surface.isoparametricCurve(fixed, value, 64);
  // Cross-check the new-kernel iso points against the preserved legacy engine.
  const crossErrors = [];
  const domain = fixed === 'v' ? surface.domainU : surface.domainV;
  curve.points.forEach((point, k) => {
    const t = domain[0] + ((domain[1] - domain[0]) * k) / Math.max(1, curve.points.length - 1);
    const legacy = fixed === 'v' ? surface.evaluateLegacy(t, value) : surface.evaluateLegacy(value, t);
    if (legacy.ok) crossErrors.push(maxComponentError(point, legacy.position));
  });
  isoparametric = {
    fixed,
    value,
    points: curve.points,
    failures: curve.failures,
    visible: $('isoVisible').checked,
    engine: curve.engine,
    crossMaxError: crossErrors.length ? Math.max(...crossErrors) : null
  };
  renderIsoReport();
}

function renderIsoReport() {
  const host = $('isoReport');
  if (!surface) {
    host.innerHTML = '<div class="bad">没有可求值的有效曲面。</div>';
    return;
  }
  const iso = isoparametric;
  if (!iso.points.length) {
    host.innerHTML = `
      <div class="bad">等参线无有效采样点。</div>
      ${iso.failures.map(f => `<div class="warn mono">${f.u?.toFixed?.(3) ?? ''} ${f.code}</div>`).join('')}
    `;
    return;
  }
  const crossCls = iso.crossMaxError === null ? 'warn' : iso.crossMaxError < 1e-9 ? 'good' : 'bad';
  const crossText = iso.crossMaxError === null ? '无对照' : iso.crossMaxError.toExponential(3);
  host.innerHTML = `
    <div class="kv"><span>驱动内核</span><span class="mono good">${iso.engine} v${surface.kernelVersion}</span></div>
    <div class="kv"><span>固定</span><span class="mono">${iso.fixed}=${formatNumber(iso.value, 4)}</span></div>
    <div class="kv"><span>有效采样</span><span class="mono">${iso.points.length}</span></div>
    <div class="kv"><span>退化采样</span><span class="mono ${iso.failures.length ? 'warn' : 'good'}">${iso.failures.length}</span></div>
    <div class="kv"><span>对旧内核最大位置误差</span><span class="mono ${crossCls}">${crossText}</span></div>
    ${iso.failures.length ? `<div class="bad">失败已回传给样条求值：${iso.failures[0].code}</div>` : '<div class="good">形态与等参线校验通过。</div>'}
  `;
}

function rebuild({ fit = false } = {}) {
  validation = validateSpec(spec);
  renderValidation();
  if (!validation.valid) {
    surface = null;
    mesh = null;
    renderTitle();
    renderKernelBadge();
    recomputeIsoparametric();
    renderViewer();
    return;
  }
  try {
    surface = compileSurface(spec);
    const segU = Math.max(1, Math.min(240, Number($('segU').value) || spec.tessellation?.uSegments || 64));
    const segV = Math.max(1, Math.min(160, Number($('segV').value) || spec.tessellation?.vSegments || 32));
    mesh = surface.tessellate(segU, segV);
    spec.tessellation = { uSegments: segU, vSegments: segV };
    if (fit) viewer.fit(spec.points.flatMap(row => row.map(pointArray)));
    renderTitle(mesh.failures.length);
    renderKernelBadge();
    recomputeIsoparametric();
    updateInspector();
  } catch (error) {
    console.error(error);
    surface = null;
    mesh = null;
    renderKernelBadge();
    $('validationBox').insertAdjacentHTML('beforeend', `<div class="msg error">${error.message}</div>`);
  }
}

function syncForm() {
  $('specName').value = spec.name || '';
  $('degreeU').value = spec.degreeU;
  $('degreeV').value = spec.degreeV;
  $('knotsU').value = `[${spec.knotsU.join(', ')}]`;
  $('knotsV').value = `[${spec.knotsV.join(', ')}]`;
  $('segU').value = spec.tessellation?.uSegments ?? 64;
  $('segV').value = spec.tessellation?.vSegments ?? 32;
  updateSelectedEditor();
  renderCPTable();
}

function updateSelectedEditor() {
  const [i, j] = selected;
  const p = spec.points[i]?.[j];
  if (!p) return;
  $('selectedLabel').textContent = `[${i},${j}]`;
  $('cpX').value = p.x;
  $('cpY').value = p.y;
  $('cpZ').value = p.z;
  $('cpW').value = spec.weights[i][j];
  $('cpXSlider').value = Math.max(-6, Math.min(6, p.x));
  $('cpYSlider').value = Math.max(-4, Math.min(4, p.y));
  $('cpZSlider').value = Math.max(-4, Math.min(4, p.z));
  $('cpWSlider').value = Math.max(0, Math.min(5, spec.weights[i][j]));
  document.querySelectorAll('.cp-row').forEach(row => row.classList.toggle('active', Number(row.dataset.i) === i && Number(row.dataset.j) === j));
}

function renderCPTable() {
  const host = $('cpTable');
  host.innerHTML = '';
  const head = document.createElement('div');
  head.className = 'cp-row head';
  head.innerHTML = '<span>i,j</span><span>X</span><span>Y</span><span>Z</span><span>w</span>';
  host.appendChild(head);

  spec.points.forEach((row, i) => row.forEach((p, j) => {
    const el = document.createElement('div');
    el.className = 'cp-row';
    el.dataset.i = i;
    el.dataset.j = j;
    if (i === selected[0] && j === selected[1]) el.classList.add('active');
    el.innerHTML = `
      <span>${i},${j}</span>
      <input data-i="${i}" data-j="${j}" data-key="x" type="number" step="0.05" value="${p.x}">
      <input data-i="${i}" data-j="${j}" data-key="y" type="number" step="0.05" value="${p.y}">
      <input data-i="${i}" data-j="${j}" data-key="z" type="number" step="0.05" value="${p.z}">
      <input data-i="${i}" data-j="${j}" data-key="w" type="number" step="0.1" value="${spec.weights[i][j]}">
    `;
    host.appendChild(el);
  }));
}

function setSelected(i, j) {
  selected = [i, j];
  updateSelectedEditor();
  renderViewer(currentMarker());
}
viewer.onSelectCP = setSelected;

function currentMarker() {
  if (!surface) return null;
  const req = currentEvalRequest();
  const r = surface.evaluate(req.u, req.v, req.options);
  if (!r.ok || !r.normal) return r.ok ? { position: r.position, normal: r.normal, du: r.derivatives.u, dv: r.derivatives.v } : { position: null };
  return { position: r.position, normal: r.normal, du: r.derivatives.u, dv: r.derivatives.v };
}

function vec(v, digits = 6) {
  return `[${v.map(x => formatNumber(x, digits)).join(', ')}]`;
}

function errorLine(label, metric) {
  if (!metric) return `<div class="kv"><span class="k">${label}</span><span class="warn">无有效对照</span></div>`;
  const cls = metric.absolute < 1e-7 ? 'good' : metric.absolute < 1e-4 ? 'warn' : 'bad';
  return `<div class="kv"><span class="k">${label}</span><span class="mono ${cls}">abs=${metric.absolute.toExponential(3)}, rel=${metric.relative.toExponential(3)}</span></div>`;
}

function updateInspector() {
  if (inspectorTimer) clearTimeout(inspectorTimer);
  inspectorTimer = setTimeout(runInspectorNow, 60);
  renderViewer(currentMarker());
}

function runInspectorNow() {
  const host = $('inspector');
  if (!surface) {
    host.innerHTML = '<div class="bad">没有可求值的有效曲面。</div>';
    return;
  }
  const req = currentEvalRequest();
  const compare = compareAtPoint(surface, req.u, req.v, activeReference(), req.options);
  const r = compare.actual;
  renderViewer(currentMarker());

  if (!r.ok) {
    host.innerHTML = `
      <div><span class="badge branch">退化分支</span><b>${r.code}</b></div>
      <div>u=${formatNumber(req.u, 8)}, v=${formatNumber(req.v, 8)}；侧 U=${r.sideU}, V=${r.sideV}</div>
      <div class="msg error">${r.message}</div>
      ${r.limitCandidates?.alongU ? `<div class="kv"><span class="k">沿 U 极限</span><span class="mono info">${vec(r.limitCandidates.alongU.position)}</span></div>` : ''}
      ${r.limitCandidates?.alongV ? `<div class="kv"><span class="k">沿 V 极限</span><span class="mono info">${vec(r.limitCandidates.alongV.position)}</span></div>` : ''}
      <div class="trace">${r.trace.map(t => `<span>${t.message}</span>`).join('')}</div>`;
    return;
  }

  host.innerHTML = `
    <div><span class="badge pass">OK</span><b>齐次有理求值</b></div>
    <div class="kv"><span class="k">参数</span><span class="mono">u=${formatNumber(r.u, 8)}, v=${formatNumber(r.v, 8)}；side=${r.sideU}/${r.sideV}</span></div>
    <div class="kv"><span class="k">位置 S</span><span class="mono good">${vec(r.position)}</span></div>
    <div class="kv"><span class="k">∂S/∂u</span><span class="mono">${vec(r.derivatives.u)}</span></div>
    <div class="kv"><span class="k">∂S/∂v</span><span class="mono">${vec(r.derivatives.v)}</span></div>
    <div class="kv"><span class="k">法线</span><span class="mono ${r.normal ? 'good' : 'bad'}">${r.normal ? vec(r.normal) : `${r.normalStatus}（未输出任意向量）`}</span></div>
    <div class="kv"><span class="k">分母条件</span><span class="mono">W=${formatNumber(r.branch.denominator, 8)}，抵消比=${r.branch.cancellationRatio.toExponential(3)}，${r.denominatorStatus}</span></div>
    ${errorLine('差分 U', compare.errors.finiteDifferenceU)}
    ${errorLine('差分 V', compare.errors.finiteDifferenceV)}
    ${errorLine('细分位置', compare.errors.subdivisionPosition)}
    ${errorLine('解析位置', compare.errors.analyticPosition)}
    ${r.normalDot !== undefined ? `<div class="kv"><span class="k">解析法线</span><span class="mono">dot=${r.normalDot.toFixed(9)}</span></div>` : ''}
    <div class="trace">${r.trace.length ? r.trace.map(t => `<span>${t.message}</span>`).join('') : '<span>参数位于光滑开区间，无额外节点分支。</span>'}</div>
  `;
}

function updateSelectedFromInputs(source) {
  const [i, j] = selected;
  if (source === 'selected') {
    for (const key of ['x', 'y', 'z']) {
      const value = Number($(`cp${key.toUpperCase()}`).value);
      if (Number.isFinite(value)) spec.points[i][j][key] = value;
    }
    const w = Number($('cpW').value);
    if (Number.isFinite(w)) spec.weights[i][j] = w;
  }
  referenceDirty = true;
  rebuild();
  updateSelectedEditor();
  document.querySelectorAll(`[data-i="${i}"][data-j="${j}"]`).forEach(input => {
    const key = input.dataset.key;
    input.value = key === 'w' ? spec.weights[i][j] : spec.points[i][j][key];
  });
}

$('cpX').addEventListener('input', () => updateSelectedFromInputs('selected'));
$('cpY').addEventListener('input', () => updateSelectedFromInputs('selected'));
$('cpZ').addEventListener('input', () => updateSelectedFromInputs('selected'));
$('cpW').addEventListener('input', () => updateSelectedFromInputs('selected'));
for (const id of ['cpXSlider', 'cpYSlider', 'cpZSlider', 'cpWSlider']) {
  $(id).addEventListener('input', () => {
    const [i, j] = selected;
    const key = id === 'cpWSlider' ? 'w' : id[2].toLowerCase();
    const value = Number($(id).value);
    if (key === 'w') {
      spec.weights[i][j] = value;
      $('cpW').value = value;
    } else {
      spec.points[i][j][key] = value;
      $(`cp${key.toUpperCase()}`).value = value;
    }
    referenceDirty = true;
    rebuild();
    const tableInput = document.querySelector(`#cpTable [data-i="${i}"][data-j="${j}"][data-key="${key}"]`);
    if (tableInput) tableInput.value = value;
  });
}

$('cpTable').addEventListener('input', event => {
  const input = event.target;
  if (!input.dataset.key) return;
  const i = Number(input.dataset.i);
  const j = Number(input.dataset.j);
  const key = input.dataset.key;
  let value = parseFloat(input.value);
  if (!Number.isFinite(value)) value = 0;
  if (key === 'w') spec.weights[i][j] = value;
  else spec.points[i][j][key] = value;
  selected = [i, j];
  referenceDirty = true;
  rebuild();
  updateSelectedEditor();
});
$('cpTable').addEventListener('click', event => {
  const row = event.target.closest('.cp-row:not(.head)');
  if (row) setSelected(Number(row.dataset.i), Number(row.dataset.j));
});

function parseKnots(text) {
  return text.replace(/[\[\]()]/g, ' ').split(/[\s,]+/).filter(Boolean).map(Number);
}

$('applyTopology').addEventListener('click', () => {
  spec.name = $('specName').value;
  spec.degreeU = Number($('degreeU').value);
  spec.degreeV = Number($('degreeV').value);
  spec.knotsU = parseKnots($('knotsU').value);
  spec.knotsV = parseKnots($('knotsV').value);
  referenceDirty = true;
  rebuild();
  syncForm();
});
$('specName').addEventListener('change', () => { spec.name = $('specName').value; renderTitle(); });

for (const id of ['segU', 'segV']) {
  $(id).addEventListener('change', () => rebuild());
}
for (const id of ['showSurface', 'showWire', 'showNet', 'showNormals']) {
  $(id).addEventListener('change', () => {
    display[id] = $(id).checked;
    renderViewer(currentMarker());
  });
}
for (const id of ['evalU', 'evalV', 'sideU', 'sideV']) {
  $(id).addEventListener('input', updateInspector);
}
for (const id of ['isoFixed', 'isoValue', 'isoVisible']) {
  $(id).addEventListener('input', () => {
    if (id === 'isoVisible') isoparametric.visible = $('isoVisible').checked;
    recomputeIsoparametric();
    renderViewer(currentMarker());
  });
}

const select = $('exampleSelect');
for (const id of EXAMPLE_IDS) {
  const ex = createExample(id);
  const option = document.createElement('option');
  option.value = id;
  option.textContent = ex.name;
  select.appendChild(option);
}
select.addEventListener('change', () => loadSpec(createExample(select.value), true, true));

function loadSpec(next, fit = true, knownReference = false) {
  spec = next;
  referenceDirty = !knownReference;
  selected = [0, 0];
  syncForm();
  rebuild({ fit });
  updateInspector();
}

$('saveProject').addEventListener('click', () => {
  const project = {
    format: 'rational-spline-studio-project',
    version: 1,
    savedAt: new Date().toISOString(),
    spec
  };
  downloadBlob(`${spec.id || 'nurbs-project'}.json`, JSON.stringify(project, null, 2), 'application/json');
});

$('openProject').addEventListener('click', () => $('projectFile').click());
$('projectFile').addEventListener('change', async event => {
  const file = event.target.files?.[0];
  event.target.value = '';
  if (!file) return;
  try {
    const json = JSON.parse(await file.text());
    const raw = json.spec ?? json;
    raw.points = raw.points.map(row => row.map(pointArray).map(p => ({ x: p[0], y: p[1], z: p[2] })));
    if (!raw.tessellation) raw.tessellation = { uSegments: 64, vSegments: 32 };
    loadSpec(raw, true, false);
  } catch (error) {
    alert(`打开工程失败：${error.message}`);
  }
});

$('exportOBJ').addEventListener('click', () => {
  if (!mesh) return alert('当前没有有效网格可导出。');
  downloadBlob(`${spec.id || 'surface'}.obj`, exportOBJ(mesh, spec.id || 'surface'), 'text/plain');
});
$('exportSTL').addEventListener('click', () => {
  if (!mesh) return alert('当前没有有效网格可导出。');
  downloadBlob(`${spec.id || 'surface'}.stl`, exportSTL(mesh), 'model/stl');
});

function histogram(values) {
  if (!values?.histogram?.count) return '<div class="hint">无数据</div>';
  const max = Math.max(...values.histogram, 1);
  return `<div class="hist">${values.histogram.map(n => `<i style="height:${Math.max(2, (n / max) * 100)}%" title="${n}"></i>`).join('')}</div>`;
}

function metricCard(label, stat) {
  if (!stat || stat.count === 0) return `<div class="stat"><b>${label}</b><span>无有效样本</span></div>`;
  return `<div class="stat"><b>${label}</b><span>max ${stat.max.toExponential(2)}<br>RMS ${stat.rms.toExponential(2)}<br>P95 ${stat.p95.toExponential(2)}</span></div>`;
}

function renderDistribution(report) {
  const labels = {
    analyticPosition: '解析位置 abs',
    analyticDerivative: '解析导数 abs',
    finiteDifferenceDerivative: '差分导数 rel',
    subdivisionPosition: '细分位置 abs',
    subdivisionDerivative: '细分导数 rel',
    basisPartition: '单位分解 abs'
  };
  $('distribution').innerHTML = `
    <div class="stat-grid">
      ${Object.entries(labels).map(([key, label]) => metricCard(label, report.stats[key])).join('')}
    </div>
    <h2 style="margin-top:8px">差分导数误差直方图</h2>
    ${histogram(report.stats.finiteDifferenceDerivative)}
    <h2 style="margin-top:8px">误差热图（绿=优，红=差）</h2>
    ${renderHeatmap(report.heatmap)}
  `;
}

function renderHeatmap(heatmap) {
  const n = heatmap.length;
  const cells = heatmap.flat().map(e => {
    let background = '#334155';
    if (Number.isFinite(e)) {
      const quality = Math.max(0, Math.min(1, (-Math.log10(Math.max(e, 1e-14)) - 5) / 4));
      const r = Math.round(220 * (1 - quality) + 34 * quality);
      const g = Math.round(70 * (1 - quality) + 197 * quality);
      const b = Math.round(70 * (1 - quality) + 94 * quality);
      background = `rgb(${r},${g},${b})`;
    }
    return `<div class="cell" title="${e === null ? '退化/无效' : e.toExponential(4)}" style="background:${background}"></div>`;
  }).join('');
  return `<div class="heatmap" style="grid-template-columns:repeat(${n},1fr)">${cells}</div>`;
}

function renderKeyChecks(records) {
  keyRecords = records;
  const host = $('keyChecks');
  const counts = records.reduce((acc, r) => { acc[r.status]++; return acc; }, { pass: 0, fail: 0, branch: 0 });
  host.innerHTML = `<div class="hint">通过 ${counts.pass}，分支 ${counts.branch}，失败 ${counts.fail}。点击“定位”可载入检查器。</div>`;
  records.forEach((r, index) => {
    const el = document.createElement('details');
    el.className = `check-record ${r.status}`;
    el.innerHTML = `
      <summary><span class="badge ${r.status}">${r.status}</span>${r.label}
        <button data-open="${index}" style="float:right;padding:1px 6px">定位</button>
      </summary>
      ${r.checks.map(c => `<div class="check-line ${c.pass ? 'pass' : 'fail'}">${c.pass ? '✓' : '×'} ${c.name}${c.detail ? ` — <span class="mono">${c.detail}</span>` : ''}</div>`).join('')}
    `;
    host.appendChild(el);
  });
  host.querySelectorAll('button[data-open]').forEach(button => button.addEventListener('click', event => {
    event.preventDefault();
    event.stopPropagation();
    const r = keyRecords[Number(button.dataset.open)];
    const [u0, u1] = surface.domainU;
    const [v0, v1] = surface.domainV;
    $('evalU').value = Math.round(((r.u - u0) / (u1 - u0)) * 1000);
    $('evalV').value = Math.round(((r.v - v0) / (v1 - v0)) * 1000);
    $('sideU').value = r.sideU === -1 ? '-1' : r.sideU === 1 ? '1' : 'auto';
    $('sideV').value = r.sideV === -1 ? '-1' : r.sideV === 1 ? '1' : 'auto';
    updateInspector();
  }));
}

async function runChecks() {
  if (!surface) return;
  const button = $('runVerification');
  button.disabled = true;
  button.textContent = '检查中…';
  await new Promise(resolve => setTimeout(resolve, 30));
  try {
    const reference = activeReference();
    const report = runErrorDistribution(surface, reference, 18);
    const records = runKeyChecks(surface, reference);
    renderDistribution(report);
    renderKeyChecks(records);
  } finally {
    button.disabled = false;
    button.textContent = '运行独立检查';
  }
}
$('runVerification').addEventListener('click', runChecks);

new ResizeObserver(() => viewer.resize()).observe($('viewer').parentElement);
window.addEventListener('resize', () => viewer.resize());

loadSpec(spec, true, true);
setTimeout(runChecks, 250);
