import { cross3, sub3, norm3 } from './math-utils.js';

export class CanvasViewer {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.camera = { yaw: 0.72, pitch: 0.42, distance: 9, target: [0, 0, 0] };
    this.onSelectCP = null;
    this.drag = null;
    this.clickPoint = null;
    this.#bindEvents();
    this.resize();
  }

  #bindEvents() {
    const c = this.canvas;
    c.addEventListener('pointerdown', e => {
      c.setPointerCapture(e.pointerId);
      this.drag = { x: e.clientX, y: e.clientY, startX: e.clientX, startY: e.clientY, button: e.button, shift: e.shiftKey };
    });
    c.addEventListener('pointermove', e => {
      if (!this.drag) return;
      const dx = e.clientX - this.drag.x;
      const dy = e.clientY - this.drag.y;
      if (this.drag.shift || this.drag.button === 2) {
        const scale = this.camera.distance * 0.0015;
        const right = [Math.cos(this.camera.yaw), 0, -Math.sin(this.camera.yaw)];
        this.camera.target[0] -= right[0] * dx * scale;
        this.camera.target[2] -= right[2] * dx * scale;
        this.camera.target[1] += dy * scale;
      } else {
        this.camera.yaw -= dx * 0.008;
        this.camera.pitch = Math.max(-1.45, Math.min(1.45, this.camera.pitch + dy * 0.008));
      }
      this.drag.x = e.clientX;
      this.drag.y = e.clientY;
      this.render();
    });
    c.addEventListener('pointerup', e => {
      const moved = Math.hypot(e.clientX - this.drag.startX, e.clientY - this.drag.startY);
      if (moved < 5 && this.onSelectCP) {
        const rect = c.getBoundingClientRect();
        this.clickPoint = { x: e.clientX - rect.left, y: e.clientY - rect.top };
        const picked = this.#pickControlPoint();
        if (picked) this.onSelectCP(picked.i, picked.j);
      }
      this.drag = null;
    });
    c.addEventListener('wheel', e => {
      e.preventDefault();
      this.camera.distance *= Math.exp(e.deltaY * 0.001);
      this.camera.distance = Math.max(0.05, Math.min(1000, this.camera.distance));
      this.render();
    }, { passive: false });
    c.addEventListener('contextmenu', e => e.preventDefault());
  }

  resize() {
    const rect = this.canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    this.width = Math.max(1, rect.width);
    this.height = Math.max(1, rect.height);
    this.canvas.width = Math.floor(this.width * dpr);
    this.canvas.height = Math.floor(this.height * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.render();
  }

  fit(points) {
    if (!points?.length) return;
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (const p of points) {
      for (let k = 0; k < 3; k++) {
        min[k] = Math.min(min[k], p[k]);
        max[k] = Math.max(max[k], p[k]);
      }
    }
    const center = min.map((x, k) => (x + max[k]) / 2);
    const radius = Math.max(...max.map((x, k) => Math.max(Math.abs(x - center[k]), Math.abs(min[k] - center[k]))));
    this.camera.target = center;
    this.camera.distance = Math.max(1.5, radius * 3.4);
    this.render();
  }

  project(p) {
    const [x0, y0, z0] = [p[0] - this.camera.target[0], p[1] - this.camera.target[1], p[2] - this.camera.target[2]];
    const cy = Math.cos(this.camera.yaw);
    const sy = Math.sin(this.camera.yaw);
    const cp = Math.cos(this.camera.pitch);
    const sp = Math.sin(this.camera.pitch);
    const x1 = x0 * cy + z0 * sy;
    const z1 = -x0 * sy + z0 * cy;
    const y2 = y0 * cp - z1 * sp;
    const z2 = y0 * sp + z1 * cp;
    const depth = this.camera.distance - z2;
    if (depth <= 0.01) return null;
    const focal = this.height / (2 * Math.tan(Math.PI / 8));
    const scale = focal / depth;
    return {
      x: this.width / 2 + x1 * scale,
      y: this.height / 2 - y2 * scale,
      depth,
      scale
    };
  }

  #pickControlPoint() {
    const controlPoints = this.scene?.net?.controlPoints;
    if (!controlPoints || !this.clickPoint) return null;
    let best = null;
    let bestDistance = 14;
    controlPoints.forEach(p => {
      const q = this.project(p.position);
      if (!q) return;
      const d = Math.hypot(q.x - this.clickPoint.x, q.y - this.clickPoint.y);
      if (d < bestDistance) {
        bestDistance = d;
        best = { i: p.i, j: p.j };
      }
    });
    return best;
  }

  render(scene = this.scene) {
    if (scene) this.scene = scene;
    const ctx = this.ctx;
    ctx.clearRect(0, 0, this.width, this.height);
    const gradient = ctx.createLinearGradient(0, 0, 0, this.height);
    gradient.addColorStop(0, '#0b1220');
    gradient.addColorStop(1, '#111827');
    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, this.width, this.height);
    this.#drawAxes();

    if (!scene) return;
    if (scene.mesh && scene.showSurface !== false) this.#drawSurface(scene.mesh, scene.showWire !== false);
    if (scene.mesh && scene.showWire && scene.showSurface === false) this.#drawGrid(scene.mesh, 'rgba(148,197,255,0.72)');
    if (scene.net && scene.showNet !== false) this.#drawControlNet(scene.net);
    if (scene.mesh && scene.showNormals) this.#drawMeshNormals(scene.mesh);
    if (scene.marker) this.#drawMarker(scene.marker);
  }

  #drawAxes() {
    const origin = this.project(this.camera.target);
    if (!origin) return;
    const axes = [
      { v: [1, 0, 0], color: '#ef4444', label: 'X' },
      { v: [0, 1, 0], color: '#22c55e', label: 'Y' },
      { v: [0, 0, 1], color: '#3b82f6', label: 'Z' }
    ];
    const ctx = this.ctx;
    ctx.font = '12px system-ui';
    for (const ax of axes) {
      const end = this.project(this.camera.target.map((x, k) => x + ax.v[k] * 0.8));
      if (!end) continue;
      ctx.strokeStyle = ax.color;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(origin.x, origin.y);
      ctx.lineTo(end.x, end.y);
      ctx.stroke();
      ctx.fillStyle = ax.color;
      ctx.fillText(ax.label, end.x + 4, end.y - 4);
    }
  }

  #validPoint(mesh, index) {
    const k = index * 3;
    return Number.isFinite(mesh.positions[k]);
  }

  #point(mesh, index) {
    const k = index * 3;
    return [mesh.positions[k], mesh.positions[k + 1], mesh.positions[k + 2]];
  }

  #drawSurface(mesh, wire) {
    const triangles = [];
    for (let t = 0; t < mesh.indices.length; t += 3) {
      const ia = mesh.indices[t];
      const ib = mesh.indices[t + 1];
      const ic = mesh.indices[t + 2];
      if (!this.#validPoint(mesh, ia) || !this.#validPoint(mesh, ib) || !this.#validPoint(mesh, ic)) continue;
      const a = this.#point(mesh, ia);
      const b = this.#point(mesh, ib);
      const c = this.#point(mesh, ic);
      const pa = this.project(a);
      const pb = this.project(b);
      const pc = this.project(c);
      if (!pa || !pb || !pc) continue;
      const n = cross3(sub3(b, a), sub3(c, a));
      const nl = norm3(n) || 1;
      const nz = n[2] / nl;
      const ny = n[1] / nl;
      const light = Math.max(0.28, 0.45 + 0.35 * Math.abs(nz) + 0.2 * Math.max(0, ny));
      triangles.push({ pa, pb, pc, depth: (pa.depth + pb.depth + pc.depth) / 3, light });
    }
    triangles.sort((a, b) => b.depth - a.depth);
    const ctx = this.ctx;
    for (const tr of triangles) {
      const l = Math.round(110 + tr.light * 95);
      ctx.fillStyle = `rgb(${Math.round(35 * tr.light + 20)},${Math.round(105 * tr.light + 45)},${Math.min(255, l + 35)})`;
      ctx.beginPath();
      ctx.moveTo(tr.pa.x, tr.pa.y);
      ctx.lineTo(tr.pb.x, tr.pb.y);
      ctx.lineTo(tr.pc.x, tr.pc.y);
      ctx.closePath();
      ctx.fill();
    }
    if (wire) this.#drawGrid(mesh, 'rgba(219,234,254,0.34)');
  }

  #drawGrid(mesh, color) {
    const ctx = this.ctx;
    const cols = mesh.cols;
    const rows = mesh.rows;
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let j = 0; j < rows; j++) {
      for (let i = 0; i < cols; i++) {
        const a = j * cols + i;
        if (i + 1 < cols) this.#segmentPath(mesh, a, a + 1);
        if (j + 1 < rows) this.#segmentPath(mesh, a, a + cols);
      }
    }
    ctx.stroke();
  }

  #segmentPath(mesh, ia, ib) {
    if (!this.#validPoint(mesh, ia) || !this.#validPoint(mesh, ib)) return;
    const a = this.project(this.#point(mesh, ia));
    const b = this.project(this.#point(mesh, ib));
    if (!a || !b) return;
    this.ctx.moveTo(a.x, a.y);
    this.ctx.lineTo(b.x, b.y);
  }

  #drawControlNet(net) {
    const ctx = this.ctx;
    const projected = net.controlPoints.map(p => ({ ...p, q: this.project(p.position) })).filter(p => p.q);
    ctx.strokeStyle = 'rgba(251,146,60,0.72)';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (const [a, b] of net.linePairs) {
      const pa = projected.find(p => p.index === a)?.q;
      const pb = projected.find(p => p.index === b)?.q;
      if (pa && pb) {
        ctx.moveTo(pa.x, pa.y);
        ctx.lineTo(pb.x, pb.y);
      }
    }
    ctx.stroke();

    for (const p of projected) {
      const selected = p.i === net.selected?.[0] && p.j === net.selected?.[1];
      const weight = p.weight ?? 1;
      const radius = selected ? 7 : Math.max(3, Math.min(7, 3 + Math.log10(Math.max(1, weight)) * 0.7));
      ctx.beginPath();
      ctx.arc(p.q.x, p.q.y, radius, 0, Math.PI * 2);
      ctx.fillStyle = selected ? '#fde047' : weight === 0 ? '#f87171' : '#fb923c';
      ctx.fill();
      ctx.strokeStyle = '#431407';
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }

  #drawMeshNormals(mesh) {
    const stepI = Math.max(1, Math.floor(mesh.cols / 12));
    const stepJ = Math.max(1, Math.floor(mesh.rows / 8));
    const ctx = this.ctx;
    ctx.strokeStyle = 'rgba(34,211,238,0.72)';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let j = 0; j < mesh.rows; j += stepJ) {
      for (let i = 0; i < mesh.cols; i += stepI) {
        const idx = j * mesh.cols + i;
        if (!mesh.valid[idx]) continue;
        const p = this.#point(mesh, idx);
        const n = [mesh.normals[idx * 3], mesh.normals[idx * 3 + 1], mesh.normals[idx * 3 + 2]];
        const length = this.camera.distance * 0.045;
        const a = this.project(p);
        const b = this.project(p.map((x, k) => x + n[k] * length));
        if (a && b) {
          ctx.moveTo(a.x, a.y);
          ctx.lineTo(b.x, b.y);
        }
      }
    }
    ctx.stroke();
  }

  #drawMarker(marker) {
    if (!marker?.position || !Number.isFinite(marker.position[0])) return;
    const p = this.project(marker.position);
    if (!p) return;
    const ctx = this.ctx;
    const length = this.camera.distance * 0.055;
    const vectors = [
      { v: marker.normal, color: '#22d3ee', label: 'n' },
      { v: marker.du, color: '#f87171', label: 'u' },
      { v: marker.dv, color: '#60a5fa', label: 'v' }
    ];
    ctx.font = 'bold 12px system-ui';
    for (const item of vectors) {
      if (!item.v) continue;
      const end = this.project(marker.position.map((x, k) => x + item.v[k] * length / Math.max(norm3(item.v), 1e-12) * (item.label === 'n' ? 1 : 0.7)));
      if (!end) continue;
      ctx.strokeStyle = item.color;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.moveTo(p.x, p.y);
      ctx.lineTo(end.x, end.y);
      ctx.stroke();
      ctx.fillStyle = item.color;
      ctx.fillText(item.label, end.x + 4, end.y - 4);
    }
    ctx.beginPath();
    ctx.arc(p.x, p.y, 6, 0, Math.PI * 2);
    ctx.fillStyle = '#fef08a';
    ctx.fill();
    ctx.strokeStyle = '#713f12';
    ctx.stroke();
  }
}
