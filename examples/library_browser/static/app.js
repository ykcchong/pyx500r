"use strict";

// ── tiny fetch helper ──────────────────────────────────────────────
async function getJSON(path) {
  const res = await fetch(path);
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status}: ${body}`);
  }
  return res.json();
}

// ── state ───────────────────────────────────────────────────────────
const state = {
  search: "",
  libraryId: "",
  offset: 0,
  limit: 100,
  total: 0,
  selectedCompound: null,
  spectra: [],
  selectedSpectrumId: null,
  kind: "centroid",
  spectrumCache: new Map(), // `${id}:${kind}` -> spectrum
};

// ── elements ────────────────────────────────────────────────────────
const el = {
  stats: document.getElementById("stats"),
  search: document.getElementById("search"),
  library: document.getElementById("library"),
  list: document.getElementById("compound-list"),
  pager: document.getElementById("pager"),
  detail: document.getElementById("detail"),
  tabs: document.getElementById("spectrum-tabs"),
  plotControls: document.getElementById("plot-controls"),
  plotInfo: document.getElementById("plot-info"),
  plot: document.getElementById("plot"),
  spectrumMeta: document.getElementById("spectrum-meta"),
};

// ── init ────────────────────────────────────────────────────────────
async function init() {
  try {
    const stats = await getJSON("/api/stats");
    el.stats.textContent = `${stats.compounds.toLocaleString()} compounds · ${stats.spectra.toLocaleString()} spectra · ${stats.libraries} libraries`;
  } catch (e) { el.stats.textContent = String(e); }

  try {
    const libs = await getJSON("/api/libraries");
    for (const lib of libs) {
      const opt = document.createElement("option");
      opt.value = lib.id;
      opt.textContent = `${lib.name} (${lib.compound_count})`;
      el.library.appendChild(opt);
    }
  } catch (e) { /* non-fatal */ }

  el.search.addEventListener("input", debounce(() => {
    state.search = el.search.value.trim();
    state.offset = 0;
    loadCompounds();
  }, 250));

  el.library.addEventListener("change", () => {
    state.libraryId = el.library.value;
    state.offset = 0;
    loadCompounds();
  });

  el.plotControls.addEventListener("change", (ev) => {
    if (ev.target.name === "kind") {
      state.kind = ev.target.value;
      if (state.selectedSpectrumId) showSpectrum(state.selectedSpectrumId);
    }
  });

  window.addEventListener("resize", () => {
    const sp = currentSpectrum();
    if (sp) drawSpectrum(sp);
  });

  loadCompounds();
}

// ── compound list ───────────────────────────────────────────────────
async function loadCompounds() {
  el.list.innerHTML = `<div class="empty-note">Loading…</div>`;
  const params = new URLSearchParams({
    offset: state.offset, limit: state.limit,
  });
  if (state.search) params.set("search", state.search);
  if (state.libraryId) params.set("library_id", state.libraryId);

  let data;
  try {
    data = await getJSON(`/api/compounds?${params}`);
  } catch (e) {
    el.list.innerHTML = `<div class="empty-note">${e}</div>`;
    return;
  }
  state.total = data.total;
  renderCompoundList(data.items);
  renderPager();
}

function renderCompoundList(items) {
  if (!items.length) {
    el.list.innerHTML = `<div class="empty-note">No compounds match.</div>`;
    return;
  }
  el.list.innerHTML = "";
  for (const c of items) {
    const row = document.createElement("div");
    row.className = "compound-row";
    if (c.id === state.selectedCompound) row.classList.add("active");
    row.innerHTML = `
      <div class="name">${escapeHtml(c.name)}</div>
      <div class="sub">
        <span>${escapeHtml(c.formula || "—")}</span>
        ${c.cas ? `<span>CAS ${escapeHtml(c.cas)}</span>` : ""}
        <span class="badge">${c.spectrum_count} spec</span>
      </div>`;
    row.addEventListener("click", () => selectCompound(c.id, row));
    el.list.appendChild(row);
  }
}

function renderPager() {
  const start = state.total ? state.offset + 1 : 0;
  const end = Math.min(state.offset + state.limit, state.total);
  el.pager.innerHTML = "";
  const prev = button("‹ Prev", state.offset === 0, () => {
    state.offset = Math.max(0, state.offset - state.limit); loadCompounds();
  });
  const count = document.createElement("span");
  count.className = "count";
  count.textContent = `${start}–${end} of ${state.total.toLocaleString()}`;
  const next = button("Next ›", end >= state.total, () => {
    state.offset += state.limit; loadCompounds();
  });
  el.pager.append(prev, count, next);
}

// ── compound detail ─────────────────────────────────────────────────
async function selectCompound(id, rowEl) {
  state.selectedCompound = id;
  document.querySelectorAll(".compound-row.active").forEach(r => r.classList.remove("active"));
  if (rowEl) rowEl.classList.add("active");

  el.detail.innerHTML = `<div class="empty-note">Loading…</div>`;
  let c;
  try {
    c = await getJSON(`/api/compounds/${id}`);
  } catch (e) {
    el.detail.innerHTML = `<div class="empty-note">${e}</div>`;
    return;
  }
  renderDetail(c);
  renderSpectrumTabs(c.spectra);
}

function renderDetail(c) {
  const names = c.names.map(n =>
    `<span class="chip ${n.is_default ? "default" : ""}">${escapeHtml(n.name)}${n.region ? ` · ${escapeHtml(n.region)}` : ""}</span>`
  ).join("");
  const libs = c.libraries.map(l => `<span class="chip">${escapeHtml(l.name)}</span>`).join("") || `<span class="muted">—</span>`;
  const rts = c.retention_times.length
    ? c.retention_times.map(rt => `<span class="chip ${rt.is_default ? "default" : ""}">${fmt(rt.value, 2)} min${rt.instrument ? ` · ${escapeHtml(rt.instrument)}` : ""}</span>`).join("")
    : `<span class="muted">—</span>`;

  el.detail.innerHTML = `
    <h2 class="detail-title">${escapeHtml(c.name)}</h2>
    <p class="detail-sub">${escapeHtml(c.formula || "")} ${c.cas ? "· CAS " + escapeHtml(c.cas) : ""}</p>

    <dl class="kv">
      <dt>Identifier</dt><dd>${escapeHtml(c.identifier || "—")}</dd>
      <dt>Molecular weight</dt><dd>${fmt(c.molecular_weight, 4)}</dd>
      <dt>Monoisotopic mass</dt><dd>${fmt(c.monoisotopic_mass, 4)}</dd>
      <dt>Structure source</dt><dd>${escapeHtml(c.molecular_structure_source || "—")}</dd>
      <dt>Active</dt><dd>${c.active === null ? "—" : c.active}</dd>
      <dt>Last updated</dt><dd>${escapeHtml(c.last_updated || "—")}</dd>
    </dl>

    <div class="threshold-bar">
      <div class="t"><div class="v">${fmt(c.purity_threshold, 0, "—")}</div><div class="l">Purity</div></div>
      <div class="t"><div class="v" style="color:var(--warn)">${fmt(c.yellow_flag_threshold, 0, "—")}</div><div class="l">Yellow flag</div></div>
      <div class="t"><div class="v" style="color:var(--neg)">${fmt(c.red_flag_threshold, 0, "—")}</div><div class="l">Red flag</div></div>
    </div>

    <div class="section-h">Names (${c.names.length})</div>
    <div>${names}</div>

    <div class="section-h">Libraries (${c.libraries.length})</div>
    <div>${libs}</div>

    <div class="section-h">Retention times (${c.retention_times.length})</div>
    <div>${rts}</div>

    ${c.comment ? `<div class="section-h">Comment</div><div>${escapeHtml(c.comment)}</div>` : ""}
  `;
}

// ── spectra tabs + plot ─────────────────────────────────────────────
function renderSpectrumTabs(spectra) {
  state.spectra = spectra || [];
  state.selectedSpectrumId = null;
  el.tabs.innerHTML = "";
  clearPlot();

  if (!state.spectra.length) {
    el.tabs.innerHTML = `<div class="empty-note">No reference spectra.</div>`;
    el.plotControls.classList.add("hidden");
    el.spectrumMeta.innerHTML = "";
    return;
  }

  state.spectra.forEach((s, i) => {
    const tab = document.createElement("div");
    const decodable = s.has_centroid || s.has_raw;
    tab.className = "tab" + (s.encrypted || !decodable ? " empty" : "");
    tab.innerHTML = `
      <span class="pol ${s.polarity}">${s.polarity}</span>
      <span class="meta">CE ${fmt(s.collision_energy, 0, "?")} · ${escapeHtml(s.instrument || "?")}</span>
      <span class="meta">m/z ${fmt(s.precursor_mz, 4, "—")}</span>`;
    tab.title = s.encrypted ? "This spectrum's data is encrypted and cannot be decoded." : "";
    tab.addEventListener("click", () => {
      document.querySelectorAll(".tab.active").forEach(t => t.classList.remove("active"));
      tab.classList.add("active");
      showSpectrum(s.id);
    });
    el.tabs.appendChild(tab);
    if (i === 0) { tab.classList.add("active"); }
  });

  // auto-select the first spectrum that actually has decodable data
  const firstGood = state.spectra.find(s => !s.encrypted && (s.has_centroid || s.has_raw)) || state.spectra[0];
  const idx = state.spectra.indexOf(firstGood);
  document.querySelectorAll(".tab").forEach((t, i) => t.classList.toggle("active", i === idx));
  el.plotControls.classList.remove("hidden");
  showSpectrum(firstGood.id);
}

function currentSpectrum() {
  return state.spectrumCache.get(`${state.selectedSpectrumId}:${state.kind}`);
}

async function showSpectrum(id) {
  state.selectedSpectrumId = id;
  const cacheKey = `${id}:${state.kind}`;
  let sp = state.spectrumCache.get(cacheKey);
  if (!sp) {
    clearPlot("Loading…");
    try {
      sp = await getJSON(`/api/spectra/${id}?kind=${state.kind}`);
    } catch (e) {
      clearPlot(String(e));
      return;
    }
    state.spectrumCache.set(cacheKey, sp);
  }
  renderSpectrumMeta(sp);
  drawSpectrum(sp);
}

function renderSpectrumMeta(sp) {
  el.plotInfo.textContent = `${sp.num_peaks.toLocaleString()} peaks · ${sp.kind}`;
  const m = (label, val) => `<div><b>${label}:</b> ${val}</div>`;
  el.spectrumMeta.innerHTML =
    m("Polarity", sp.polarity) +
    m("Precursor m/z", fmt(sp.precursor_mz, 4, "—")) +
    m("Charge", fmt(sp.charge_state, 0, "—")) +
    m("Collision energy", fmt(sp.collision_energy, 1, "—")) +
    m("CE spread", fmt(sp.collision_energy_spread, 1, "—")) +
    m("Scan type", sp.scan_type || sp.type || "—") +
    m("Instrument", sp.instrument || "—") +
    m("Ion source", sp.ion_source || "—") +
    m("RT window", `${fmt(sp.start_rt, 2, "—")}–${fmt(sp.end_rt, 2, "—")}`);
}

// ── canvas spectrum plot (stick / centroid + profile line) ──────────
function clearPlot(message) {
  const ctx = el.plot.getContext("2d");
  resizeCanvas();
  ctx.clearRect(0, 0, el.plot.width, el.plot.height);
  if (message) {
    ctx.fillStyle = "#8593a3";
    ctx.font = "14px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(message, el.plot.width / 2 / DPR, el.plot.height / 2 / DPR);
  }
}

let DPR = window.devicePixelRatio || 1;
function resizeCanvas() {
  DPR = window.devicePixelRatio || 1;
  const rect = el.plot.getBoundingClientRect();
  el.plot.width = Math.max(300, rect.width) * DPR;
  el.plot.height = Math.max(200, rect.height) * DPR;
}

function drawSpectrum(sp) {
  resizeCanvas();
  const ctx = el.plot.getContext("2d");
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  const W = el.plot.width / DPR, H = el.plot.height / DPR;
  ctx.clearRect(0, 0, W, H);

  const mz = sp.mz, inten = sp.intensity;
  if (!mz || !mz.length) {
    ctx.fillStyle = "#8593a3";
    ctx.font = "14px sans-serif"; ctx.textAlign = "center";
    ctx.fillText(sp.encrypted ? "Spectrum data is encrypted (cannot decode)."
                              : "No peak data for this spectrum.", W / 2, H / 2);
    return;
  }

  const pad = { l: 60, r: 14, t: 14, b: 34 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;

  let mzMin = Infinity, mzMax = -Infinity, iMax = -Infinity;
  for (let i = 0; i < mz.length; i++) {
    if (mz[i] < mzMin) mzMin = mz[i];
    if (mz[i] > mzMax) mzMax = mz[i];
    if (inten[i] > iMax) iMax = inten[i];
  }
  // pad m/z range a touch
  const span = (mzMax - mzMin) || 1;
  mzMin = Math.max(0, mzMin - span * 0.02);
  mzMax = mzMax + span * 0.02;
  if (iMax <= 0) iMax = 1;

  const xToPx = (m) => pad.l + ((m - mzMin) / (mzMax - mzMin)) * plotW;
  const yToPx = (v) => pad.t + plotH - (v / iMax) * plotH;

  // axes
  ctx.strokeStyle = "#2d3845"; ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + plotH); ctx.lineTo(pad.l + plotW, pad.t + plotH);
  ctx.stroke();

  // gridlines + y labels (relative %)
  ctx.fillStyle = "#8593a3"; ctx.font = "10px sans-serif"; ctx.textAlign = "right";
  for (let f = 0; f <= 1.0001; f += 0.25) {
    const y = pad.t + plotH - f * plotH;
    ctx.strokeStyle = "#1d2630"; ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + plotW, y); ctx.stroke();
    ctx.fillText(`${Math.round(f * 100)}%`, pad.l - 6, y + 3);
  }

  // x labels
  ctx.textAlign = "center";
  const ticks = 6;
  for (let t = 0; t <= ticks; t++) {
    const m = mzMin + (mzMax - mzMin) * (t / ticks);
    const x = xToPx(m);
    ctx.fillText(m.toFixed(1), x, pad.t + plotH + 16);
  }
  ctx.fillText("m/z", pad.l + plotW / 2, H - 4);

  const isProfile = sp.kind === "raw" && mz.length > 600;

  if (isProfile) {
    // continuous line for dense profile data
    ctx.strokeStyle = "#4ea1ff"; ctx.lineWidth = 1; ctx.beginPath();
    for (let i = 0; i < mz.length; i++) {
      const x = xToPx(mz[i]), y = yToPx(inten[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  } else {
    // centroid sticks
    ctx.strokeStyle = sp.polarity === "NEG" ? "#ff7b72" : "#4ea1ff";
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 0; i < mz.length; i++) {
      const x = xToPx(mz[i]);
      ctx.moveTo(x, pad.t + plotH);
      ctx.lineTo(x, yToPx(inten[i]));
    }
    ctx.stroke();

    // label the top few peaks
    const idx = inten.map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]).slice(0, 8);
    ctx.fillStyle = "#d7dee6"; ctx.font = "10px sans-serif"; ctx.textAlign = "center";
    for (const [v, i] of idx) {
      if (v < iMax * 0.05) continue;
      ctx.fillText(mz[i].toFixed(4), xToPx(mz[i]), yToPx(v) - 4);
    }
  }
}

// ── utilities ───────────────────────────────────────────────────────
function button(label, disabled, onClick) {
  const b = document.createElement("button");
  b.textContent = label; b.disabled = disabled;
  b.addEventListener("click", onClick);
  return b;
}
function fmt(v, digits, dash = "—") {
  return (v === null || v === undefined || Number.isNaN(v)) ? dash : Number(v).toFixed(digits);
}
function escapeHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, c =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}
function debounce(fn, ms) {
  let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
}

init();
