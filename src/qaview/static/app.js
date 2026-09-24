// qaview page: loci table, plotting controls, and a Bokeh report sized to its cell.
// All layout is CSS; the server only renders figures (docs/qaview.md §13.4).

const $ = (sel) => document.querySelector(sel);

const state = {
  project: null,
  rows: new Map(),          // locus id -> row
  current: null,            // open locus id
  sort: { key: "delta_s", asc: true },
  filter: "",
  stars: new Set(),         // starred seed names (server-side, per project)
  starOnly: false,
  viewAbort: null,
  views: [],
};

// ---- formatting --------------------------------------------------------------------

function fmtBp(bp) {
  if (bp >= 1e6) return (bp / 1e6).toFixed(bp >= 1e7 ? 0 : 1) + " Mb";
  if (bp >= 1e3) return Math.round(bp / 1e3) + " kb";
  return bp + " bp";
}

function fmtDs(x) {
  if (x === null || x === undefined) return "";
  const a = Math.abs(x), s = x > 0 ? "+" : x < 0 ? "−" : "";
  if (a >= 1e6) return s + (a / 1e6).toFixed(2) + "M";
  if (a >= 1e3) return s + (a / 1e3).toFixed(a >= 1e5 ? 0 : a >= 1e4 ? 1 : 2) + "k";
  return s + a.toFixed(1);
}

// Table cells are narrow: 301k, 1.2M.
function fmtBpShort(bp) {
  if (bp >= 1e6) return (bp / 1e6).toFixed(bp >= 1e7 ? 0 : 1) + "M";
  if (bp >= 1e3) return Math.round(bp / 1e3) + "k";
  return String(bp);
}

function dsColor(x) {
  const c = state.project.colors;
  if (x > 0) return c.ratio_pos;
  if (x < 0) return c.ratio_neg;
  return "#d2d2ca";
}

function statusText(r) {
  if (r.status === "scoring") return `${Math.round(r.progress * 100)}%`;
  return r.status === "scored" ? "done" : r.status;
}

// ---- project ----------------------------------------------------------------------

async function loadProject() {
  const p = await (await fetch("/api/project")).json();
  state.project = p;
  document.title = `qaview — ${p.name}`;
  $("#project-name").textContent = p.name;
  $("#project-pair").textContent = `${p.asm1.name} vs ${p.asm2.name}  ·  ΔS > 0 favours ${p.asm1.name}`;

  const rec = $("#recruitment");
  rec.replaceChildren();
  for (const [k, v] of Object.entries({ ...p.recruitment, ...p.resources })) {
    rec.append(el("dt", k), el("dd", String(v)));
  }
  const files = $("#files");
  files.replaceChildren();
  const f = {
    [`${p.asm1.name} fasta`]: p.asm1.fasta, [`${p.asm1.name} bam`]: p.asm1.bam,
    [`${p.asm2.name} fasta`]: p.asm2.fasta, [`${p.asm2.name} bam`]: p.asm2.bam,
    db: p.db, output: p.output_dir,
  };
  p.seed_sources.forEach((s, i) => { f[`seeds ${i + 1} (${s.format}, ${s.seeds})`] = s.path; });
  for (const [k, v] of Object.entries(f)) files.append(el("dt", k), el("dd", v));

  const layers = $("#layers");
  layers.replaceChildren();
  for (const name of p.layers) {
    const box = el("input"); box.type = "checkbox"; box.checked = true; box.value = name;
    box.addEventListener("change", scheduleView);
    const lab = el("label"); lab.append(box, name);
    layers.append(lab);
  }
}

function el(tag, text) {
  const e = document.createElement(tag);
  if (text !== undefined) e.textContent = text;
  return e;
}

// ---- loci table -------------------------------------------------------------------

async function loadLoci() {
  const d = await (await fetch("/api/loci")).json();
  for (const r of d.rows) state.rows.set(r.id, r);
  showProgress(d.progress, d.failed_seeds);
  renderTable();
}

function showProgress(p, failed) {
  const done = p.seeds ? (p.recruited + p.recruit_failed) / p.seeds : 0;
  $("#bar-recruit").style.width = `${Math.round(done * 100)}%`;
  $("#progress-text").textContent =
    `seeds recruited ${p.recruited}/${p.seeds}  ·  loci ${p.loci}, scored ${p.scored}` +
    (p.scoring ? `, scoring ${p.scoring}` : "");
  const f = $("#failed-seeds");
  const n = (failed || []).length + (p.score_failed || 0);
  f.textContent = n ? `${(failed || []).length} seeds failed recruitment, ${p.score_failed} loci failed scoring` : "";
  f.title = (failed || []).map((x) => `${x.seed}: ${x.error}`).join("\n");
}

function isStarred(r) {
  return r.seeds.some((s) => state.stars.has(s));
}

function sortValue(r, key) {
  switch (key) {
    case "star": return isStarred(r) ? 0 : 1;
    case "span": return Object.values(r.span)[0] || 0;
    case "seeds": return r.seeds[0] || "";
    case "status": return r.status;
    default: return r[key];
  }
}

function renderTable() {
  const { key, asc } = state.sort;
  const q = state.filter.toLowerCase();
  let rows = [...state.rows.values()];
  if (q) rows = rows.filter((r) => r.id.includes(q) || r.seeds.some((s) => s.toLowerCase().includes(q)));
  if (state.starOnly) rows = rows.filter(isStarred);
  rows.sort((a, b) => {
    const x = sortValue(a, key), y = sortValue(b, key);
    if (x === null || x === undefined) return 1;           // unscored last either way
    if (y === null || y === undefined) return -1;
    return (x < y ? -1 : x > y ? 1 : 0) * (asc ? 1 : -1);
  });
  for (const th of document.querySelectorAll("#loci th")) {
    th.classList.toggle("sorted", th.dataset.key === key);
    th.classList.toggle("asc", th.dataset.key === key && asc);
  }
  const names = [state.project.asm1.name, state.project.asm2.name];
  const body = $("#loci tbody");
  body.replaceChildren(...rows.map((r) => {
    const tr = el("tr");
    tr.dataset.id = r.id;
    if (r.id === state.current) tr.classList.add("selected");
    const ds = el("td"); ds.className = "num";
    if (r.delta_s !== null) {
      const dot = el("span"); dot.className = "dot"; dot.style.background = dsColor(r.delta_s);
      ds.append(dot, fmtDs(r.delta_s));
      ds.title = `ΔS ${r.delta_s.toFixed(2)} nats`;
    }
    const reads = el("td", r.reads.toLocaleString()); reads.className = "num";
    const span = el("td", `${fmtBpShort(r.span[names[0]])} / ${fmtBpShort(r.span[names[1]])}`);
    span.className = "num";
    span.title = `${names[0]}: ${r.segments[names[0]]} segments, ${r.span[names[0]].toLocaleString()} bp\n` +
                 `${names[1]}: ${r.segments[names[1]]} segments, ${r.span[names[1]].toLocaleString()} bp`;
    const seeds = el("td", r.seeds[0] + (r.seeds.length > 1 ? `  +${r.seeds.length - 1}` : ""));
    seeds.title = r.seeds.join("\n");
    const st = el("td", statusText(r)); st.className = `status-${r.status}`;
    if (r.error) st.title = r.error;
    const star = el("td"); star.className = "star";
    const sb = el("button", isStarred(r) ? "\u2605" : "\u2606");
    sb.className = isStarred(r) ? "on" : "";
    sb.title = isStarred(r) ? "unstar" : "star this locus";
    star.append(sb);
    tr.append(star, ds, reads, span, seeds, st);
    return tr;
  }));
}

function onUpdate(d) {
  showProgress(d.progress, d.failed_seeds);
  let reload = false;
  for (const r of d.rows) {
    const old = state.rows.get(r.id);
    if (r.id === state.current && old && old.status !== "scored" && r.status === "scored") reload = true;
    state.rows.set(r.id, r);
  }
  renderTable();
  if (reload) loadView();
}

// ---- controls ---------------------------------------------------------------------

function params() {
  const cell = $("#plot-cell");
  const cs = getComputedStyle(cell);
  const w = cell.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight);
  const h = cell.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom);
  const layers = [...document.querySelectorAll("#layers input:checked")].map((b) => b.value);
  return new URLSearchParams({
    w: Math.max(300, Math.floor(w)), h: Math.max(300, Math.floor(h)),
    pid: $("#pid").value, plen: $("#plen").value,
    best_pid_delta: $("#best_pid_delta").value, pid_tie_eps: $("#pid_tie_eps").value,
    binsize: $("#binsize").value || 0,
    color_by: $("#color_by .on").dataset.v, C: $("#C").value, layers: layers.join(","),
  });
}

function checkTie() {
  $("#tie-warning").hidden = !(parseFloat($("#pid_tie_eps").value) >= parseFloat($("#best_pid_delta").value));
}

let viewTimer = null;
function scheduleView() {
  checkTie();
  clearTimeout(viewTimer);
  viewTimer = setTimeout(loadView, 350);
}

// ---- the report -------------------------------------------------------------------

function clearPlot() {
  for (const v of state.views) { try { v.remove(); } catch (_) { /* already gone */ } }
  state.views = [];
  if (window.Bokeh && Bokeh.documents) {
    for (const doc of Bokeh.documents) { try { doc.clear(); } catch (_) { /* ignore */ } }
    Bokeh.documents.length = 0;
  }
  $("#plot").replaceChildren();
}

async function openLocus(id) {
  state.current = id;
  renderTable();
  fetch(`/api/loci/${id}/prioritize`, { method: "POST" });
  const r = state.rows.get(id);
  $("#locus-title").textContent = `${r.seeds[0]}${r.seeds.length > 1 ? `  (+${r.seeds.length - 1} seeds)` : ""}  ·  ${r.reads.toLocaleString()} reads`;
  $("#locus-title").title = r.seeds.join("\n");
  $("#locus-title").classList.remove("muted");
  history.replaceState(null, "", `#${id}`);
  clearPlot();
  loadView();
}

async function loadView() {
  const id = state.current;
  if (!id) return;
  if (state.viewAbort) state.viewAbort.abort();
  const ctl = new AbortController();
  state.viewAbort = ctl;
  $("#plot-status").textContent = "Loading…";
  let d;
  try {
    const res = await fetch(`/api/loci/${id}/view?${params()}`, { signal: ctl.signal });
    if (!res.ok) throw new Error(`${res.status} ${await res.text()}`);
    d = await res.json();
  } catch (e) {
    if (e.name === "AbortError") return;
    $("#plot-status").textContent = `Failed to load: ${e.message}`;
    return;
  }
  if (id !== state.current) return;
  clearPlot();
  state.views = await Bokeh.embed.embed_item(d.item, "plot");
  $("#plot-status").textContent = "";
  renderScore(d);
}

function renderScore(d) {
  const card = $("#score-card");
  card.replaceChildren();
  const r = d.locus;
  if (!d.score) {
    const s = r.status === "scoring" ? `scoring ${Math.round(r.progress * 100)}%` : r.status;
    card.append(el("span", `ΔS: ${s}`));
    card.lastChild.className = "muted";
    return;
  }
  const s = d.score, p = state.project;
  const ds = el("span", `ΔS = ${fmtDs(s.delta_s)} nats`); ds.className = "ds";
  ds.title = `${s.delta_s.toFixed(2)} nats at C = ${s.C}`;
  const parts = el("span",
    `emission ${fmtDs(s.emission)} (${s.at_plus_C} at +C, ${s.at_minus_C} at −C)  ·  ` +
    `exposure ${fmtDs(s.exposure)}  ·  ${s.scored.toLocaleString()} scored, ${s.dropped} dropped`);
  parts.className = "parts";
  const ramp = el("span"); ramp.className = "legend-ramp";
  const bar = el("i");
  bar.style.background = `linear-gradient(to right, ${p.colors.ratio_neg}, #d2d2ca, ${p.colors.ratio_pos})`;
  ramp.append(el("span", p.asm2.name), bar, el("span", p.asm1.name));
  card.append(ds, parts, ramp);
}

// ---- wiring -----------------------------------------------------------------------

function wire() {
  $("#loci tbody").addEventListener("click", (e) => {
    if (e.target.closest("td.star button")) {
      toggleStar(e.target.closest("tr").dataset.id);
      return;                                     // starring does not open the locus
    }
    const tr = e.target.closest("tr");
    if (tr) openLocus(tr.dataset.id);
  });
  for (const th of document.querySelectorAll("#loci th")) {
    th.addEventListener("click", () => {
      const k = th.dataset.key;
      state.sort = { key: k, asc: state.sort.key === k ? !state.sort.asc : true };
      renderTable();
    });
  }
  $("#filter").addEventListener("input", (e) => { state.filter = e.target.value; renderTable(); });
  for (const id of ["pid", "plen", "best_pid_delta", "pid_tie_eps", "binsize", "C"]) {
    $(`#${id}`).addEventListener("change", scheduleView);
  }
  for (const b of document.querySelectorAll("#color_by button")) {
    b.addEventListener("click", () => {
      document.querySelectorAll("#color_by button").forEach((x) => x.classList.toggle("on", x === b));
      scheduleView();
    });
  }
  $("#star-only").addEventListener("click", (e) => {
    state.starOnly = !state.starOnly;
    e.currentTarget.classList.toggle("on", state.starOnly);
    renderTable();
  });
  restoreWidths();
  wireResizers();
  for (const b of document.querySelectorAll(".collapse")) {
    b.addEventListener("click", () => {
      $("#app").classList.toggle(`${b.dataset.side}-collapsed`);
    });
  }
  // The plot follows its cell: window resizes and sidebar collapses both land here.
  let last = "";
  new ResizeObserver(() => {
    const cell = $("#plot-cell");
    const size = `${cell.clientWidth}x${cell.clientHeight}`;
    if (size !== last) { last = size; scheduleView(); }
  }).observe($("#plot-cell"));

  const es = new EventSource("/api/events");
  es.addEventListener("update", (e) => onUpdate(JSON.parse(e.data)));
  es.addEventListener("progress", (e) => showProgress(JSON.parse(e.data), null));
}

// ---- stars -----------------------------------------------------------------------

async function loadStars() {
  try {
    const d = await (await fetch("/api/stars")).json();
    state.stars = new Set(d.seeds);
  } catch (_) { /* the table still works without stars */ }
}

async function toggleStar(id) {
  const r = state.rows.get(id);
  if (!r) return;
  const res = await fetch("/api/stars", {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ seeds: r.seeds, starred: !isStarred(r) }),
  });
  if (res.ok) state.stars = new Set((await res.json()).seeds);
  renderTable();
}

// ---- resizable sidebars -----------------------------------------------------------

const WIDTHS_KEY = "qaview.sidebar-widths";

function setWidth(side, px) {
  const max = Math.max(260, window.innerWidth * 0.45);
  const w = Math.round(Math.min(max, Math.max(200, px)));
  $("#app").style.setProperty(`--${side}-w`, `${w}px`);
  return w;
}

function restoreWidths() {
  try {
    const saved = JSON.parse(localStorage.getItem(WIDTHS_KEY) || "{}");
    for (const side of ["left", "right"]) if (saved[side]) setWidth(side, saved[side]);
  } catch (_) { /* private window or blocked storage: keep the defaults */ }
}

function saveWidths() {
  try {
    const cs = getComputedStyle($("#app"));
    localStorage.setItem(WIDTHS_KEY, JSON.stringify({
      left: parseFloat(cs.getPropertyValue("--left-w")),
      right: parseFloat(cs.getPropertyValue("--right-w")),
    }));
  } catch (_) { /* not critical */ }
}

function wireResizers() {
  for (const handle of document.querySelectorAll(".resizer")) {
    handle.addEventListener("pointerdown", (e) => {
      e.preventDefault();
      const side = handle.dataset.side;
      handle.setPointerCapture(e.pointerId);
      handle.classList.add("active");
      document.body.classList.add("resizing");
      const move = (ev) => setWidth(side, side === "left" ? ev.clientX : window.innerWidth - ev.clientX);
      const up = () => {
        handle.removeEventListener("pointermove", move);
        handle.removeEventListener("pointerup", up);
        handle.removeEventListener("pointercancel", up);
        handle.classList.remove("active");
        document.body.classList.remove("resizing");
        saveWidths();
      };
      handle.addEventListener("pointermove", move);
      handle.addEventListener("pointerup", up);
      handle.addEventListener("pointercancel", up);
    });
    // Double-click restores the default width.
    handle.addEventListener("dblclick", () => {
      $("#app").style.removeProperty(`--${handle.dataset.side}-w`);
      saveWidths();
    });
  }
}

async function main() {
  await loadProject();
  wire();
  checkTie();
  await loadStars();
  await loadLoci();
  const id = location.hash.slice(1);
  if (id && state.rows.has(id)) openLocus(id);
}

main();
