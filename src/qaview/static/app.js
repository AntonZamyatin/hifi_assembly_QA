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
  alignAbort: null,
  lastRect: null,           // the rectangle the alignment panel shows
  showReads: false,
  lastAlignment: null,      // the result the panel shows, for "Remember"
  lastSettings: null,
  savedId: null,            // set when the panel shows a remembered alignment
  rects: [],                // remembered alignments on this locus, in plot coordinates
  boxes: new Map(),         // id -> its BoxAnnotation on the current dotplot
  boxPlot: null,
  views: [],
  recruitSeeds: [],         // the open locus's seeds and their recruitment parameters
  recruitSeedsOf: null,
  recruitWatch: new Map(),  // seed being re-recruited -> {from: locus id, seq}
  seq: -1,                  // the server state the rows reflect
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
  if (r.recruiting && r.recruiting.length) return "recruiting";
  if (r.status === "scoring") return `${Math.round(r.progress * 100)}%`;
  return r.status === "scored" ? "done" : r.status;
}

// ---- project ----------------------------------------------------------------------

async function loadProject(p) {
  state.project = p;
  document.title = `qaview — ${p.name}`;
  $("#project-name").textContent = p.name;
  $("#project-pair").textContent = `${p.asm1.name} vs ${p.asm2.name}  ·  ΔS > 0 favours ${p.asm1.name}`;

  buildRecruitForm(p);
  const res = $("#resources");
  res.replaceChildren();
  for (const [k, v] of Object.entries(p.resources)) res.append(el("dt", k), el("dd", String(v)));
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
  state.seq = d.seq;
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
    if (r.custom && r.custom.length) {
      const own = el("span", "*"); own.className = "own";
      seeds.append(own);
      seeds.title += `\n\nrecruited with parameters of its own: ${r.custom.join(", ")}`;
    }
    const st = el("td", statusText(r));
    st.className = r.recruiting && r.recruiting.length ? "status-scoring" : `status-${r.status}`;
    if (r.error || r.recruit_error) st.title = [r.error, r.recruit_error].filter(Boolean).join("\n");
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
  let reload = false, gone = null, seedsChanged = false;
  if (d.seq !== undefined) state.seq = Math.max(state.seq, d.seq);
  for (const r of d.rows) {
    if (r.removed) {                         // its seeds were re-recruited into other loci
      state.rows.delete(r.id);
      if (r.id === state.current) gone = r;
      continue;
    }
    const old = state.rows.get(r.id);
    if (r.id === state.current && old && old.status !== "scored" && r.status === "scored") reload = true;
    if (r.id === state.current && old &&
        String([old.seeds, old.recruiting]) !== String([r.seeds, r.recruiting])) seedsChanged = true;
    state.rows.set(r.id, r);
  }
  renderTable();
  if (checkRecruitWatch(gone)) return;       // it opened another locus
  if (seedsChanged) loadRecruitSeeds(state.current);
  if (reload) loadView();
}

// ---- controls ---------------------------------------------------------------------

const MIN_CELL = 160;                              // the server's floor for w and h

function params() {
  const cell = $("#plot-cell");
  const cs = getComputedStyle(cell);
  const w = cell.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight);
  const h = cell.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom);
  const layers = [...document.querySelectorAll("#layers input:checked")].map((b) => b.value);
  return new URLSearchParams({
    w: Math.max(MIN_CELL, Math.floor(w)), h: Math.max(MIN_CELL, Math.floor(h)),
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

async function openLocus(id, keepMsg = false) {
  if (!keepMsg) recruitMsg("");
  state.current = id;
  renderTable();
  fetch(`/api/loci/${id}/prioritize`, { method: "POST" });
  const r = state.rows.get(id);
  $("#locus-title").textContent = `${r.seeds[0]}${r.seeds.length > 1 ? `  (+${r.seeds.length - 1} seeds)` : ""}  ·  ${r.reads.toLocaleString()} reads`;
  $("#locus-title").title = r.seeds.join("\n");
  $("#locus-title").classList.remove("muted");
  history.replaceState(null, "", `#${id}`);
  clearAlignment();
  clearPlot();
  state.rects = [];
  loadRects(id);
  loadRecruitSeeds(id);
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
  drawRects();
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

// ---- bin alignment (right sidebar) ------------------------------------------------

function clearAlignment() {
  if (state.alignAbort) state.alignAbort.abort();
  $("#align-meta").textContent = "Shift+drag a rectangle on the dotplot to align the stretches under it.";
  $("#align-strip").hidden = true;
  $("#align-legend").hidden = true;
  state.lastRect = null;
  state.lastAlignment = null;
  setSaved(null);
  clearPiles();
}

function clearPiles() {
  for (const id of ["#pile-top", "#pile-bottom"]) { $(id).replaceChildren(); $(id).hidden = true; }
  $("#reads-legend").hidden = true;
}

async function alignRect(x0, x1, y0, y1) {
  const id = state.current;
  if (!id) return;
  if (state.alignAbort) state.alignAbort.abort();
  const ctl = new AbortController();
  state.alignAbort = ctl;
  $("#align-meta").textContent = "Aligning\u2026";
  state.lastRect = [x0, x1, y0, y1];
  const q = params();
  q.set("x0", x0); q.set("x1", x1); q.set("y0", y0); q.set("y1", y1);
  if (state.showReads) q.set("reads", "1");
  $("#align-meta").textContent = state.showReads ? "Aligning and fetching reads\u2026" : "Aligning\u2026";
  let a;
  try {
    const res = await fetch(`/api/loci/${id}/align?${q}`, { signal: ctl.signal });
    if (!res.ok) {
      let msg = await res.text();
      try { msg = JSON.parse(msg).detail || msg; } catch (_) { /* plain text */ }
      throw new Error(msg);
    }
    a = await res.json();
  } catch (e) {
    if (e.name !== "AbortError") {
      $("#align-meta").textContent = `Alignment failed: ${e.message}`;
      $("#align-strip").hidden = true;
      $("#align-legend").hidden = true;
    }
    return;
  }
  if (id !== state.current) return;
  const q2 = params();
  state.lastAlignment = a;
  state.lastSettings = Object.fromEntries(
    ["pid", "plen", "best_pid_delta", "pid_tie_eps", "binsize", "color_by"].map((k) => [k, q2.get(k)]));
  setSaved(null);
  renderAlignment(a);
}

function renderAlignment(a, saved = null) {
  const n = (v) => v.toLocaleString();
  const meta = $("#align-meta");
  meta.replaceChildren();
  if (saved) {
    const d = el("div", `remembered ${saved.created.replace("T", " ")}` + (saved.note ? `  \u00b7  ${saved.note}` : ""));
    d.id = "align-saved";
    meta.append(d);
  }
  const line = (label, text) => {
    const d = el("div"); const b = el("b", label); d.append(b, " ", text); return d;
  };
  meta.append(
    line(`${a.top_name}`, `${a.top_contig}:${n(a.top_start + 1)}\u2013${n(a.top_end)}  (top` +
         (a.clipped[0] ? ", kept to one segment)" : ")")),
    line(`${a.bottom_name}`, `${a.bottom_contig}:${n(a.bottom_start + 1)}\u2013${n(a.bottom_end)}  (bottom` +
         (a.strand === "-" ? ", reverse complement" : "") + (a.clipped[1] ? ", kept to one segment)" : ")")),
    el("div", `compared ${n(a.compared)} columns \u00b7 identity ${(100 * a.identity).toFixed(2)}% \u00b7 ` +
              `${n(a.mismatches)} mismatches \u00b7 ${n(a.ins)} bp in ${a.ins_opens} insertions \u00b7 ` +
              `${n(a.del)} bp in ${a.del_opens} deletions`),
    el("div", `outside the other bin: ${n(a.outside_top)} bp top, ${n(a.outside_bottom)} bp bottom`),
  );
  // Runs of equal column kind become one span per row, so a long bin stays cheap.
  // A blank in the top or bottom row marks bin sequence whose counterpart lies
  // outside the other bin: shown, but not a difference.
  const kindAt = (i) => {
    const c = a.mid[i];
    if (c === "|") return "";
    if (c === ".") return "mm";
    return a.top[i] === " " || a.bottom[i] === " " ? "out" : "gap";
  };
  const rows = [a.top, a.mid, a.bottom].map(() => []);
  let start = 0;
  for (let i = 1; i <= a.mid.length; i++) {
    if (i < a.mid.length && kindAt(i) === kindAt(start)) continue;
    const k = kindAt(start);
    [a.top, a.mid, a.bottom].forEach((row, r) => {
      const text = row.slice(start, i);
      if (!k) { rows[r].push(document.createTextNode(text)); return; }
      const sp = el("span", text); sp.className = k; rows[r].push(sp);
    });
    start = i;
  }
  const pre = $("#align-rows");
  pre.replaceChildren();
  rows.forEach((parts, r) => { pre.append(...parts); if (r < 2) pre.append("\n"); });
  const strip = $("#align-strip");
  strip.hidden = false;
  // Open on the comparison, not on sequence that lies outside the other stretch.
  const k = [...a.mid].findIndex((c, i) => c === "|" || c === "." || a.top[i] === "-" || a.bottom[i] === "-");
  const grid = measureGrid(pre);
  strip.scrollLeft = k > 0 ? Math.max(0, grid.x0 + k * grid.cw - 24) : 0;
  if (a.reads) renderPiles(a.reads, grid);
  else clearPiles();
  const top = pre.getBoundingClientRect().top - strip.getBoundingClientRect().top + strip.scrollTop;
  strip.scrollTop = Math.max(0, top - (strip.clientHeight - pre.offsetHeight) / 2);
  $("#align-legend").hidden = false;
}

// ---- reads under the alignment ------------------------------------------------------

const ROW_H = 5;                                   // 4 px bar + 1 px gap

// The alignment text's column grid, measured rather than derived: a run of
// letters placed at the start of the text gives the exact advance of one column
// and where column 0 begins, relative to the strip's left edge -- which the
// pileups share. (Padding and scrollWidth differ between browsers.)
function measureGrid(pre) {
  const probe = document.createElement("span");
  probe.textContent = "A".repeat(200);
  pre.prepend(probe);
  const r = probe.getBoundingClientRect();
  const x0 = r.left - pre.getBoundingClientRect().left;
  probe.remove();
  return { x0, cw: r.width / 200 };
}

function renderPiles(reads, grid) {
  // Both pileups grow away from the alignment: band 0 (the marginal's bottom band)
  // sits next to the sequence rows, above the top row and below the bottom row.
  const { x0, cw } = grid;
  const X = (c) => x0 + c * cw;
  for (const [side, id] of [["top", "#pile-top"], ["bottom", "#pile-bottom"]]) {
    const pile = $(id);
    const rows = reads[side];
    const n = rows.length;
    pile.style.height = `${n * ROW_H + (n ? 3 : 0)}px`;
    const parts = [];
    rows.forEach((r, i) => {
      const y = side === "top" ? (n - 1 - i) * ROW_H + 1 : i * ROW_H + 2;
      const band = reads.bands[r.band];
      // One box per read carries the tooltip; its runs are drawn inside it.
      const row = el("div");
      row.className = "read";
      row.style.cssText = `left:${X(r.col0)}px;width:${(r.col1 - r.col0 + 1) * cw}px;top:${y}px`;
      row.title = `${r.name}\n${band.name}\npid ${r.pid}  ${r.strand}${r.supplementary ? "  supplementary" : ""}\n` +
                  `${(r.start + 1).toLocaleString()}\u2013${r.end.toLocaleString()}` +
                  (r.sequence ? "" : "\n(no stored sequence: differences not shown)");
      for (const [c0, c1, kind] of r.runs) {
        const seg = el("div");
        seg.className = `run ${kind}`;
        seg.style.left = `${(c0 - r.col0) * cw}px`;
        seg.style.width = `${(c1 - c0 + 1) * cw}px`;
        if (kind === "b" || kind === "g") seg.style.background = band.color;
        row.append(seg);
      }
      parts.push(row);
      for (const [c, len] of r.ticks) {
        const t = el("div");
        t.className = "tick";
        t.title = `${len} base${len > 1 ? "s" : ""} the read has and its assembly lacks, ` +
                  "with no gap column in the alignment to put them in";
        t.style.cssText = `left:${X(c) - 1}px;top:${y - 1}px`;
        parts.push(t);
      }
    });
    pile.replaceChildren(...parts);
    pile.hidden = n === 0;
  }
  const legend = $("#reads-legend");
  const used = new Set([...reads.top, ...reads.bottom].map((r) => r.band));
  const item = (cls, color, text) => {
    const sp = el("span"); const i = el("i"); i.className = `sw ${cls}`;
    if (color) i.style.background = color;
    sp.append(i, text); return sp;
  };
  legend.replaceChildren(
    el("span", `reads: ${reads.top.length} on top, ${reads.bottom.length} below`),
    ...reads.bands.map((b, k) => used.has(k) ? item("", b.color, b.name) : null).filter(Boolean),
    item("mkx", null, "mismatch"),
    item("mki", null, "bases the read has and its assembly lacks"),
    item("mkd", null, "bases its assembly has and the read lacks"),
    item("mkg", null, "gap in its assembly too"),
    item("mkt", null, "extra bases in the read, no gap column for them"));
  legend.hidden = false;
}

function toggleReads() {
  state.showReads = !state.showReads;
  const b = $("#reads-toggle");
  b.textContent = state.showReads ? "Hide reads" : "Show reads";
  b.classList.toggle("on", state.showReads);
  if (!state.showReads) clearPiles();
  else if (state.lastRect) alignRect(...state.lastRect);
}

// ---- remembered alignments ---------------------------------------------------------

async function loadRects(id) {
  try {
    const res = await fetch(`/api/loci/${id}/remembered`);
    if (!res.ok || id !== state.current) return;
    state.rects = await res.json();
    drawRects();
  } catch (_) { /* the plot works without them */ }
}

function drawRects() {
  // Bokeh builds views for new annotations asynchronously, so boxes are changed
  // in place and added or removed in one assignment -- rebuilding them all on
  // every change raced that build ("cannot find a view for BoxAnnotation").
  const pv = dotplotView();
  if (!pv) return;
  const plot = pv.model;
  if (state.boxPlot !== plot) {                    // a new dotplot: its boxes start afresh
    state.boxPlot = plot;
    state.boxes = new Map();
  }
  const style = (id) => ({
    fill_alpha: id === state.savedId ? 0.38 : 0.2, line_width: id === state.savedId ? 2 : 1,
  });
  const want = new Map(state.rects.map((r) => [r.id, r]));
  const removed = [...state.boxes].filter(([id]) => !want.has(id)).map(([, b]) => b);
  for (const [id] of [...state.boxes]) if (!want.has(id)) state.boxes.delete(id);
  const BA = Bokeh.Models.get("BoxAnnotation");
  const added = [];
  for (const [id, r] of want) {
    const box = state.boxes.get(id);
    if (box) { box.setv(style(id)); continue; }
    const b = new BA({ left: r.x0, right: r.x1, bottom: r.y0, top: r.y1, editable: false,
                       fill_color: "#f2c200", line_color: "#c99a00", line_alpha: 0.9, ...style(id) });
    state.boxes.set(id, b);
    added.push(b);
  }
  if (removed.length || added.length) {
    plot.center = [...plot.center.filter((x) => !removed.includes(x)), ...added];
  }
}

function rectsAt(clientX, clientY) {
  const pv = dotplotView();
  if (!pv || !state.rects.length) return [];
  const cv = pv.canvas_view.el.getBoundingClientRect();
  const f = pv.frame.bbox;
  const sx = clientX - cv.left, sy = clientY - cv.top;
  if (sx < f.left || sx > f.right || sy < f.top || sy > f.bottom) return [];
  const x = pv.frame.x_scale.invert(sx), y = pv.frame.y_scale.invert(sy);
  return state.rects.filter((r) => x >= r.x0 && x <= r.x1 && y >= r.y0 && y <= r.y1)
    .sort((p, q) => (p.x1 - p.x0) * (p.y1 - p.y0) - (q.x1 - q.x0) * (q.y1 - q.y0));
}

function setSaved(id) {
  state.savedId = id;
  const b = $("#remember-toggle"), note = $("#remember-note");
  b.textContent = id ? "Forget" : "Remember";
  b.disabled = !state.lastAlignment;
  note.disabled = !state.lastAlignment || !!id;
  if (!id) note.value = "";
  drawRects();
}

async function restoreRemembered(ident) {
  let e;
  try {
    const res = await fetch(`/api/remembered/${encodeURIComponent(ident)}`);
    if (!res.ok) throw new Error(`${res.status}`);
    e = await res.json();
  } catch (err) {
    $("#align-meta").textContent = `Could not load the remembered alignment: ${err.message}`;
    return;
  }
  if (state.alignAbort) state.alignAbort.abort();
  state.lastAlignment = e.alignment;
  state.lastSettings = e.settings;
  state.lastRect = null;                           // "Show reads" does not recompute a snapshot
  renderAlignment(e.alignment, e);
  setSaved(e.id);
  $("#remember-note").value = e.note || "";
}

async function toggleRemember() {
  const b = $("#remember-toggle");
  b.disabled = true;
  try {
    if (state.savedId) {
      const res = await fetch(`/api/remembered/${encodeURIComponent(state.savedId)}`, { method: "DELETE" });
      if (!res.ok) throw new Error(`${res.status}`);
      state.rects = state.rects.filter((r) => r.id !== state.savedId);
      setSaved(null);
    } else {
      const res = await fetch("/api/remembered", {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ lid: state.current, alignment: state.lastAlignment,
                               note: $("#remember-note").value, settings: state.lastSettings || {} }),
      });
      if (!res.ok) throw new Error(`${res.status} ${await res.text()}`);
      const meta = await res.json();
      if (meta.rect) state.rects.push(meta.rect);
      renderAlignment(state.lastAlignment, meta);
      setSaved(meta.id);
      $("#remember-note").value = meta.note || "";
    }
  } catch (err) {
    $("#align-meta").prepend(el("div", `Remember failed: ${err.message}`));
  } finally {
    b.disabled = !state.lastAlignment;
  }
}

function wireRects() {
  const cell = $("#plot-cell");
  const tip = $("#rect-tip");
  let down = null;
  cell.addEventListener("pointerdown", (e) => { down = [e.clientX, e.clientY]; });
  // A click -- not the end of a pan or a drag -- inside a rectangle restores it.
  cell.addEventListener("click", (e) => {
    if (e.shiftKey || !down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4) return;
    const hit = rectsAt(e.clientX, e.clientY);
    if (hit.length) restoreRemembered(hit[0].id);
  });
  cell.addEventListener("pointermove", (e) => {
    const hit = e.buttons ? [] : rectsAt(e.clientX, e.clientY);
    if (!hit.length) { tip.hidden = true; return; }
    const r = hit[0];
    tip.textContent = `remembered ${r.created.replace("T", " ")}` + (r.note ? ` \u00b7 ${r.note}` : "") +
                      "  \u2014 click to open";
    const c = cell.getBoundingClientRect();
    tip.hidden = false;
    // Beside the pointer, flipped to the other side near the cell's right or bottom edge.
    let x = e.clientX - c.left + 14, y = e.clientY - c.top + 14;
    if (x + tip.offsetWidth > c.width) x = Math.max(0, e.clientX - c.left - 14 - tip.offsetWidth);
    if (y + tip.offsetHeight > c.height) y = Math.max(0, e.clientY - c.top - 14 - tip.offsetHeight);
    tip.style.left = `${x}px`;
    tip.style.top = `${y}px`;
  });
  cell.addEventListener("pointerleave", () => { tip.hidden = true; });
}

// ---- shift+drag a rectangle on the dotplot -----------------------------------------

function dotplotView() {
  if (!window.Bokeh) return null;
  const views = [];
  const walk = (v) => { if (!v) return; views.push(v); for (const c of (v.child_views || [])) walk(c); };
  for (const v of Object.values(Bokeh.index)) walk(v);
  // The dotplot is the figure owning the wheel zoom. The grid shows the same
  // toolbar at its edge, so match figures only (they have a canvas view).
  return views.find((v) => v.canvas_view && v.model && v.model.toolbar &&
    v.model.toolbar.tools.some((t) => t.constructor.__name__ === "WheelZoomTool")) || null;
}

function wireRectAlign() {
  const cell = $("#plot-cell");
  const box = el("div");
  box.id = "align-rect";
  box.hidden = true;
  cell.append(box);

  // Capture phase on an ancestor: a shift+press on the dotplot is ours and never
  // reaches Bokeh, so its pan tool does not move the plot. Without shift Bokeh
  // gets everything as before.
  const grab = (e) => {
    if (!e.shiftKey || e.button !== 0 || !state.current) return;
    const pv = dotplotView();
    if (!pv) return;
    const cv = pv.canvas_view.el.getBoundingClientRect();
    const f = pv.frame.bbox;
    const fx0 = cv.left + f.left, fx1 = cv.left + f.right, fy0 = cv.top + f.top, fy1 = cv.top + f.bottom;
    if (e.clientX < fx0 || e.clientX > fx1 || e.clientY < fy0 || e.clientY > fy1) return;
    e.preventDefault();
    e.stopPropagation();
    if (e.type !== "pointerdown") return;

    const clampX = (x) => Math.min(fx1, Math.max(fx0, x));
    const clampY = (y) => Math.min(fy1, Math.max(fy0, y));
    const x0 = e.clientX, y0 = e.clientY;
    let x1 = x0, y1 = y0;
    const draw = () => {
      const c = cell.getBoundingClientRect();
      box.style.left = `${Math.min(x0, x1) - c.left}px`;
      box.style.top = `${Math.min(y0, y1) - c.top}px`;
      box.style.width = `${Math.abs(x1 - x0)}px`;
      box.style.height = `${Math.abs(y1 - y0)}px`;
    };
    draw();
    box.hidden = false;
    const move = (ev) => {
      ev.preventDefault(); ev.stopPropagation();
      x1 = clampX(ev.clientX); y1 = clampY(ev.clientY);
      draw();
    };
    const up = (ev) => {
      ev.preventDefault(); ev.stopPropagation();
      window.removeEventListener("pointermove", move, true);
      window.removeEventListener("pointerup", up, true);
      box.hidden = true;
      const dx = (sx) => pv.frame.x_scale.invert(sx - cv.left);
      const dy = (sy) => pv.frame.y_scale.invert(sy - cv.top);
      alignRect(dx(x0), dx(x1), dy(y0), dy(y1));     // a shift+click aligns the bins under it
    };
    window.addEventListener("pointermove", move, true);
    window.addEventListener("pointerup", up, true);
  };
  cell.addEventListener("pointerdown", grab, true);
  cell.addEventListener("mousedown", grab, true);
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
  wireRectAlign();
  $("#reads-toggle").addEventListener("click", toggleReads);
  $("#remember-toggle").addEventListener("click", toggleRemember);
  $("#recruit-run").addEventListener("click", runRecruit);
  $("#recruit-reset").addEventListener("click", () => fillRecruitForm(state.project.recruitment));
  wireRects();
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

// ---- recruitment -----------------------------------------------------------------

// Arrow steps; any value in range can be typed.
const REC_STEPS = { pid_threshold: 0.0005, plen_threshold: 0.0005, best_pid_delta: 0.0001,
                    merge_gap: 1000, extension_window: 10, ext_max: 10000 };

function buildRecruitForm(p) {
  const form = $("#recruit-form");
  form.replaceChildren();
  for (const [k, [kind, lo, hi]] of Object.entries(p.recruitment_ranges)) {
    const input = el("input");
    input.type = "number"; input.dataset.key = k;
    input.min = lo; input.max = hi;
    input.step = REC_STEPS[k] || (kind === "int" ? 1 : "any");
    input.value = p.recruitment[k];
    input.addEventListener("input", markRecruitChanges);
    const lab = el("label", k); lab.append(input);
    form.append(lab);
  }
  markRecruitChanges();
}

function recruitInputs() {
  return [...document.querySelectorAll("#recruit-form [data-key]")];
  for (const [k, names] of Object.entries(p.recruitment_choices || {})) {
    const sel = el("select"); sel.dataset.key = k;
    for (const n of names) sel.append(el("option", n));
    sel.value = p.recruitment[k];
    sel.addEventListener("change", markRecruitChanges);
    const lab = el("label", k); lab.append(sel);
    form.append(lab);
  }
}

function fillRecruitForm(params) {
  for (const i of recruitInputs()) i.value = params[i.dataset.key];
  markRecruitChanges();
}

// A value unlike the project file's is marked, with the file's value in its tooltip.
function markRecruitChanges() {
  const def = state.project.recruitment;
  for (const i of recruitInputs()) {
    const d = def[i.dataset.key];
    const changed = i.tagName === "SELECT" ? i.value !== d : i.value !== "" && Number(i.value) !== d;
    i.classList.toggle("changed", changed);
    i.title = changed ? `project file: ${d}` : "as in the project file";
  }
}

function recruitMsg(text, warn = false) {
  const m = $("#recruit-msg");
  m.replaceChildren();
  if (text instanceof Node) m.append(text); else m.textContent = text || "";
  m.classList.toggle("warn", warn);
  m.hidden = !text;
}

function updateRecruitButtons() {
  const open = !!state.current && state.recruitSeeds.length > 0;
  const ticked = document.querySelectorAll("#recruit-seeds input:checked").length;
  $("#recruit-run").disabled = !open || ticked === 0;
  $("#recruit-reset").disabled = !open;
}

async function loadRecruitSeeds(id) {
  // Reloading the same locus keeps what was ticked; another locus starts all ticked.
  const same = state.recruitSeeds.length && state.recruitSeedsOf === id;
  const unticked = new Set(same ? [...document.querySelectorAll("#recruit-seeds input:not(:checked)")]
    .map((b) => b.value) : []);
  state.recruitSeeds = [];
  updateRecruitButtons();
  let d;
  try {
    const res = await fetch(`/api/loci/${id}/seeds`);
    if (!res.ok) return;
    d = await res.json();
  } catch (_) { return; }
  if (id !== state.current) return;
  state.recruitSeeds = d.seeds;
  state.recruitSeedsOf = id;
  const several = d.seeds.length > 1;
  $("#recruit-seeds").replaceChildren(...d.seeds.map((s) => {
    const lab = el("label");
    const cb = el("input"); cb.type = "checkbox"; cb.checked = !unticked.has(s.name); cb.value = s.name;
    cb.hidden = !several;                     // one seed: nothing to choose
    cb.addEventListener("change", updateRecruitButtons);
    const name = el("span", s.name); name.className = "seed-name";
    const where = el("span", `${s.contig}:${(s.start + 1).toLocaleString()}-${s.end.toLocaleString()}`);
    where.className = "muted";
    lab.append(cb, name, where);
    lab.title = `${s.assembly} ${s.contig}:${(s.start + 1).toLocaleString()}-${s.end.toLocaleString()}` +
                ` (${(s.end - s.start).toLocaleString()} bp), from ${s.source}`;
    const own = Object.entries(s.own);
    if (s.pending) {
      const b = el("span", "recruiting\u2026"); b.className = "badge busy"; lab.append(b);
    } else if (own.length) {
      const b = el("span", "own parameters"); b.className = "badge own";
      b.title = own.map(([k, v]) => `${k} = ${v}`).join("\n");
      lab.append(b);
    }
    return lab;
  }));
  const first = d.seeds[0];
  fillRecruitForm(first.pending || first.params);
  $("#recruit-hint").textContent = several
    ? "Tick the seeds to recruit again; the values shown are the first seed's." : "";
  $("#recruit-hint").hidden = !several;
  updateRecruitButtons();
}

async function runRecruit() {
  const seeds = [...document.querySelectorAll("#recruit-seeds input:checked")].map((b) => b.value);
  const params = {};
  for (const i of recruitInputs()) params[i.dataset.key] = i.value;
  const from = state.current;
  $("#recruit-run").disabled = true;
  try {
    const res = await fetch("/api/recruit", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ seeds, params }),
    });
    if (!res.ok) {
      let msg = await res.text();
      try { msg = JSON.parse(msg).detail || msg; } catch (_) { /* plain text */ }
      recruitMsg(msg, true);
      return;
    }
    const d = await res.json();
    for (const s of d.queued) state.recruitWatch.set(s, { from, seq: d.seq });
    const parts = [];
    if (d.queued.length) parts.push(`Recruiting ${d.queued.join(", ")}\u2026`);
    if (d.unchanged.length) parts.push(`Already recruited with these parameters: ${d.unchanged.join(", ")}.`);
    recruitMsg(parts.join(" "));
    // A recruitment already on disk can be through before this reply: settle now.
    // The seed list follows the row updates (recruiting, then moved or back).
    checkRecruitWatch(null);
  } finally {
    updateRecruitButtons();
  }
}

// Seeds being re-recruited settle when no row lists them as recruiting: each is
// then in the row of its new locus, or still in its old one (same result, or a
// failure). Only rows at least as new as the request count: before that a row
// shows the seed as it was. Returns true when it opened another locus.
function checkRecruitWatch(gone) {
  if (!state.recruitWatch.size) {
    if (gone && gone.moved_to && state.rows.has(gone.moved_to)) { openLocus(gone.moved_to, true); return true; }
    return false;
  }
  const busy = new Set(), home = new Map();
  for (const r of state.rows.values()) {
    for (const s of r.recruiting || []) busy.add(s);
    for (const s of r.seeds) home.set(s, r);
  }
  const notes = [];
  let open = null, warn = false;
  for (const [seed, { from, seq }] of [...state.recruitWatch]) {
    const r = home.get(seed);
    if (state.seq < seq || busy.has(seed) || !r) continue;
    state.recruitWatch.delete(seed);
    if (r.id === from) {
      const err = (r.recruit_error || "").startsWith(`${seed}:`) ? r.recruit_error : null;
      if (err) { notes.push(`Re-recruiting failed: ${err}`); warn = true; }
      else notes.push(`${seed}: the same locus as before.`);
    } else {
      notes.push(`${seed}: moved to another locus, ${r.reads.toLocaleString()} reads.`);
      if (from === state.current && gone) open = open || r.id;
      else if (from === state.current) {
        const a = el("a", `open ${seed}'s locus`); a.href = "#";
        a.addEventListener("click", (e) => { e.preventDefault(); openLocus(r.id); });
        notes.push(a);
      }
    }
  }
  if (!open && gone && gone.moved_to && state.rows.has(gone.moved_to)) open = gone.moved_to;
  if (notes.length) {
    const frag = document.createDocumentFragment();
    notes.forEach((n, i) => { if (i) frag.append(" "); frag.append(n); });
    recruitMsg(frag, warn);
  }
  if (open) { openLocus(open, true); return true; }
  return false;
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
  // The CSS caps a sidebar at 40% / 30% of the window as well; keep the two in step.
  const max = Math.max(200, window.innerWidth * (side === "left" ? 0.4 : 0.3));
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

// ---- opening a project ------------------------------------------------------------

function showOpen(show, p, cancellable) {
  $("#open-project").hidden = !show;
  $("#open-cancel").hidden = !cancellable;
  $("#open-error").hidden = true;
  if (!show) return;
  const input = $("#open-path");
  if (p && p.open && !input.value) input.value = p.path;
  const list = $("#recent");
  list.replaceChildren();
  const recent = (p && p.recent) || [];
  for (const r of recent) {
    const li = el("li"), b = el("button");
    b.type = "button";
    const cut = r.path.lastIndexOf("/");
    const name = el("b", r.path.slice(cut + 1));
    b.append(name, "  ", el("span", r.path.slice(0, cut)));
    b.title = r.exists ? r.path : `${r.path}\n(no longer there)`;
    if (r.exists) b.addEventListener("click", () => openProject(r.path));
    else b.className = "missing";
    li.append(b);
    list.append(li);
  }
  $("#recent-wrap").hidden = recent.length === 0;
  input.focus();
}

async function openProject(path) {
  const err = $("#open-error");
  err.hidden = true;
  $("#open-path").value = path;
  const submit = $("#open-form button[type=submit]");
  submit.disabled = true;
  submit.textContent = "Opening\u2026";
  try {
    const res = await fetch("/api/open", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path }),
    });
    if (!res.ok) {
      let msg = await res.text();
      try { msg = JSON.parse(msg).detail || msg; } catch (_) { /* plain text */ }
      throw new Error(msg);
    }
    // A fresh page for the new project: every table, view and stream starts clean.
    history.replaceState(null, "", location.pathname);
    location.reload();
  } catch (e) {
    err.textContent = `Cannot open this project:\n${e.message}`;
    err.hidden = false;
    submit.disabled = false;
    submit.textContent = "Open";
  }
}

function wireOpen(p) {
  $("#open-form").addEventListener("submit", (e) => {
    e.preventDefault();
    const path = $("#open-path").value.trim();
    if (path) openProject(path);
  });
  $("#open-cancel").addEventListener("click", () => showOpen(false, p, false));
  $("#change-project").addEventListener("click", (e) => {
    e.preventDefault();
    showOpen($("#open-project").hidden, p, true);
  });
}

async function main() {
  const p = await (await fetch("/api/project")).json();
  wireOpen(p);
  if (!p.open) {
    // No project yet: the form is the page. Nothing else is wired or fetched.
    $("#app").classList.add("no-project");
    document.title = "qaview";
    $("#locus-title").textContent = "Open a project in the left sidebar to begin.";
    showOpen(true, p, false);
    return;
  }
  await loadProject(p);
  wire();
  checkTie();
  await loadStars();
  await loadLoci();
  const id = location.hash.slice(1);
  if (id && state.rows.has(id)) openLocus(id);
}

main();
