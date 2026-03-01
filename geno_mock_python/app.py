#!/usr/bin/env python3
import csv
import json
import os
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


def clean_text(value):
    if value is None:
        return ""
    text = str(value)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("_x000D_", "")
    return text


def read_csv_rows(path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cleaned = {k: clean_text(v) for k, v in row.items()}
            rows.append(cleaned)
    return rows


ROOT_DIR = Path(__file__).resolve().parent.parent
GENE_CSV = ROOT_DIR / "Gene_table.csv"
ALLELE_CSV = ROOT_DIR / "Allele_table.csv"
VARIANT_CSV = ROOT_DIR / "Variant_table.csv"
EXON_CSV = ROOT_DIR / "Exon_table.csv"
BRIDGE_CSV = ROOT_DIR / "Bridge_table.csv"

gene_table = read_csv_rows(GENE_CSV)
allele_table = read_csv_rows(ALLELE_CSV)
variant_table = read_csv_rows(VARIANT_CSV)
exon_table = read_csv_rows(EXON_CSV)
bridge_table = read_csv_rows(BRIDGE_CSV)

for row in allele_table:
    row["Allele_id"] = clean_text(row.get("Allele_id", ""))
for row in variant_table:
    row["Variant_id"] = clean_text(row.get("Variant_id", ""))
for row in bridge_table:
    row["Allele_id"] = clean_text(row.get("Allele_id", ""))
    row["Variant_id"] = clean_text(row.get("Variant_id", ""))

bridge_by_allele = {}
bridge_by_variant = {}
for row in bridge_table:
    aid = row.get("Allele_id", "")
    vid = row.get("Variant_id", "")
    if aid and vid:
        bridge_by_allele.setdefault(aid, []).append(vid)
        bridge_by_variant.setdefault(vid, []).append(aid)


def unique_keep_order(items):
    out = []
    seen = set()
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


for row in allele_table:
    aid = row.get("Allele_id", "")
    linked = unique_keep_order(bridge_by_allele.get(aid, []))
    row["Variants"] = "\n".join(linked)

all_genes_set = set()
for source in (gene_table, allele_table, variant_table, exon_table):
    for row in source:
        gene = clean_text(row.get("Gene", "")).strip()
        if gene:
            all_genes_set.add(gene)
all_genes = sorted(all_genes_set)

gene_view_cols = list(gene_table[0].keys()) if gene_table else []
if "Gene_id" in gene_view_cols and "Gene" in gene_view_cols:
    gene_view_cols = [col for col in gene_view_cols if col != "Gene"]
    gene_id_idx = gene_view_cols.index("Gene_id")
    gene_view_cols.insert(gene_id_idx + 1, "Gene")
allele_view_cols = list(allele_table[0].keys()) if allele_table else []
if "Variants" not in allele_view_cols:
    allele_view_cols.append("Variants")
variant_view_cols = list(variant_table[0].keys()) if variant_table else []
exon_view_cols = list(exon_table[0].keys()) if exon_table else []

INIT_PAYLOAD = {
    "gene_table": gene_table,
    "allele_table": allele_table,
    "variant_table": variant_table,
    "exon_table": exon_table,
    "bridge_by_allele": bridge_by_allele,
    "bridge_by_variant": bridge_by_variant,
    "all_genes": all_genes,
    "counts": {
        "gene": len(gene_table),
        "allele": len(allele_table),
        "variant": len(variant_table),
        "exon": len(exon_table),
        "bridge": len(bridge_table),
    },
    "columns": {
        "gene_table": gene_view_cols,
        "allele_view": allele_view_cols,
        "variant_view": variant_view_cols,
        "exon_view": exon_view_cols,
    },
}

FEEDBACK_LOG = []

INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Geno Database (Python)</title>
  <style>
    html, body {
      margin: 0;
      padding: 0;
      height: 100%;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: #f3f4f6;
      color: #1f2937;
    }
    * { box-sizing: border-box; }
    .app-shell {
      min-height: 100vh;
      padding: 12px;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .app-root {
      display: flex;
      gap: 12px;
      min-height: calc(100vh - 36px);
    }
    .panel-a {
      flex: 0 0 80%;
      min-width: 0;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .panel-b {
      flex: 0 0 20%;
      min-width: 0;
      min-height: 0;
      height: 100%;
      display: flex;
      flex-direction: column;
      gap: 12px;
      overflow: hidden;
    }
    .app-root.panel-b-collapsed .panel-a {
      flex: 1 1 100%;
    }
    .app-root.panel-b-collapsed .panel-b {
      display: none;
    }
    .section-box {
      border: 1px solid #c8c8c8;
      border-radius: 8px;
      background: #fff;
      padding: 10px;
      min-width: 0;
    }
    .a1-box {
      min-height: 78px;
      display: flex;
      align-items: center;
    }
    .a2-box {
      flex: 1 1 auto;
      min-height: 0;
      min-width: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .mode-row {
      display: flex;
      gap: 8px;
      width: 100%;
      align-items: center;
      flex-wrap: wrap;
    }
    .sidebar-toggle-wrap {
      margin-left: auto;
    }
    .toggle-btn {
      border: 1px solid #b8b8b8;
      border-radius: 6px;
      background: #f7f7f7;
      padding: 10px 14px;
      font-size: 15px;
      cursor: pointer;
      min-width: 120px;
    }
    .toggle-btn.active {
      background: #dbeaf8;
      border-color: #6a9ecf;
      font-weight: 600;
    }
    .download-actions {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
    }
    .download-note {
      margin-top: 6px;
      font-size: 0.92em;
      color: #555;
    }
    .b1-box {
      flex: 0 0 26%;
      min-height: 120px;
      max-height: 220px;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .b1-box h4 {
      margin: 0 0 6px 0;
      font-size: 14px;
    }
    #latest-updates {
      font-size: 12px;
      line-height: 1.3;
      overflow-y: auto;
      min-height: 0;
      padding-right: 4px;
    }
    #latest-updates ul {
      margin: 0;
      padding-left: 16px;
    }
    #latest-updates li {
      margin-bottom: 4px;
    }
    .b2-box {
      flex: 1 1 auto;
      min-height: 0;
      overflow-y: auto;
    }
    .a2-allele-layout {
      display: flex;
      gap: 12px;
      flex: 1 1 auto;
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .a2-allele-a {
      flex: 0 0 24%;
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .gene-search-wrap {
      flex: 0 0 auto;
      min-height: 0;
    }
    .a2-allele-b {
      flex: 1 1 76%;
      min-width: 0;
      height: 100%;
      min-height: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .gene-search-wrap label,
    .detail-controls-grid label {
      display: block;
      margin-bottom: 4px;
      font-weight: 600;
    }
    .input {
      width: 100%;
      border: 1px solid #c5c9cf;
      border-radius: 6px;
      padding: 8px 10px;
      font-size: 14px;
      background: #fff;
    }
    .gene-list-wrap {
      margin-top: 6px;
      border: 1px solid #e2e2e2;
      border-radius: 6px;
      background: #fbfbfb;
      overflow: auto;
      display: flex;
      flex-direction: column;
      flex: 1 1 auto;
      min-height: 0;
    }
    .detail-controls {
      flex: 0 0 auto;
      margin-bottom: 8px;
    }
    .detail-controls-grid {
      display: grid;
      grid-template-columns: minmax(220px, 1fr) 140px 140px 170px;
      gap: 8px;
      align-items: end;
    }
    .a2-allele-layout.gene-selector-collapsed .a2-allele-a {
      display: none;
    }
    .a2-allele-layout.gene-selector-collapsed .a2-allele-b {
      flex: 1 1 100%;
    }
    .detail-table-wrap {
      flex: 1 1 auto;
      min-height: 0;
      overflow: hidden;
    }
    .table-shell {
      border: 1px solid #d9dce1;
      border-radius: 6px;
      background: #fff;
      overflow-x: auto;
      overflow-y: auto;
      width: 100%;
      min-width: 0;
    }
    table.data-table {
      width: max-content;
      min-width: 100%;
      border-collapse: collapse;
      font-size: 14px;
      background: #fff;
    }
    table.data-table th,
    table.data-table td {
      border-bottom: 1px solid #e5e7eb;
      padding: 8px 10px;
      text-align: left;
      vertical-align: top;
      white-space: nowrap;
    }
    table.data-table thead th {
      position: sticky;
      top: 0;
      z-index: 2;
      background: #f8fafc;
      font-weight: 700;
      border-bottom: 1px solid #d1d5db;
    }
    table.data-table td.multiline {
      white-space: pre-line;
    }
    table.data-table th.allele-name-col,
    table.data-table td.allele-name-col {
      width: 25ch;
      min-width: 25ch;
      max-width: 25ch;
      white-space: normal;
      overflow-wrap: anywhere;
      word-break: break-word;
    }
    table.data-table tr.selected-row {
      background: #dbeafe;
    }
    table.data-table tr:hover {
      background: #f3f6fb;
    }
    .check-cell {
      width: 34px;
      text-align: center;
    }
    .hidden {
      display: none !important;
    }
    .last-feedback {
      margin-top: 10px;
      border-top: 1px solid #e6e6e6;
      padding-top: 8px;
      font-size: 0.92em;
    }
    .muted {
      color: #6b7280;
    }
    .modal-overlay {
      position: fixed;
      inset: 0;
      background: rgba(31, 41, 55, 0.45);
      display: none;
      align-items: center;
      justify-content: center;
      z-index: 1000;
      padding: 18px;
    }
    .modal-overlay.open {
      display: flex;
    }
    .modal-box {
      width: min(1200px, 95vw);
      max-height: 90vh;
      background: #fff;
      border-radius: 10px;
      border: 1px solid #d1d5db;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .modal-head {
      padding: 12px 14px;
      border-bottom: 1px solid #e5e7eb;
      font-weight: 700;
    }
    .modal-body {
      padding: 12px;
      min-height: 0;
      overflow: hidden;
      flex: 1 1 auto;
    }
    .modal-footer {
      padding: 10px 12px;
      border-top: 1px solid #e5e7eb;
      display: flex;
      gap: 8px;
      justify-content: flex-end;
      flex-wrap: wrap;
    }
    .btn {
      border: 1px solid #b8b8b8;
      border-radius: 6px;
      background: #f7f7f7;
      padding: 8px 12px;
      font-size: 14px;
      cursor: pointer;
    }
    .btn-primary {
      background: #2563eb;
      color: #fff;
      border-color: #1d4ed8;
    }
    .status-ok {
      color: #047857;
      font-size: 13px;
      margin-top: 8px;
    }
    @media (max-width: 1100px) {
      .app-root {
        flex-direction: column;
      }
      .panel-a, .panel-b {
        flex: 1 1 auto;
      }
      .a2-allele-layout {
        flex-direction: column;
      }
      .a2-allele-a, .a2-allele-b {
        max-height: none;
      }
      .detail-controls-grid {
        grid-template-columns: 1fr;
      }
      .gene-list-wrap {
        max-height: 40vh;
      }
    }
  </style>
</head>
<body>
  <div class="app-shell">
    <div class="app-root">
      <div class="panel-a">
        <div class="section-box a1-box">
          <div class="mode-row">
            <button id="mode-gene" class="toggle-btn active">Gene Table</button>
            <button id="mode-allele" class="toggle-btn">Allele Table</button>
            <button id="mode-exon" class="toggle-btn">Exon Table</button>
            <div class="sidebar-toggle-wrap">
              <button id="toggle-panel-b" class="btn">Hide Sidebar</button>
            </div>
          </div>
        </div>
        <div class="section-box a2-box">
          <div id="gene-mode"></div>
          <div id="exon-mode" class="hidden"></div>
          <div id="allele-mode" class="hidden">
            <div id="allele-layout" class="a2-allele-layout">
              <div class="a2-allele-a">
                <div class="gene-search-wrap">
                  <label for="gene-search">Search genes</label>
                  <input id="gene-search" class="input" type="text" placeholder="Type to filter genes" />
                </div>
                <div class="gene-list-wrap">
                  <div id="gene-selector-table"></div>
                </div>
              </div>
              <div class="a2-allele-b">
                <div class="detail-controls">
                  <div class="detail-controls-grid">
                    <div>
                      <label for="table-search">Search table</label>
                      <input id="table-search" class="input" type="text" placeholder="Filter rows in active table" />
                    </div>
                    <button id="gran-allele" class="toggle-btn active">Allele</button>
                    <button id="gran-variant" class="toggle-btn">Variant</button>
                    <button id="toggle-gene-selector" class="btn">Hide Gene Selector</button>
                  </div>
                </div>
                <div class="detail-table-wrap">
                  <div id="detail-table"></div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div class="panel-b">
        <div class="section-box">
          <div class="download-actions">
            <button id="download-main-btn" class="btn">Download CSV</button>
            <button id="clear-main-btn" class="btn">Clear All</button>
          </div>
          <div id="download-context" class="download-note"></div>
        </div>
        <div class="section-box b1-box">
          <h4>Latest updates</h4>
          <div id="latest-updates"></div>
        </div>
        <div class="section-box b2-box">
          <h4>Feedback form</h4>
          <label for="feedback-state">State</label>
          <input id="feedback-state" class="input" type="text" readonly />
          <div style="height: 8px;"></div>
          <label for="feedback-email">Email</label>
          <input id="feedback-email" class="input" type="text" />
          <div style="height: 8px;"></div>
          <label for="feedback-comment">Comment</label>
          <textarea id="feedback-comment" class="input" rows="5"></textarea>
          <div style="height: 10px;"></div>
          <button id="feedback-submit" class="btn btn-primary">Submit</button>
          <div id="feedback-flash" class="status-ok"></div>
          <div id="feedback-status" class="last-feedback muted">No submissions yet.</div>
        </div>
      </div>
    </div>
  </div>

  <div id="relation-modal" class="modal-overlay">
    <div class="modal-box">
      <div id="modal-title" class="modal-head"></div>
      <div class="modal-body">
        <div id="modal-table"></div>
      </div>
      <div class="modal-footer">
        <button id="modal-download-btn" class="btn">Download CSV</button>
        <button id="modal-clear-btn" class="btn">Clear All</button>
        <button id="modal-close-btn" class="btn">Close</button>
      </div>
    </div>
  </div>

  <script>
    let DATA = null;

    const state = {
      mode: "GT",
      selectedGene: "",
      geneSearch: "",
      granularity: "Allele",
      tableSearch: "",
      panelBCollapsed: false,
      alleleSelectorCollapsed: false,
      activeTable: "gene",
      checked: {
        gene: new Set(),
        exon: new Set(),
        detail: new Set(),
        modalVariants: new Set(),
        modalAlleles: new Set()
      },
      modal: {
        open: false,
        type: null,
        rows: [],
        title: ""
      }
    };

    function text(v) {
      if (v === null || v === undefined) return "";
      return String(v);
    }

    function ciContains(haystack, needle) {
      return text(haystack).toLowerCase().includes(text(needle).toLowerCase());
    }

    function parseSearchGroups(rawQuery) {
      const query = text(rawQuery).trim();
      if (!query) return [];
      return query
        .split(";")
        .map((andPart) => andPart.trim())
        .filter(Boolean)
        .map((andPart) => andPart.split(",")
          .map((orTerm) => orTerm.trim())
          .filter(Boolean));
    }

    function matchesSearchGroups(haystack, rawQuery) {
      const groups = parseSearchGroups(rawQuery);
      if (groups.length === 0) return true;
      const source = text(haystack).toLowerCase();
      return groups.every((orGroup) => orGroup.some((term) => source.includes(term.toLowerCase())));
    }

    function uniqueOrdered(list) {
      const out = [];
      const seen = new Set();
      list.forEach((item) => {
        const key = text(item);
        if (!seen.has(key)) {
          seen.add(key);
          out.push(key);
        }
      });
      return out;
    }

    function normalizeRow(row, columns) {
      const out = {};
      columns.forEach((col) => {
        out[col] = row[col] === undefined || row[col] === null ? "" : row[col];
      });
      return out;
    }

    function getFilteredGenes() {
      if (!state.geneSearch.trim()) return DATA.all_genes.slice();
      return DATA.all_genes.filter((g) => matchesSearchGroups(g, state.geneSearch));
    }

    function ensureSelectedGene() {
      const genes = getFilteredGenes();
      if (genes.length === 0) {
        state.selectedGene = "";
        return;
      }
      if (!genes.includes(state.selectedGene)) {
        state.selectedGene = genes[0];
        state.checked.detail.clear();
      }
    }

    function detailBaseRows() {
      if (!state.selectedGene) return [];
      const source = state.granularity === "Allele" ? DATA.allele_table : DATA.variant_table;
      const rows = source.filter((r) => text(r.Gene) === state.selectedGene);
      if (!state.tableSearch.trim()) return rows;
      return rows.filter((row) => {
        const joined = Object.values(row).map(text).join(" || ");
        return matchesSearchGroups(joined, state.tableSearch);
      });
    }

    function mapRowsByIds(sourceRows, idField, ids) {
      const byId = new Map();
      sourceRows.forEach((row) => byId.set(text(row[idField]), row));
      const out = [];
      ids.forEach((id) => {
        const row = byId.get(text(id));
        if (row) out.push(row);
      });
      return out;
    }

    function setActiveButtons() {
      document.getElementById("mode-gene").classList.toggle("active", state.mode === "GT");
      document.getElementById("mode-allele").classList.toggle("active", state.mode === "AT");
      document.getElementById("mode-exon").classList.toggle("active", state.mode === "ET");
      document.getElementById("gran-allele").classList.toggle("active", state.granularity === "Allele");
      document.getElementById("gran-variant").classList.toggle("active", state.granularity === "Variant");
    }

    function currentVisibleStateText() {
      const modeLabel = state.mode === "GT" ? "Gene" : (state.mode === "AT" ? "Allele" : "Exon");
      const geneLabel = state.mode === "AT" ? (state.selectedGene || "") : "N/A";
      const viewLabel = state.mode === "AT"
        ? state.granularity + "_Table"
        : (state.mode === "ET" ? "Exon_Table" : "N/A");
      return "Mode: " + modeLabel
        + " | Gene: " + geneLabel
        + " | B2 View: " + viewLabel
        + " | GeneSearch: " + state.geneSearch
        + " | TableSearch: " + state.tableSearch;
    }

    function updateStateField() {
      document.getElementById("feedback-state").value = currentVisibleStateText();
    }

    function activeTableInfo() {
      if (state.activeTable === "gene") {
        return { label: "Gene Table", checked: state.checked.gene };
      }
      if (state.activeTable === "exon") {
        return { label: "Exon Table", checked: state.checked.exon };
      }
      if (state.activeTable === "detail") {
        const label = state.granularity === "Allele" ? "Allele Detail Table" : "Variant Detail Table";
        return { label: label, checked: state.checked.detail };
      }
      if (state.activeTable === "modalVariants") {
        return { label: "Variants Popup Table", checked: state.checked.modalVariants };
      }
      if (state.activeTable === "modalAlleles") {
        return { label: "Alleles Popup Table", checked: state.checked.modalAlleles };
      }
      return { label: "Unknown", checked: new Set() };
    }

    function updateDownloadContext() {
      const info = activeTableInfo();
      document.getElementById("download-context").textContent =
        "Active table: " + info.label + " | Checked rows: " + info.checked.size;
    }

    function applyLayoutState() {
      const appRoot = document.querySelector(".app-root");
      const panelBtn = document.getElementById("toggle-panel-b");
      const alleleLayout = document.getElementById("allele-layout");
      const alleleBtn = document.getElementById("toggle-gene-selector");

      if (appRoot) {
        appRoot.classList.toggle("panel-b-collapsed", state.panelBCollapsed);
      }
      if (panelBtn) {
        panelBtn.textContent = state.panelBCollapsed ? "Show Sidebar" : "Hide Sidebar";
      }
      if (alleleLayout) {
        alleleLayout.classList.toggle("gene-selector-collapsed", state.alleleSelectorCollapsed);
      }
      if (alleleBtn) {
        alleleBtn.textContent = state.alleleSelectorCollapsed ? "Show Gene Selector" : "Hide Gene Selector";
      }
    }

    function renderTable(containerId, config) {
      const container = document.getElementById(containerId);
      container.innerHTML = "";

      const rows = config.rows || [];
      const columns = config.columns || [];
      const checkedSet = config.checkedSet;
      const maxHeight = config.maxHeight || "52vh";
      const multilineCols = new Set(config.multilineCols || []);
      const truncateTooltipCols = new Set(config.truncateTooltipCols || []);
      const truncateLength = Number.isInteger(config.truncateLength) ? config.truncateLength : 10;
      const wrapCols = new Set(["Allele_name", "Reference_No_PMID", "Accession_number", "Phenotype"]);
      const rowKey = typeof config.rowKey === "function"
        ? config.rowKey
        : ((row, idx) => String(idx + 1));

      const validKeys = new Set(rows.map((row, idx) => text(rowKey(row, idx))));
      Array.from(checkedSet).forEach((key) => {
        if (!validKeys.has(text(key))) checkedSet.delete(key);
      });

      const shell = document.createElement("div");
      shell.className = "table-shell";
      shell.style.maxHeight = maxHeight;
      shell.style.height = maxHeight;

      if (rows.length === 0) {
        const empty = document.createElement("div");
        empty.className = "muted";
        empty.style.padding = "10px";
        empty.textContent = "No rows to display.";
        shell.appendChild(empty);
        container.appendChild(shell);
        return;
      }

      const table = document.createElement("table");
      table.className = "data-table";

      const thead = document.createElement("thead");
      const htr = document.createElement("tr");
      const hCheck = document.createElement("th");
      hCheck.className = "check-cell";
      htr.appendChild(hCheck);
      columns.forEach((col) => {
        const th = document.createElement("th");
        th.textContent = col;
        if (wrapCols.has(col)) th.classList.add("allele-name-col");
        htr.appendChild(th);
      });
      thead.appendChild(htr);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      rows.forEach((row, idx) => {
        const key = text(rowKey(row, idx));
        const tr = document.createElement("tr");

        const checkTd = document.createElement("td");
        checkTd.className = "check-cell";
        const cb = document.createElement("input");
        cb.type = "checkbox";
        cb.checked = checkedSet.has(key);
        cb.addEventListener("change", (ev) => {
          ev.stopPropagation();
          if (cb.checked) {
            checkedSet.add(key);
          } else {
            checkedSet.delete(key);
          }
          state.activeTable = config.tableKey;
          updateDownloadContext();
        });
        checkTd.appendChild(cb);
        tr.appendChild(checkTd);

        columns.forEach((col) => {
          const td = document.createElement("td");
          const fullText = text(row[col]);
          if (truncateTooltipCols.has(col) && fullText.length > truncateLength) {
            td.textContent = "...";
            td.title = fullText;
          } else {
            td.textContent = fullText;
          }
          if (multilineCols.has(col)) td.classList.add("multiline");
          if (wrapCols.has(col)) td.classList.add("allele-name-col");
          tr.appendChild(td);
        });

        tr.addEventListener("click", (ev) => {
          if (ev.target && ev.target.tagName === "INPUT") return;
          if (typeof config.onRowClick === "function") config.onRowClick(row, idx, key);
        });
        tr.addEventListener("dblclick", (ev) => {
          if (ev.target && ev.target.tagName === "INPUT") return;
          if (typeof config.onRowDblClick === "function") config.onRowDblClick(row, idx, key);
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);

      shell.appendChild(table);
      container.appendChild(shell);
    }

    function renderGeneSelectorTable() {
      const container = document.getElementById("gene-selector-table");
      container.innerHTML = "";

      const genes = getFilteredGenes();
      const shell = document.createElement("div");
      shell.className = "table-shell";
      shell.style.height = "100%";
      shell.style.maxHeight = "100%";

      if (genes.length === 0) {
        const empty = document.createElement("div");
        empty.className = "muted";
        empty.style.padding = "10px";
        empty.textContent = "No genes match the current search.";
        shell.appendChild(empty);
        container.appendChild(shell);
        return;
      }

      const table = document.createElement("table");
      table.className = "data-table";
      const thead = document.createElement("thead");
      const htr = document.createElement("tr");
      const th = document.createElement("th");
      th.textContent = "Gene";
      htr.appendChild(th);
      thead.appendChild(htr);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      genes.forEach((gene) => {
        const tr = document.createElement("tr");
        if (gene === state.selectedGene) tr.classList.add("selected-row");
        const td = document.createElement("td");
        td.textContent = gene;
        tr.appendChild(td);
        tr.addEventListener("click", () => {
          if (state.selectedGene !== gene) {
            state.selectedGene = gene;
            state.checked.detail.clear();
            renderAlleleMode();
            updateStateField();
            updateDownloadContext();
          }
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
      shell.appendChild(table);
      container.appendChild(shell);
    }

    function setActiveMainTableForMode() {
      if (state.mode === "GT") {
        state.activeTable = "gene";
      } else if (state.mode === "ET") {
        state.activeTable = "exon";
      } else {
        state.activeTable = "detail";
      }
    }

    function openRelationModal(type, title, rows) {
      state.modal.open = true;
      state.modal.type = type;
      state.modal.rows = rows.slice();
      state.modal.title = title;
      if (type === "variants") {
        state.checked.modalVariants.clear();
        state.activeTable = "modalVariants";
      } else if (type === "alleles") {
        state.checked.modalAlleles.clear();
        state.activeTable = "modalAlleles";
      }
      document.getElementById("modal-title").textContent = title;
      document.getElementById("relation-modal").classList.add("open");
      renderModalTable();
      updateDownloadContext();
    }

    function closeRelationModal() {
      state.modal.open = false;
      state.modal.type = null;
      state.modal.rows = [];
      state.modal.title = "";
      document.getElementById("relation-modal").classList.remove("open");
      setActiveMainTableForMode();
      updateDownloadContext();
    }

    function renderModalTable() {
      if (!state.modal.open) {
        document.getElementById("modal-table").innerHTML = "";
        return;
      }
      if (state.modal.type === "variants") {
        const displayRows = state.modal.rows.map((r) => normalizeRow(r, DATA.columns.variant_view));
        renderTable("modal-table", {
          rows: displayRows,
          columns: DATA.columns.variant_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.modalVariants,
          tableKey: "modalVariants",
          truncateTooltipCols: ["Ref_allele_curated", "Alt_allele_curated"],
          truncateLength: 10,
          maxHeight: "45vh"
        });
      } else if (state.modal.type === "alleles") {
        const displayRows = state.modal.rows.map((r) => normalizeRow(r, DATA.columns.allele_view));
        renderTable("modal-table", {
          rows: displayRows,
          columns: DATA.columns.allele_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.modalAlleles,
          tableKey: "modalAlleles",
          multilineCols: ["Variants"],
          maxHeight: "45vh"
        });
      }
    }

    function renderGeneMode() {
      const displayRows = DATA.gene_table.map((r) => normalizeRow(r, DATA.columns.gene_table));
      renderTable("gene-mode", {
        rows: displayRows,
        columns: DATA.columns.gene_table,
        rowKey: (row, idx) => String(idx + 1),
        checkedSet: state.checked.gene,
        tableKey: "gene",
        maxHeight: "100%"
      });
    }

    function renderExonMode() {
      const displayRows = DATA.exon_table.map((r) => normalizeRow(r, DATA.columns.exon_view));
      renderTable("exon-mode", {
        rows: displayRows,
        columns: DATA.columns.exon_view,
        rowKey: (row, idx) => String(idx + 1),
        checkedSet: state.checked.exon,
        tableKey: "exon",
        maxHeight: "100%"
      });
    }

    function renderAlleleMode() {
      renderGeneSelectorTable();

      const detailRows = detailBaseRows();
      if (state.granularity === "Allele") {
        const displayRows = detailRows.map((r) => normalizeRow(r, DATA.columns.allele_view));
        renderTable("detail-table", {
          rows: displayRows,
          columns: DATA.columns.allele_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.detail,
          tableKey: "detail",
          multilineCols: ["Variants"],
          maxHeight: "100%",
          onRowDblClick: (row) => {
            const alleleId = text(row.Allele_id);
            const linked = uniqueOrdered(DATA.bridge_by_allele[alleleId] || []);
            const rows = mapRowsByIds(DATA.variant_table, "Variant_id", linked);
            const label = text(row.Allele_name) || alleleId;
            openRelationModal("variants", "Variants associated with allele " + label, rows);
          }
        });
      } else {
        const displayRows = detailRows.map((r) => normalizeRow(r, DATA.columns.variant_view));
        renderTable("detail-table", {
          rows: displayRows,
          columns: DATA.columns.variant_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.detail,
          tableKey: "detail",
          truncateTooltipCols: ["Ref_allele_curated", "Alt_allele_curated"],
          truncateLength: 10,
          maxHeight: "100%",
          onRowClick: (row) => {
            const variantId = text(row.Variant_id);
            const linked = uniqueOrdered(DATA.bridge_by_variant[variantId] || []);
            const rows = mapRowsByIds(DATA.allele_table, "Allele_id", linked);
            const note = text(row.Nucleotide_change);
            const label = note ? (variantId + " (" + note + ")") : variantId;
            openRelationModal("alleles", "Alleles associated with variant " + label, rows);
          }
        });
      }
    }

    function refreshMainPanels() {
      applyLayoutState();
      setActiveButtons();
      const geneMode = document.getElementById("gene-mode");
      const exonMode = document.getElementById("exon-mode");
      const alleleMode = document.getElementById("allele-mode");
      if (state.mode === "GT") {
        geneMode.classList.remove("hidden");
        exonMode.classList.add("hidden");
        alleleMode.classList.add("hidden");
        renderGeneMode();
      } else if (state.mode === "ET") {
        geneMode.classList.add("hidden");
        exonMode.classList.remove("hidden");
        alleleMode.classList.add("hidden");
        renderExonMode();
      } else {
        geneMode.classList.add("hidden");
        exonMode.classList.add("hidden");
        alleleMode.classList.remove("hidden");
        renderAlleleMode();
      }
      updateStateField();
      updateDownloadContext();
    }

    function clearAllCheckedRows() {
      state.checked.gene.clear();
      state.checked.exon.clear();
      state.checked.detail.clear();
      state.checked.modalVariants.clear();
      state.checked.modalAlleles.clear();
      refreshMainPanels();
      if (state.modal.open) renderModalTable();
    }

    function toCsv(rows, columns) {
      const cols = columns.slice();
      const esc = (val) => {
        const s = text(val).replace(/"/g, '""');
        if (s.includes(",") || s.includes("\\n") || s.includes('"')) {
          return '"' + s + '"';
        }
        return s;
      };
      const lines = [];
      lines.push(cols.map(esc).join(","));
      rows.forEach((row) => {
        lines.push(cols.map((c) => esc(row[c])).join(","));
      });
      return lines.join("\\r\\n");
    }

    function saveCsv(rows, columns, filename) {
      const csvText = toCsv(rows, columns);
      const blob = new Blob([csvText], { type: "text/csv;charset=utf-8;" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    function checkedRowsByIndex(sourceRows, checkedSet) {
      const idx = Array.from(checkedSet)
        .map((v) => parseInt(v, 10))
        .filter((n) => Number.isInteger(n) && n >= 1 && n <= sourceRows.length)
        .sort((a, b) => a - b);
      return idx.map((n) => sourceRows[n - 1]);
    }

    function downloadMainSelection() {
      let rows = [];
      let columns = [];
      let filename = "checked_rows.csv";
      const dateTag = new Date().toISOString().slice(0, 10).replace(/-/g, "");

      if (state.activeTable === "gene") {
        columns = DATA.columns.gene_table;
        rows = checkedRowsByIndex(DATA.gene_table, state.checked.gene).map((r) => normalizeRow(r, columns));
        filename = "gene_table_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "exon") {
        columns = DATA.columns.exon_view;
        rows = checkedRowsByIndex(DATA.exon_table, state.checked.exon).map((r) => normalizeRow(r, columns));
        filename = "exon_table_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "detail") {
        if (state.granularity === "Allele") {
          columns = DATA.columns.allele_view;
          rows = checkedRowsByIndex(detailBaseRows(), state.checked.detail).map((r) => normalizeRow(r, columns));
          filename = "allele_detail_table_checked_rows_" + dateTag + ".csv";
        } else {
          columns = DATA.columns.variant_view;
          rows = checkedRowsByIndex(detailBaseRows(), state.checked.detail).map((r) => normalizeRow(r, columns));
          filename = "variant_detail_table_checked_rows_" + dateTag + ".csv";
        }
      } else if (state.activeTable === "modalVariants") {
        columns = DATA.columns.variant_view;
        rows = checkedRowsByIndex(state.modal.rows, state.checked.modalVariants).map((r) => normalizeRow(r, columns));
        filename = "variants_popup_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "modalAlleles") {
        columns = DATA.columns.allele_view;
        rows = checkedRowsByIndex(state.modal.rows, state.checked.modalAlleles).map((r) => normalizeRow(r, columns));
        filename = "alleles_popup_checked_rows_" + dateTag + ".csv";
      }

      saveCsv(rows, columns, filename);
    }

    function downloadModalSelection() {
      const dateTag = new Date().toISOString().slice(0, 10).replace(/-/g, "");
      if (state.modal.type === "variants") {
        const rows = checkedRowsByIndex(state.modal.rows, state.checked.modalVariants)
          .map((r) => normalizeRow(r, DATA.columns.variant_view));
        saveCsv(rows, DATA.columns.variant_view, "variants_popup_checked_rows_" + dateTag + ".csv");
      } else if (state.modal.type === "alleles") {
        const rows = checkedRowsByIndex(state.modal.rows, state.checked.modalAlleles)
          .map((r) => normalizeRow(r, DATA.columns.allele_view));
        saveCsv(rows, DATA.columns.allele_view, "alleles_popup_checked_rows_" + dateTag + ".csv");
      }
    }

    function clearModalSelection() {
      if (state.modal.type === "variants") {
        state.checked.modalVariants.clear();
      } else if (state.modal.type === "alleles") {
        state.checked.modalAlleles.clear();
      }
      renderModalTable();
      updateDownloadContext();
    }

    function renderLatestUpdates() {
      const root = document.getElementById("latest-updates");
      root.innerHTML = "";
      const ul = document.createElement("ul");
      const lines = [
        "Loaded Gene_table: " + DATA.counts.gene.toLocaleString() + " rows",
        "Loaded Allele_table: " + DATA.counts.allele.toLocaleString() + " rows",
        "Loaded Variant_table: " + DATA.counts.variant.toLocaleString() + " rows",
        "Loaded Exon_table: " + DATA.counts.exon.toLocaleString() + " rows",
        "Loaded Bridge_table: " + DATA.counts.bridge.toLocaleString() + " rows",
        "Use A1 buttons to switch between Gene, Allele, and Exon exploration modes.",
        "Allele rows: double-click to inspect linked variants.",
        "Variant rows: single-click to inspect linked alleles."
      ];
      lines.forEach((line) => {
        const li = document.createElement("li");
        li.textContent = line;
        ul.appendChild(li);
      });
      root.appendChild(ul);
    }

    function renderFeedbackStatus(entry) {
      const box = document.getElementById("feedback-status");
      if (!entry) {
        box.className = "last-feedback muted";
        box.textContent = "No submissions yet.";
        return;
      }
      box.className = "last-feedback";
      box.innerHTML = "";
      const title = document.createElement("strong");
      title.textContent = "Last submission";
      box.appendChild(title);
      const rows = [
        entry.submitted_at,
        "Email: " + (entry.email || "(blank)"),
        "State: " + (entry.state || ""),
        "Comment length: " + text(entry.comment || "").length
      ];
      rows.forEach((line) => {
        const div = document.createElement("div");
        div.textContent = line;
        box.appendChild(div);
      });
    }

    async function submitFeedback() {
      const payload = {
        state: document.getElementById("feedback-state").value,
        email: document.getElementById("feedback-email").value,
        comment: document.getElementById("feedback-comment").value
      };
      const res = await fetch("/api/feedback", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload)
      });
      const data = await res.json();
      if (data.ok) {
        document.getElementById("feedback-email").value = "";
        document.getElementById("feedback-comment").value = "";
        document.getElementById("feedback-flash").textContent = "Feedback submitted.";
        renderFeedbackStatus(data.last);
        setTimeout(() => {
          document.getElementById("feedback-flash").textContent = "";
        }, 2200);
      }
    }

    function bindEvents() {
      document.getElementById("mode-gene").addEventListener("click", () => {
        state.mode = "GT";
        setActiveMainTableForMode();
        refreshMainPanels();
      });
      document.getElementById("mode-allele").addEventListener("click", () => {
        state.mode = "AT";
        setActiveMainTableForMode();
        refreshMainPanels();
      });
      document.getElementById("mode-exon").addEventListener("click", () => {
        state.mode = "ET";
        setActiveMainTableForMode();
        refreshMainPanels();
      });
      document.getElementById("gran-allele").addEventListener("click", () => {
        if (state.granularity !== "Allele") {
          state.granularity = "Allele";
          state.checked.detail.clear();
          state.activeTable = "detail";
          refreshMainPanels();
        }
      });
      document.getElementById("gran-variant").addEventListener("click", () => {
        if (state.granularity !== "Variant") {
          state.granularity = "Variant";
          state.checked.detail.clear();
          state.activeTable = "detail";
          refreshMainPanels();
        }
      });
      document.getElementById("gene-search").addEventListener("input", (ev) => {
        state.geneSearch = ev.target.value || "";
        const prev = state.selectedGene;
        ensureSelectedGene();
        if (prev !== state.selectedGene) {
          state.checked.detail.clear();
        }
        renderGeneSelectorTable();
        renderAlleleMode();
        updateStateField();
        updateDownloadContext();
      });
      document.getElementById("table-search").addEventListener("input", (ev) => {
        state.tableSearch = ev.target.value || "";
        state.checked.detail.clear();
        renderAlleleMode();
        updateStateField();
        updateDownloadContext();
      });
      document.getElementById("toggle-panel-b").addEventListener("click", () => {
        state.panelBCollapsed = !state.panelBCollapsed;
        applyLayoutState();
      });
      document.getElementById("toggle-gene-selector").addEventListener("click", () => {
        state.alleleSelectorCollapsed = !state.alleleSelectorCollapsed;
        applyLayoutState();
      });
      document.getElementById("download-main-btn").addEventListener("click", downloadMainSelection);
      document.getElementById("clear-main-btn").addEventListener("click", clearAllCheckedRows);
      document.getElementById("feedback-submit").addEventListener("click", submitFeedback);

      document.getElementById("modal-close-btn").addEventListener("click", closeRelationModal);
      document.getElementById("modal-clear-btn").addEventListener("click", clearModalSelection);
      document.getElementById("modal-download-btn").addEventListener("click", downloadModalSelection);
      document.getElementById("relation-modal").addEventListener("click", (ev) => {
        if (ev.target.id === "relation-modal") closeRelationModal();
      });
    }

    async function init() {
      const response = await fetch("/api/init");
      DATA = await response.json();
      state.selectedGene = DATA.all_genes.length > 0 ? DATA.all_genes[0] : "";
      bindEvents();
      renderLatestUpdates();
      renderFeedbackStatus(null);
      refreshMainPanels();
    }

    init();
  </script>
</body>
</html>
"""


class AppHandler(BaseHTTPRequestHandler):
    def _send_json(self, payload, status=200):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self, html_text, status=200):
        body = html_text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path in ("/", "/index.html"):
            self._send_html(INDEX_HTML)
            return
        if parsed.path == "/api/init":
            self._send_json(INIT_PAYLOAD)
            return
        if parsed.path == "/api/feedback/latest":
            latest = FEEDBACK_LOG[-1] if FEEDBACK_LOG else None
            self._send_json({"ok": True, "last": latest, "count": len(FEEDBACK_LOG)})
            return
        self._send_json({"ok": False, "error": "Not found"}, status=404)

    def do_POST(self):
        parsed = urlparse(self.path)
        if parsed.path != "/api/feedback":
            self._send_json({"ok": False, "error": "Not found"}, status=404)
            return

        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length > 0 else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except Exception:
            self._send_json({"ok": False, "error": "Invalid JSON"}, status=400)
            return

        entry = {
            "submitted_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "state": clean_text(payload.get("state", "")),
            "email": clean_text(payload.get("email", "")),
            "comment": clean_text(payload.get("comment", "")),
        }
        FEEDBACK_LOG.append(entry)
        self._send_json({"ok": True, "last": entry, "count": len(FEEDBACK_LOG)})

    def log_message(self, format, *args):
        return


def main():
    host = "0.0.0.0"
    port = int(os.environ.get("PORT", "10000"))
    server = ThreadingHTTPServer((host, port), AppHandler)
    print(f"Serving app at http://{host}:{port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
