#!/usr/bin/env python3
import csv
import json
import os
import re
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen


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

ISBT_SYSTEMS_URL = "https://api-blooddatabase.isbtweb.org/system"
ISBT_ANTIGEN_URL = "https://api-blooddatabase.isbtweb.org/antigen/search?system_symbol="
ISBT_ALLELE_URL = "https://api-blooddatabase.isbtweb.org/allele/search?system_symbol="
ISBT_CACHE_PATH = ROOT_DIR / "isbt_cache.json"
ISBT_CACHE_MAX_DAYS = 30
ISBT_FLAG_COLS = ["sv_allele", "null_allele", "mod_allele", "partial_allele", "weak_allele", "el_allele"]
ISBT_GROUPED_VIEW_COLS = [
    "Group",
    "Allele_id",
    "ISBT_Allele",
    "Phenotype",
    "Gene_name",
    "DNA Change",
    "Exon/Intron",
    "HGVS Transcript",
    "sv_allele",
    "null_allele",
    "mod_allele",
    "partial_allele",
    "weak_allele",
    "el_allele",
]


def now_utc_epoch():
    return int(datetime.utcnow().timestamp())


def extract_records(payload):
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        for key in ("results", "data", "items", "systems", "antigens", "alleles", "variants"):
            value = payload.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
        return [payload]
    return []


def flatten_record(record, prefix=""):
    out = {}
    if not isinstance(record, dict):
        return out
    for key, value in record.items():
        clean_key = clean_text(key).strip()
        if not clean_key:
            continue
        full_key = f"{prefix}_{clean_key}" if prefix else clean_key
        if isinstance(value, dict):
            nested = flatten_record(value, full_key)
            out.update(nested)
        elif isinstance(value, list):
            if all(not isinstance(item, (dict, list)) for item in value):
                out[full_key] = ", ".join(clean_text(item) for item in value)
            else:
                out[full_key] = json.dumps(value, ensure_ascii=False)
        else:
            out[full_key] = clean_text(value)
    return out


def first_nonempty_ci(record, candidate_keys):
    if not isinstance(record, dict):
        return ""
    for wanted in candidate_keys:
        wanted_low = wanted.lower()
        for key, value in record.items():
            if clean_text(key).lower() == wanted_low:
                txt = clean_text(value).strip()
                if txt:
                    return txt
    return ""


def to_bool_flag_text(value):
    txt = clean_text(value).strip()
    low = txt.lower()
    if low in ("1", "true", "t", "yes", "y"):
        return "Yes"
    if low in ("0", "false", "f", "no", "n"):
        return "No"
    return txt


def extract_isbt_allele_phenotype(record):
    phenotype = first_nonempty_ci(record, ["Phenotype", "isbt_phenotype", "phenotype"])
    if phenotype:
        return phenotype

    raw_phenotypes = first_nonempty_ci(record, ["phenotypes"])
    if not raw_phenotypes:
        return ""

    try:
        parsed = json.loads(raw_phenotypes)
    except (TypeError, ValueError, json.JSONDecodeError):
        return ""

    if not isinstance(parsed, list):
        return ""

    labels = []
    seen = set()
    for item in parsed:
        if not isinstance(item, dict):
            continue
        antigen = item.get("antigen") if isinstance(item.get("antigen"), dict) else {}
        antigen_name = clean_text(antigen.get("display_name") or antigen.get("isbt_code")).strip()
        antigen_name = re.sub(r"<[^>]+>", "", antigen_name).strip()
        phenotype_name = clean_text(item.get("phenotype")).strip()
        if phenotype_name == "-":
            phenotype_name = ""
        label = f"{antigen_name} {phenotype_name}".strip()
        if label and label not in seen:
            seen.add(label)
            labels.append(label)
    return "; ".join(labels)


def normalize_dna_change(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return "-"
    if ":" in txt:
        txt = txt.split(":", 1)[1].strip()
    txt = re.sub(r"^(c\.[0-9_+\-*]+)[A-Za-z]+=$", r"\1=", txt)
    return txt or "-"


def normalize_exon_label(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return ""
    low = txt.lower()
    if low.startswith("exon"):
        rest = txt[4:].strip()
        return f"Exon {rest}" if rest else "Exon"
    if re.fullmatch(r"-?\d+", txt):
        return f"Exon {txt}"
    return f"Exon {txt}"


def normalize_intron_label(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return ""
    low = txt.lower()
    if low.startswith("intron"):
        rest = txt[6:].strip()
        return f"Intron {rest}" if rest else "Intron"
    return f"Intron {txt}"


def normalize_exon_intron(exon_value, intron_value):
    exon_label = normalize_exon_label(exon_value)
    intron_label = normalize_intron_label(intron_value)
    if exon_label and intron_label:
        return f"{exon_label} / {intron_label}"
    if exon_label:
        return exon_label
    if intron_label:
        return intron_label
    return "-"


def parse_c_position(dna_change):
    txt = clean_text(dna_change)
    match = re.search(r"c\.(-?\d+)", txt)
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return 10**9
    return 10**9


def parse_isbt_number(record):
    txt = first_nonempty_ci(record, ["isbt_number", "system_number", "number"])
    match = re.search(r"\d+", txt)
    if match:
        try:
            return int(match.group(0))
        except ValueError:
            return 10**9
    return 10**9


def fetch_json(url, timeout_sec=25):
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "geno-mock-python/1.0 (+https://github.com/RazaS/geno_mock_python)",
        },
    )
    with urlopen(request, timeout=timeout_sec) as response:
        body = response.read().decode("utf-8")
    return json.loads(body)


def variant_like_record(record):
    check_fields = [
        "input",
        "dna_change",
        "nucleotide_change",
        "hgvs",
        "hgvs_transcript",
        "exon",
        "intron",
        "variant_id",
    ]
    return bool(first_nonempty_ci(record, check_fields))


def dedupe_rows(rows, key_fields):
    seen = set()
    out = []
    for row in rows:
        key = tuple(clean_text(row.get(field, "")) for field in key_fields)
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def variant_matches_group(variant_row, selected_group):
    group = clean_text(selected_group).upper().strip()
    symbol = first_nonempty_ci(variant_row, ["system_symbol"]).upper()
    if group in ("RHD", "RHCE"):
        isbt_allele = first_nonempty_ci(variant_row, ["isbt_allele"]).upper()
        gene_name = first_nonempty_ci(variant_row, ["gene_name", "gene"]).upper()
        return isbt_allele.startswith(group + "*") or gene_name == group
    return symbol == group


def build_isbt_allele_meta_by_id(alleles_rows):
    allele_meta_by_id = {}
    for row in alleles_rows:
        allele_id = clean_text(row.get("allele_id", "")).strip()
        if not allele_id:
            continue
        meta = allele_meta_by_id.setdefault(allele_id, {})
        phenotype = extract_isbt_allele_phenotype(row)
        if phenotype and not clean_text(meta.get("Phenotype", "")).strip():
            meta["Phenotype"] = phenotype
        for flag_col in ISBT_FLAG_COLS:
            flag_value = to_bool_flag_text(first_nonempty_ci(row, [flag_col]))
            if flag_value and not clean_text(meta.get(flag_col, "")).strip():
                meta[flag_col] = flag_value
    return allele_meta_by_id


def normalize_isbt_variant_row(row, allele_meta_by_id, idx):
    out = dict(row)
    allele_id = first_nonempty_ci(out, ["allele_id"])
    if not allele_id:
        allele_id = f"missing:{idx}"
    out["allele_id"] = allele_id
    out["isbt_allele"] = first_nonempty_ci(out, ["isbt_allele", "allele_name"])
    out["gene_name"] = first_nonempty_ci(out, ["gene_name", "gene", "gene_symbol"])
    out["variant_id"] = first_nonempty_ci(out, ["variant_id", "id"]) or f"v:{idx}"

    dna_change = normalize_dna_change(first_nonempty_ci(out, ["input", "dna_change", "nucleotide_change", "hgvs", "hgvs_c"]))
    exon_value = first_nonempty_ci(out, ["exon", "exon_number", "exon_no"])
    intron_value = first_nonempty_ci(out, ["intron", "intron_number", "intron_no"])
    hgvs_transcript = clean_text(first_nonempty_ci(out, ["input", "hgvs_transcript", "hgvs"])).strip() or "-"

    out["DNA Change"] = dna_change
    out["Exon/Intron"] = normalize_exon_intron(exon_value, intron_value)
    out["HGVS Transcript"] = hgvs_transcript

    allele_meta = allele_meta_by_id.get(allele_id, {})
    phenotype = extract_isbt_allele_phenotype(out) or clean_text(allele_meta.get("Phenotype", "")).strip()
    out["Phenotype"] = phenotype

    for flag_col in ISBT_FLAG_COLS:
        if not first_nonempty_ci(out, [flag_col]):
            out[flag_col] = allele_meta.get(flag_col, "")
        out[flag_col] = to_bool_flag_text(out.get(flag_col, ""))

    return out


def hydrate_isbt_dataset(dataset):
    if not isinstance(dataset, dict):
        return empty_isbt_dataset("Invalid ISBT dataset.")

    alleles = [row for row in dataset.get("alleles", []) if isinstance(row, dict)]
    variants = [row for row in dataset.get("variants", []) if isinstance(row, dict)]
    antigens = [row for row in dataset.get("antigens", []) if isinstance(row, dict)]
    systems = [row for row in dataset.get("systems", []) if isinstance(row, dict)]

    for row in alleles:
        phenotype = extract_isbt_allele_phenotype(row)
        if phenotype:
            row["Phenotype"] = phenotype
        for flag_col in ISBT_FLAG_COLS:
            row[flag_col] = to_bool_flag_text(first_nonempty_ci(row, [flag_col]))

    allele_meta_by_id = build_isbt_allele_meta_by_id(alleles)
    normalized_variants = [
        normalize_isbt_variant_row(row, allele_meta_by_id, idx)
        for idx, row in enumerate(variants, start=1)
    ]
    normalized_variants = dedupe_rows(
        normalized_variants,
        ["system_symbol", "allele_id", "variant_id", "DNA Change", "Exon/Intron", "HGVS Transcript"],
    )
    grouped_rows = make_isbt_grouped_rows(normalized_variants)
    raw_export_columns = sorted({key for row in normalized_variants for key in row.keys() if not key.startswith("__")})

    dataset["systems"] = systems
    dataset["antigens"] = antigens
    dataset["alleles"] = alleles
    dataset["variants"] = normalized_variants
    dataset["grouped_rows"] = grouped_rows

    columns = dataset.setdefault("columns", {})
    columns["grouped_view"] = ISBT_GROUPED_VIEW_COLS
    columns["raw_export"] = raw_export_columns

    metadata = dataset.setdefault("metadata", {})
    row_counts = metadata.setdefault("row_counts", {})
    row_counts["systems"] = len(systems)
    row_counts["antigens"] = len(antigens)
    row_counts["alleles"] = len(alleles)
    row_counts["variants"] = len(normalized_variants)
    row_counts["grouped_rows"] = len(grouped_rows)

    return dataset


def make_isbt_grouped_rows(variants_rows):
    by_group = {}
    for row in variants_rows:
        symbol = first_nonempty_ci(row, ["system_symbol"]).upper()
        if not symbol:
            continue
        target_groups = [symbol]
        if symbol == "RH":
            target_groups = ["RHD", "RHCE"]
        for group in target_groups:
            if variant_matches_group(row, group):
                by_group.setdefault(group, []).append(row)

    grouped_rows = []
    for group, group_rows in by_group.items():
        by_allele = {}
        for row in group_rows:
            allele_id = first_nonempty_ci(row, ["allele_id"])
            isbt_allele = first_nonempty_ci(row, ["isbt_allele"])
            if not allele_id:
                allele_id = isbt_allele or f"{group}:{len(by_allele) + 1}"
            row_key = f"{group}|{allele_id}"
            by_allele.setdefault(row_key, []).append(row)

        for row_key, rows in by_allele.items():
            rows_sorted = sorted(
                rows,
                key=lambda r: (
                    parse_c_position(first_nonempty_ci(r, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"])),
                    first_nonempty_ci(r, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"]),
                ),
            )

            dna_lines = []
            exon_intron_lines = []
            hgvs_lines = []
            for raw in rows_sorted:
                dna_lines.append(first_nonempty_ci(raw, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"]) or "-")
                exon_intron_lines.append(first_nonempty_ci(raw, ["Exon/Intron", "exon_intron"]) or "-")
                hgvs_lines.append(first_nonempty_ci(raw, ["HGVS Transcript", "hgvs_transcript", "input", "hgvs"]) or "-")

            first_row = rows_sorted[0]
            grouped_row = {
                "Group": group,
                "Allele_id": first_nonempty_ci(first_row, ["allele_id"]),
                "ISBT_Allele": first_nonempty_ci(first_row, ["isbt_allele", "allele_name"]),
                "Phenotype": first_nonempty_ci(first_row, ["Phenotype", "isbt_phenotype", "phenotype"]),
                "Gene_name": first_nonempty_ci(first_row, ["gene_name", "gene"]),
                "DNA Change": "\n".join(dna_lines),
                "Exon/Intron": "\n".join(exon_intron_lines),
                "HGVS Transcript": "\n".join(hgvs_lines),
                "__group": group,
                "__row_key": row_key,
                "__raw_variant_rows": rows_sorted,
            }
            for flag_col in ISBT_FLAG_COLS:
                grouped_row[flag_col] = first_nonempty_ci(first_row, [flag_col])
            grouped_rows.append(grouped_row)

    grouped_rows.sort(
        key=lambda row: (
            clean_text(row.get("Group", "")),
            clean_text(row.get("ISBT_Allele", "")),
            clean_text(row.get("Allele_id", "")),
        )
    )
    return grouped_rows


def build_isbt_dataset():
    systems_payload = fetch_json(ISBT_SYSTEMS_URL)
    systems_raw = extract_records(systems_payload)

    systems = []
    symbols = []
    for record in systems_raw:
        flat = flatten_record(record)
        symbol = first_nonempty_ci(flat, ["system_symbol", "symbol"]).upper()
        if not symbol:
            continue
        flat["system_symbol"] = symbol
        systems.append(flat)
        symbols.append(symbol)

    systems = sorted(
        dedupe_rows(systems, ["system_symbol"]),
        key=lambda row: (parse_isbt_number(row), clean_text(row.get("system_symbol", ""))),
    )
    symbols = [clean_text(row.get("system_symbol", "")).upper() for row in systems]

    antigens = []
    alleles = []
    variants = []
    pull_errors = []

    for symbol in symbols:
        if not symbol:
            continue

        try:
            antigen_payload = fetch_json(ISBT_ANTIGEN_URL + quote(symbol))
            antigen_rows = extract_records(antigen_payload)
        except Exception as exc:
            pull_errors.append(f"antigen:{symbol}:{clean_text(exc)}")
            antigen_rows = []
        for record in antigen_rows:
            flat = flatten_record(record)
            flat["system_symbol"] = symbol
            antigens.append(flat)

        try:
            allele_payload = fetch_json(ISBT_ALLELE_URL + quote(symbol))
            allele_rows = extract_records(allele_payload)
        except Exception as exc:
            pull_errors.append(f"allele:{symbol}:{clean_text(exc)}")
            allele_rows = []
        for record in allele_rows:
            if not isinstance(record, dict):
                continue
            flat_allele = flatten_record(record)
            flat_allele["system_symbol"] = symbol

            allele_id = first_nonempty_ci(flat_allele, ["allele_id", "id"])
            if not allele_id:
                allele_id = f"{symbol}:{len(alleles) + 1}"
            flat_allele["allele_id"] = allele_id
            flat_allele["isbt_allele"] = first_nonempty_ci(flat_allele, ["isbt_allele", "allele_name", "name"])
            flat_allele["gene_name"] = first_nonempty_ci(flat_allele, ["gene_name", "gene", "gene_symbol"])
            flat_allele["Phenotype"] = extract_isbt_allele_phenotype(flat_allele)

            for flag_col in ISBT_FLAG_COLS:
                flat_allele[flag_col] = to_bool_flag_text(first_nonempty_ci(flat_allele, [flag_col]))

            alleles.append(flat_allele)

            nested_variants = []
            for key in ("variants", "variant", "variant_list", "allele_variants"):
                value = record.get(key)
                if isinstance(value, list):
                    nested_variants = [item for item in value if isinstance(item, dict)]
                    if nested_variants:
                        break

            if nested_variants:
                for variant_record in nested_variants:
                    flat_variant = flatten_record(variant_record)
                    flat_variant["system_symbol"] = symbol
                    flat_variant["allele_id"] = allele_id
                    if not first_nonempty_ci(flat_variant, ["isbt_allele"]):
                        flat_variant["isbt_allele"] = flat_allele.get("isbt_allele", "")
                    if not first_nonempty_ci(flat_variant, ["gene_name", "gene"]):
                        flat_variant["gene_name"] = flat_allele.get("gene_name", "")
                    variants.append(flat_variant)
            elif variant_like_record(flat_allele):
                variants.append(dict(flat_allele))

    allele_by_id = {}
    for row in alleles:
        allele_id = clean_text(row.get("allele_id", ""))
        if not allele_id:
            continue
        if allele_id not in allele_by_id:
            allele_by_id[allele_id] = dict(row)
            continue
        existing = allele_by_id[allele_id]
        for key, value in row.items():
            if value and not clean_text(existing.get(key, "")).strip():
                existing[key] = value
    alleles = list(allele_by_id.values())
    antigens = dedupe_rows(antigens, ["system_symbol", "antigen_id", "id", "name"])

    group_symbols = [clean_text(row.get("system_symbol", "")).upper() for row in systems if clean_text(row.get("system_symbol", ""))]
    groups = []
    for symbol in group_symbols:
        if symbol == "RH":
            groups.extend(["RHD", "RHCE"])
        else:
            groups.append(symbol)
    groups = unique_keep_order(groups)

    meta = {
        "updated_at": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "updated_epoch": now_utc_epoch(),
        "source_urls": [ISBT_SYSTEMS_URL, ISBT_ANTIGEN_URL + "<SYMBOL>", ISBT_ALLELE_URL + "<SYMBOL>"],
        "row_counts": {
            "systems": len(systems),
            "antigens": len(antigens),
            "alleles": len(alleles),
            "variants": len(variants),
            "grouped_rows": 0,
        },
        "pull_errors": pull_errors,
    }

    return hydrate_isbt_dataset({
        "metadata": meta,
        "systems": systems,
        "antigens": antigens,
        "alleles": alleles,
        "variants": variants,
        "groups": groups,
        "grouped_rows": [],
        "columns": {
            "grouped_view": ISBT_GROUPED_VIEW_COLS,
            "raw_export": [],
        },
    })


def empty_isbt_dataset(error_message=""):
    return {
        "metadata": {
            "updated_at": "Not available",
            "updated_epoch": 0,
            "source_urls": [ISBT_SYSTEMS_URL, ISBT_ANTIGEN_URL + "<SYMBOL>", ISBT_ALLELE_URL + "<SYMBOL>"],
            "row_counts": {"systems": 0, "antigens": 0, "alleles": 0, "variants": 0, "grouped_rows": 0},
            "pull_errors": [],
            "error": clean_text(error_message),
        },
        "systems": [],
        "antigens": [],
        "alleles": [],
        "variants": [],
        "groups": [],
        "grouped_rows": [],
        "columns": {"grouped_view": ISBT_GROUPED_VIEW_COLS, "raw_export": []},
    }


def load_isbt_dataset():
    cached = None
    if ISBT_CACHE_PATH.exists():
        try:
            with ISBT_CACHE_PATH.open("r", encoding="utf-8") as handle:
                cached = json.load(handle)
        except Exception:
            cached = None

    if cached:
        updated_epoch = int(cached.get("metadata", {}).get("updated_epoch", 0) or 0)
        if updated_epoch > 0:
            age_days = (now_utc_epoch() - updated_epoch) / 86400.0
            if age_days <= ISBT_CACHE_MAX_DAYS:
                return hydrate_isbt_dataset(cached)

    try:
        fresh = build_isbt_dataset()
    except (URLError, TimeoutError, json.JSONDecodeError, OSError, ValueError) as exc:
        if cached:
            cached_meta = cached.setdefault("metadata", {})
            cached_meta["refresh_error"] = clean_text(exc)
            return hydrate_isbt_dataset(cached)
        return empty_isbt_dataset(clean_text(exc))

    # Cache write failures should not discard freshly pulled data.
    try:
        with ISBT_CACHE_PATH.open("w", encoding="utf-8") as handle:
            json.dump(fresh, handle, ensure_ascii=False)
    except OSError as exc:
        fresh_meta = fresh.setdefault("metadata", {})
        fresh_meta["cache_write_error"] = clean_text(exc)
    return hydrate_isbt_dataset(fresh)


isbt_data = load_isbt_dataset()

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
        "isbt_grouped": len(isbt_data.get("grouped_rows", [])),
    },
    "columns": {
        "gene_table": gene_view_cols,
        "allele_view": allele_view_cols,
        "variant_view": variant_view_cols,
        "exon_view": exon_view_cols,
    },
    "isbt": isbt_data,
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
    .a2-isbt-layout {
      display: flex;
      gap: 12px;
      flex: 1 1 auto;
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .a2-isbt-a {
      flex: 0 0 24%;
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .a2-isbt-b {
      flex: 1 1 76%;
      min-width: 0;
      height: 100%;
      min-height: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .isbt-controls-grid {
      display: grid;
      grid-template-columns: minmax(220px, 1fr) minmax(220px, 1fr) 190px;
      gap: 8px;
      align-items: end;
    }
    .a2-isbt-layout.group-selector-collapsed .a2-isbt-a {
      display: none;
    }
    .a2-isbt-layout.group-selector-collapsed .a2-isbt-b {
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
      .a2-allele-a, .a2-allele-b, .a2-isbt-a, .a2-isbt-b {
        max-height: none;
      }
      .detail-controls-grid, .isbt-controls-grid {
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
            <button id="mode-isbt" class="toggle-btn">ISBT Table</button>
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
          <div id="isbt-mode" class="hidden">
            <div id="isbt-layout" class="a2-isbt-layout">
              <div class="a2-isbt-a">
                <div class="gene-search-wrap">
                  <label for="isbt-group-search">Search groups</label>
                  <input id="isbt-group-search" class="input" type="text" placeholder="Type to filter ISBT groups" />
                </div>
                <div class="gene-list-wrap">
                  <div id="isbt-group-selector-table"></div>
                </div>
              </div>
              <div class="a2-isbt-b">
                <div class="detail-controls">
                  <div class="isbt-controls-grid">
                    <div>
                      <label for="isbt-table-search">Search ISBT table</label>
                      <input id="isbt-table-search" class="input" type="text" placeholder="Filter grouped allele rows" />
                    </div>
                    <div id="isbt-last-updated" class="download-note"></div>
                    <button id="toggle-isbt-selector" class="btn">Hide Group Selector</button>
                  </div>
                </div>
                <div class="detail-table-wrap">
                  <div id="isbt-detail-table"></div>
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
      selectedIsbtGroup: "",
      isbtGroupSearch: "",
      isbtTableSearch: "",
      panelBCollapsed: false,
      alleleSelectorCollapsed: false,
      isbtSelectorCollapsed: false,
      activeTable: "gene",
      checked: {
        gene: new Set(),
        exon: new Set(),
        detail: new Set(),
        isbtDetail: new Set(),
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

    function getFilteredIsbtGroups() {
      const groups = (DATA.isbt && Array.isArray(DATA.isbt.groups)) ? DATA.isbt.groups : [];
      if (!state.isbtGroupSearch.trim()) return groups.slice();
      return groups.filter((group) => matchesSearchGroups(group, state.isbtGroupSearch));
    }

    function ensureSelectedIsbtGroup() {
      const groups = getFilteredIsbtGroups();
      if (groups.length === 0) {
        state.selectedIsbtGroup = "";
        return;
      }
      if (!groups.includes(state.selectedIsbtGroup)) {
        state.selectedIsbtGroup = groups[0];
        state.checked.isbtDetail.clear();
      }
    }

    function isbtDetailBaseRows() {
      if (!state.selectedIsbtGroup || !DATA.isbt || !Array.isArray(DATA.isbt.grouped_rows)) return [];
      const rows = DATA.isbt.grouped_rows.filter((row) => text(row.__group) === state.selectedIsbtGroup);
      if (!state.isbtTableSearch.trim()) return rows;
      return rows.filter((row) => {
        const joined = Object.values(row)
          .filter((value) => !Array.isArray(value) && typeof value !== "object")
          .map(text)
          .join(" || ");
        return matchesSearchGroups(joined, state.isbtTableSearch);
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
      document.getElementById("mode-isbt").classList.toggle("active", state.mode === "IT");
      document.getElementById("gran-allele").classList.toggle("active", state.granularity === "Allele");
      document.getElementById("gran-variant").classList.toggle("active", state.granularity === "Variant");
    }

    function currentVisibleStateText() {
      const modeLabel = state.mode === "GT"
        ? "Gene"
        : (state.mode === "AT" ? "Allele" : (state.mode === "IT" ? "ISBT" : "Exon"));
      const geneLabel = state.mode === "AT"
        ? (state.selectedGene || "")
        : (state.mode === "IT" ? (state.selectedIsbtGroup || "") : "N/A");
      const viewLabel = state.mode === "AT"
        ? state.granularity + "_Table"
        : (state.mode === "ET" ? "Exon_Table" : (state.mode === "IT" ? "ISBT_Grouped_Table" : "N/A"));
      return "Mode: " + modeLabel
        + " | Gene: " + geneLabel
        + " | B2 View: " + viewLabel
        + " | GeneSearch: " + state.geneSearch
        + " | TableSearch: " + state.tableSearch
        + " | ISBTGroupSearch: " + state.isbtGroupSearch
        + " | ISBTTableSearch: " + state.isbtTableSearch;
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
      if (state.activeTable === "isbtDetail") {
        return { label: "ISBT Grouped Table", checked: state.checked.isbtDetail };
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
      const isbtLayout = document.getElementById("isbt-layout");
      const isbtBtn = document.getElementById("toggle-isbt-selector");

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
      if (isbtLayout) {
        isbtLayout.classList.toggle("group-selector-collapsed", state.isbtSelectorCollapsed);
      }
      if (isbtBtn) {
        isbtBtn.textContent = state.isbtSelectorCollapsed ? "Show Group Selector" : "Hide Group Selector";
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

    function renderIsbtGroupSelectorTable() {
      const container = document.getElementById("isbt-group-selector-table");
      container.innerHTML = "";

      const groups = getFilteredIsbtGroups();
      const shell = document.createElement("div");
      shell.className = "table-shell";
      shell.style.height = "100%";
      shell.style.maxHeight = "100%";

      if (groups.length === 0) {
        const empty = document.createElement("div");
        empty.className = "muted";
        empty.style.padding = "10px";
        const isbtMeta = (DATA.isbt && DATA.isbt.metadata) ? DATA.isbt.metadata : {};
        const loadError = text(isbtMeta.error || isbtMeta.refresh_error || "");
        empty.textContent = loadError
          ? ("ISBT data unavailable: " + loadError)
          : "No ISBT groups match the current search.";
        shell.appendChild(empty);
        container.appendChild(shell);
        return;
      }

      const table = document.createElement("table");
      table.className = "data-table";
      const thead = document.createElement("thead");
      const htr = document.createElement("tr");
      const th = document.createElement("th");
      th.textContent = "Group";
      htr.appendChild(th);
      thead.appendChild(htr);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      groups.forEach((group) => {
        const tr = document.createElement("tr");
        if (group === state.selectedIsbtGroup) tr.classList.add("selected-row");
        const td = document.createElement("td");
        td.textContent = group;
        tr.appendChild(td);
        tr.addEventListener("click", () => {
          if (state.selectedIsbtGroup !== group) {
            state.selectedIsbtGroup = group;
            state.checked.isbtDetail.clear();
            renderIsbtMode();
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
      } else if (state.mode === "IT") {
        state.activeTable = "isbtDetail";
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

    function renderIsbtMode() {
      ensureSelectedIsbtGroup();
      renderIsbtGroupSelectorTable();

      const lastUpdated = document.getElementById("isbt-last-updated");
      const isbtMeta = (DATA.isbt && DATA.isbt.metadata) ? DATA.isbt.metadata : {};
      const updatedText = text(isbtMeta.updated_at || "Not available");
      const groupedCount = Number((isbtMeta.row_counts || {}).grouped_rows || 0).toLocaleString();
      const metaNotes = [];
      if (isbtMeta.refresh_error) metaNotes.push("refresh error");
      if (isbtMeta.cache_write_error) metaNotes.push("cache not persisted");
      if (lastUpdated) {
        const noteSuffix = metaNotes.length > 0 ? (" | Notes: " + metaNotes.join(", ")) : "";
        lastUpdated.textContent = "Last pulled: " + updatedText + " | Grouped rows: " + groupedCount + noteSuffix;
      }

      const detailRows = isbtDetailBaseRows();
      const columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.grouped_view))
        ? DATA.isbt.columns.grouped_view
        : [];
      const displayRows = detailRows.map((row) => normalizeRow(row, columns));

      renderTable("isbt-detail-table", {
        rows: displayRows,
        columns: columns,
        rowKey: (row, idx) => String(idx + 1),
        checkedSet: state.checked.isbtDetail,
        tableKey: "isbtDetail",
        multilineCols: ["DNA Change", "Exon/Intron", "HGVS Transcript"],
        maxHeight: "100%"
      });
    }

    function refreshMainPanels() {
      applyLayoutState();
      setActiveButtons();
      const geneMode = document.getElementById("gene-mode");
      const exonMode = document.getElementById("exon-mode");
      const alleleMode = document.getElementById("allele-mode");
      const isbtMode = document.getElementById("isbt-mode");
      if (state.mode === "GT") {
        geneMode.classList.remove("hidden");
        exonMode.classList.add("hidden");
        alleleMode.classList.add("hidden");
        isbtMode.classList.add("hidden");
        renderGeneMode();
      } else if (state.mode === "ET") {
        geneMode.classList.add("hidden");
        exonMode.classList.remove("hidden");
        alleleMode.classList.add("hidden");
        isbtMode.classList.add("hidden");
        renderExonMode();
      } else if (state.mode === "IT") {
        geneMode.classList.add("hidden");
        exonMode.classList.add("hidden");
        alleleMode.classList.add("hidden");
        isbtMode.classList.remove("hidden");
        renderIsbtMode();
      } else {
        geneMode.classList.add("hidden");
        exonMode.classList.add("hidden");
        alleleMode.classList.remove("hidden");
        isbtMode.classList.add("hidden");
        renderAlleleMode();
      }
      updateStateField();
      updateDownloadContext();
    }

    function clearAllCheckedRows() {
      state.checked.gene.clear();
      state.checked.exon.clear();
      state.checked.detail.clear();
      state.checked.isbtDetail.clear();
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

    function expandIsbtGroupedRows(groupedRows) {
      const rawRows = [];
      groupedRows.forEach((row) => {
        const variants = Array.isArray(row.__raw_variant_rows) ? row.__raw_variant_rows : [];
        variants.forEach((raw) => rawRows.push(raw));
      });
      return rawRows;
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
      } else if (state.activeTable === "isbtDetail") {
        columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.raw_export))
          ? DATA.isbt.columns.raw_export
          : [];
        const groupedRows = checkedRowsByIndex(isbtDetailBaseRows(), state.checked.isbtDetail);
        const rawRows = expandIsbtGroupedRows(groupedRows);
        rows = rawRows.map((r) => normalizeRow(r, columns));
        filename = "isbt_grouped_selection_expanded_" + dateTag + ".csv";
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
        "Loaded ISBT grouped rows: " + DATA.counts.isbt_grouped.toLocaleString() + " rows",
        "Use A1 buttons to switch between Gene, Allele, and Exon exploration modes.",
        "ISBT mode groups rows by allele and can export expanded raw variants.",
        "Allele rows: double-click to inspect linked variants.",
        "Variant rows: single-click to inspect linked alleles."
      ];
      if (DATA.isbt && DATA.isbt.metadata) {
        lines.push("ISBT last pulled: " + text(DATA.isbt.metadata.updated_at || "Not available"));
        const pullErrors = Array.isArray(DATA.isbt.metadata.pull_errors) ? DATA.isbt.metadata.pull_errors : [];
        if (pullErrors.length > 0) {
          lines.push("ISBT refresh warnings: " + pullErrors.length.toLocaleString() + " symbol-level fetch errors.");
        }
        if (DATA.isbt.metadata.refresh_error) {
          lines.push("ISBT refresh error: " + text(DATA.isbt.metadata.refresh_error));
        }
        if (DATA.isbt.metadata.error) {
          lines.push("ISBT load error: " + text(DATA.isbt.metadata.error));
        }
        if (DATA.isbt.metadata.cache_write_error) {
          lines.push("ISBT cache write warning: " + text(DATA.isbt.metadata.cache_write_error));
        }
      }
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
      document.getElementById("mode-isbt").addEventListener("click", () => {
        state.mode = "IT";
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
      document.getElementById("isbt-group-search").addEventListener("input", (ev) => {
        state.isbtGroupSearch = ev.target.value || "";
        ensureSelectedIsbtGroup();
        state.checked.isbtDetail.clear();
        renderIsbtMode();
        updateStateField();
        updateDownloadContext();
      });
      document.getElementById("isbt-table-search").addEventListener("input", (ev) => {
        state.isbtTableSearch = ev.target.value || "";
        state.checked.isbtDetail.clear();
        renderIsbtMode();
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
      document.getElementById("toggle-isbt-selector").addEventListener("click", () => {
        state.isbtSelectorCollapsed = !state.isbtSelectorCollapsed;
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
      state.selectedIsbtGroup = (DATA.isbt && Array.isArray(DATA.isbt.groups) && DATA.isbt.groups.length > 0)
        ? DATA.isbt.groups[0]
        : "";
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
