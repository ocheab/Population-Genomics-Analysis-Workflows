#!/usr/bin/env python3
"""Read *one* native Zarr genotype chunk for Pf8 PI4K (Southwest-5 pilot).

Safe staged retrieval:
 - reuses the verified Step-14 site-inventory and the five state sample lists;
 - downloads one source genotype chunk ONLY (one native variant chunk and
   one native sample chunk), never a whole chromosome or 282-sample matrix;
 - stores a small NPZ subset and provenance report; later work may reuse it.

This pilot does NOT create a complete 282-sample dataset, and the resulting
sample genotype calls are *unfiltered*. Network byte usage is not measured.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import malariagen_data

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://pf8-release.cog.sanger.ac.uk/snp-only"
GENE, CHROM, START, END = "PF3D7_0509800", "Pf3D7_05_v3", 409403, 415758
STATES = {"Oyo": 71, "Ogun": 38, "Lagos": 129, "Osun": 15, "Ondo": 29}
LISTDIR = ROOT / "results/Southwest5/sample_lists"
RES = ROOT / "results/Southwest5"
CACHE = ROOT / "variants/Southwest5/Sanger_SNP_chunk_cache"
INVENTORY = RES / "PI4K_Sanger_SNP_site_inventory.tsv"
PLAN = RES / "PI4K_Sanger_SNP_retrieval_plan.json"
REPORT = RES / "PI4K_Sanger_SNP_pilot_report.json"


def fail(message):
    raise RuntimeError(message)


def read_ids(file):
    items = [line.strip() for line in file.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if len(items) != len(set(items)):
        fail(f"Duplicate IDs in {file}")
    return items


def str_id(value):
    if isinstance(value, bytes):
        return value.decode("utf-8").rstrip("\x00")
    return str(value).rstrip("\x00")


def get_chunk(index, lengths):
    return int(np.searchsorted(np.cumsum(lengths), index, side="right"))


def run():
    print("=== STEP 15: Pf8 PI4K ONE-CHUNK SNP PILOT ===", flush=True)
    print("No full VCF download; only ONE native genotype chunk will be read.", flush=True)
    if not INVENTORY.is_file() or not PLAN.is_file():
        fail("Missing Step 14 site inventory or JSON plan; run Step 14 first")
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    if (plan.get("selected_sample_ids") != 282
        or plan.get("selected_snp_positions") != 2471
        or plan.get("selected_filter_pass_sites") != 1824
        or plan.get("genotype_matrix_downloaded") is not False):
        fail("Step 14 plan differs from verified 282-sample, 2471-SNP scope")

    all_ids = []
    states_by_id = {}
    for state, expected in STATES.items():
        ids = read_ids(LISTDIR / f"{state}_sample_IDs.txt")
        if len(ids) != expected:
            fail(f"{state}: expected {expected} sample IDs but found {len(ids)}")
        all_ids.extend(ids)
        for sid in ids:
            if sid in states_by_id:
                fail(f"Duplicate ID across states: {sid}")
            states_by_id[sid] = state
    if len(all_ids) != 282:
        fail("Expected 282 samples in total")

    sites = list(csv.DictReader(INVENTORY.open(encoding="utf-8-sig", newline=""), delimiter="\t"))
    if len(sites) != 2471:
        fail(f"Expected 2471 SNP positions, found {len(sites)}")
    if any(s["chromosome"] != CHROM or not (START <= int(s["position_1based"]) <= END) for s in sites):
        fail("Site inventory has a chromosome or region mismatch")
    v = np.array([int(s["global_variant_index_0based"]) for s in sites], dtype=np.int64)
    positions = np.array([int(s["position_1based"]) for s in sites], dtype=np.int64)
    passed = np.array([s["filter_pass"].strip().lower() == "true" for s in sites], dtype=bool)
    if len(np.unique(v)) != len(v) or not np.all(np.diff(v) > 0) or not np.all(np.diff(positions) >= 0):
        fail("Site variant indices are not strictly increasing and unique")
    if passed.sum() != 1824:
        fail("Site PASS count does not match the verified Step 14 inventory")

    print("Opening public SNP-only Pf8 Zarr metadata...", flush=True)
    pf8 = malariagen_data.Pf8(BASE)
    ds = pf8.variant_calls()
    gt = ds["call_genotype"].data
    if gt.ndim != 3 or gt.chunks is None:
        fail("Unexpected genotype array format")
    cv, cs, cp = gt.chunks
    if len(cp) != 1 or sum(cp) != 2:
        fail(f"Unexpected ploidy chunks: {cp}")
    v_chunks = {get_chunk(int(x), cv) for x in v}
    if len(v_chunks) != 1:
        fail(f"Expected PI4K positions in one variant chunk, got {len(v_chunks)}")

    source_ids = [str_id(s) for s in ds["sample_id"].data.compute()]
    if len(set(source_ids)) != len(source_ids):
        fail("Duplicate source sample IDs")
    index_by_id = {sid: i for i, sid in enumerate(source_ids)}
    missing = set(all_ids) - set(index_by_id)
    if missing:
        fail(f"Not found in SNP callset: {list(sorted(missing))[:5]}")
    sample_groups = {}
    for sid in all_ids:
        i = index_by_id[sid]
        k = get_chunk(i, cs)
        sample_groups.setdefault(k, []).append((sid, i))
    if len(sample_groups) != 13:
        fail(f"Expected 13 sample chunks from Step 14, found {len(sample_groups)}")

    # Use the chunk with the MOST target samples for a representative pilot.
    pilot_chunk = min(sample_groups, key=lambda k: (-len(sample_groups[k]), k))
    selected = sorted(sample_groups[pilot_chunk], key=lambda x: x[1])
    s0 = sum(cs[:pilot_chunk])
    s1 = s0 + cs[pilot_chunk]
    var_chunk = next(iter(v_chunks))
    var0 = sum(cv[:var_chunk])
    var1 = var0 + cv[var_chunk]
    if not (var0 <= v[0] <= v[-1] < var1):
        fail("Variant slice not covered by one native chunk")

    CACHE.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE / f"PI4K_native_sample_chunk_{pilot_chunk:03d}.npz"
    cache_sha = hashlib.sha256(INVENTORY.read_bytes()).hexdigest()
    sample_ids_selected = [name for name, _ in selected]
    n_native_uncompressed = int(cv[var_chunk] * cs[pilot_chunk] * sum(cp) * np.dtype(gt.dtype).itemsize)
    print(f"Target interval: {CHROM}:{START}-{END}", flush=True)
    print(f"Sites: {len(v)} (PASS={passed.sum()})", flush=True)
    print(f"Pilot chunk index: {pilot_chunk}; target samples in this chunk: {len(selected)}", flush=True)
    print(f"Source chunk UNCOMPRESSED size: {n_native_uncompressed / (1024**2):.2f} MiB", flush=True)
    print("Actual compressed network bytes are UNKNOWN.", flush=True)
    print("Reading exactly one source genotype chunk...", flush=True)
    t0 = time.monotonic()

    # Important: slice [sites, full native 100-sample chunk, ploidy].
    # This causes only one genotype Zarr chunk to be fetched and decompressed.
    raw = gt[v[0]:v[-1] + 1, s0:s1, :].compute(scheduler="single-threaded")
    if raw.ndim != 3 or raw.shape[2] != 2:
        fail(f"Unexpected genotype shape: {raw.shape}")
    row_offsets = v - v[0]
    col_offsets = np.array([i - s0 for _, i in selected], dtype=np.int64)
    subset = np.asarray(raw[row_offsets][:, col_offsets, :], dtype=np.int8)
    if subset.shape != (len(v), len(selected), 2):
        fail(f"Unexpected selected genotype shape: {subset.shape}")
    if np.any(subset < -1):
        fail("Unexpected genotype allele code below -1")
    elapsed = time.monotonic() - t0

    temp_file = cache_file.with_suffix(".npz.partial")
    with temp_file.open("wb") as fh:
        np.savez_compressed(
            fh, genotype=subset,
            variant_index_0based=v, position_1based=positions,
            filter_pass=passed,
            sample_id=np.array(sample_ids_selected, dtype="U32"),
            sample_index_0based=np.array([i for _, i in selected], dtype=np.int64),
            sample_chunk_index=np.array(pilot_chunk),
            site_inventory_sha256=np.array(cache_sha),
        )
    temp_file.replace(cache_file)
    report = {
        "status": "ONE_CHUNK_PILOT_SUCCESS",
        "source": BASE,
        "gene": GENE,
        "region": f"{CHROM}:{START}-{END}",
        "pilot_source_sample_chunk_index": pilot_chunk,
        "source_chunk_uncompressed_bytes": n_native_uncompressed,
        "pilot_target_samples": len(selected),
        "pilot_target_sample_by_state": {s: sum(states_by_id[x] == s for x in sample_ids_selected) for s in STATES},
        "sites_total": len(v),
        "sites_Pf8_PASS": int(passed.sum()),
        "sample_genotype_shape": list(subset.shape),
        "missing_allele_proportion_in_pilot": float((subset < 0).mean()),
        "elapsed_seconds": round(elapsed, 2),
        "compressed_network_download_bytes": "not measured",
        "cache_file": str(cache_file),
        "cache_file_bytes": cache_file.stat().st_size,
        "site_inventory_sha256": cache_sha,
        "warning": "PILOT ONLY: not a complete 282-sample dataset and not yet QC-filtered; contains no indel calls",
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("SUCCESS: one genotype chunk read, verified, and cached.", flush=True)
    print(f"Pilot sample genotypes shape: {subset.shape}", flush=True)
    print(f"Elapsed compute time: {elapsed:.1f} seconds", flush=True)
    print(f"Cache file: {cache_file} ({cache_file.stat().st_size:,} bytes)", flush=True)
    print(f"Report: {REPORT}", flush=True)
    print("NOTE: 12 other genotype chunks were NOT requested.", flush=True)


if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        print(f"PILOT ERROR: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
        raise
