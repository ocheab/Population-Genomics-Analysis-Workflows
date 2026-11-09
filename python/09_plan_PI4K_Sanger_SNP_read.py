#!/usr/bin/env python3
"""Metadata-only planning of a bounded Pf8 PI4K SNP genotype retrieval.

This script NEVER accesses or computes call_genotype, call_AD, or any other
sample-level genotype matrix. It fetches variant coordinates and lightweight
variant filters to identify the PI4K records and estimates native Zarr chunks.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import sys
import numpy as np
import malariagen_data

ROOT = Path('/mnt/c/PI4K_Nigeria_Study')
BASE = 'https://pf8-release.cog.sanger.ac.uk/snp-only'
CHROM = 'Pf3D7_05_v3'
START, END = 409403, 415758  # 1-based, inclusive
STATES = {'Oyo': 71, 'Ogun': 38, 'Lagos': 129, 'Osun': 15, 'Ondo': 29}
SAMPLE_DIR = ROOT / 'results/Southwest5/sample_lists'
CLEAN = ROOT / 'results/Southwest5/clean_sample_lists/Southwest5_282_sample_IDs.txt'
OUTPUT = ROOT / 'results/Southwest5'


def normalize(value):
    if isinstance(value, bytes):
        return value.decode('utf-8').rstrip('\x00')
    return str(value).rstrip('\x00')


def ids_from_file(path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing sample list: {path}')
    xs = [s.strip() for s in path.read_text(encoding='utf-8-sig').splitlines() if s.strip()]
    if len(xs) != len(set(xs)):
        raise ValueError(f'Duplicate sample IDs in {path}')
    return xs


def chunk_numbers(indices, sizes):
    ends = np.cumsum(np.asarray(sizes, dtype=np.int64))
    return np.unique(np.searchsorted(ends, indices, side='right'))


def run():
    print('=== Pf8 PI4K SNP retrieval plan (NO GENOTYPE DOWNLOAD) ===', flush=True)
    if not ROOT.is_dir():
        raise FileNotFoundError('Study directory not found: ' + str(ROOT))
    all_ids = []
    for state, expected in STATES.items():
        ids = ids_from_file(SAMPLE_DIR / f'{state}_sample_IDs.txt')
        if len(ids) != expected:
            raise ValueError(f'{state}: expected {expected} IDs, got {len(ids)}')
        all_ids.extend(ids)
        print(f'OK {state}: {len(ids)} IDs', flush=True)
    if len(all_ids) != 282 or len(set(all_ids)) != 282:
        raise ValueError('Expected exactly 282 unique IDs in five state lists')
    if CLEAN.exists() and set(ids_from_file(CLEAN)) != set(all_ids):
        raise ValueError('Five-state union differs from the previously cleaned combined ID file')

    pf8 = malariagen_data.Pf8(BASE)
    ds = pf8.variant_calls()
    nvar = int(ds.sizes['variants'])
    if 'variant_chrom' not in ds or 'variant_position' not in ds or 'sample_id' not in ds:
        raise ValueError('Missing expected Pf8 variant coordinates or sample_id')
    genotype = ds['call_genotype'].data  # Metadata reference only: NO .compute()!
    if genotype.ndim != 3:
        raise ValueError(f'Unexpected genotype dimensions: {genotype.shape}')
    if genotype.chunks is None:
        raise ValueError('Call genotype is not Dask-chunked as expected')
    chunks_v, chunks_s, chunks_p = genotype.chunks

    sample_names = [normalize(x) for x in ds['sample_id'].data.compute()]
    if len(sample_names) != len(set(sample_names)):
        raise ValueError('Source sample identifiers are not unique')
    lookup = {name: i for i, name in enumerate(sample_names)}
    not_found = sorted(set(all_ids) - lookup.keys())
    if not_found:
        raise ValueError(f'{len(not_found)} local IDs missing from SNP dataset; e.g. {not_found[:5]}')
    sample_idx = np.array(sorted(lookup[x] for x in all_ids), dtype=np.int64)
    print('Matched all 282 samples in SNP genotype array.', flush=True)

    # Small integer-stride probe. This samples one chromosome label per native
    # variant chunk, rather than materialising 10.8 million chromosome strings.
    first_chunk_len = int(chunks_v[0])
    print(f'Scanning chromosome labels at {first_chunk_len:,}-variant intervals...', flush=True)
    probes = [normalize(x) for x in ds['variant_chrom'].data[::first_chunk_len].compute()]
    matches = [i for i, label in enumerate(probes) if label == CHROM]
    if not matches:
        raise RuntimeError('Chromosome not found at sampled chunk starts; do not assume positions')
    if matches != list(range(matches[0], matches[-1] + 1)):
        raise RuntimeError('Chromosome appears in discontinuous blocks; unsafe to infer interval')

    # Cover one preceding and one following probe interval for chromosome edges.
    lo = max(0, (matches[0] - 1) * first_chunk_len)
    hi = min(nvar, (matches[-1] + 2) * first_chunk_len)
    print(f'Reading variant coordinate metadata in window [{lo:,}, {hi:,})...', flush=True)
    chr_window = ds['variant_chrom'].data[lo:hi].compute()
    pos_window = np.asarray(ds['variant_position'].data[lo:hi].compute(), dtype=np.int64)
    in_gene = np.fromiter(
        (normalize(c) == CHROM and START <= int(p) <= END
         for c, p in zip(chr_window, pos_window)),
        dtype=np.bool_, count=len(pos_window)
    )
    relative = np.flatnonzero(in_gene)
    if not len(relative):
        raise RuntimeError('No PI4K SNP sites found. Verify chromosome and coordinates.')
    genomic_idx = relative + lo
    positions = pos_window[relative]
    if not np.all(positions[:-1] <= positions[1:]):
        raise ValueError('PI4K variant positions are not ordered')
    if np.any(positions < START) or np.any(positions > END):
        raise AssertionError('PI4K positions extend outside gene bounds')
    print(f'PI4K SNP positions found: {len(positions):,}', flush=True)

    var_chunk_idx = chunk_numbers(genomic_idx, chunks_v)
    sample_chunk_idx = chunk_numbers(sample_idx, chunks_s)
    native_raw_bytes = (sum(chunks_v[i] for i in var_chunk_idx)
                        * sum(chunks_s[i] for i in sample_chunk_idx)
                        * sum(chunks_p) * np.dtype(genotype.dtype).itemsize)
    n_chunks = len(var_chunk_idx) * len(sample_chunk_idx) * len(chunks_p)
    print(f'Native genotype chunks intersected: {n_chunks:,} '
          f'({len(var_chunk_idx)} variant x {len(sample_chunk_idx)} sample x {len(chunks_p)} ploidy)', flush=True)
    print(f'Combined uncompressed size of those genotype chunks: '
          f'{native_raw_bytes / 2**30:.2f} GiB', flush=True)
    print('IMPORTANT: This is NOT an estimate of internet bytes; '
          'remote chunks are compressed and transfer size is unknown.', flush=True)

    if 'variant_filter_pass' in ds:
        pass_flags = np.asarray(ds['variant_filter_pass'].data[genomic_idx].compute(), dtype=bool)
    else:
        pass_flags = np.full(len(genomic_idx), False, dtype=bool)
        print('WARNING: variant_filter_pass not available; flags set to unknown in export')

    OUTPUT.mkdir(parents=True, exist_ok=True)
    site_path = OUTPUT / 'PI4K_Sanger_SNP_site_inventory.tsv'
    with site_path.open('w', encoding='utf-8', newline='') as f:
        wr = csv.writer(f, delimiter='\t')
        wr.writerow(['chromosome', 'position_1based', 'global_variant_index_0based', 'filter_pass'])
        for p, idx, passed in zip(positions, genomic_idx, pass_flags):
            wr.writerow([CHROM, int(p), int(idx), bool(passed)])

    data = {
        'data_source': BASE, 'gene': 'PF3D7_0509800',
        'region_1based_inclusive': f'{CHROM}:{START}-{END}',
        'genotype_matrix_downloaded': False,
        'pf8_snp_positions_genomewide': nvar,
        'selected_sample_ids': len(sample_idx),
        'selected_snp_positions': len(genomic_idx),
        'selected_filter_pass_sites': int(pass_flags.sum()),
        'native_variant_chunks_needed': len(var_chunk_idx),
        'native_sample_chunks_needed': len(sample_chunk_idx),
        'native_genotype_chunks_intersected': int(n_chunks),
        'uncompressed_native_genotype_chunk_footprint_bytes': int(native_raw_bytes),
        'uncompressed_native_genotype_chunk_footprint_GiB': round(native_raw_bytes / 2**30, 4),
        'network_transfer_bytes': 'unknown; Zarr chunks compressed, no genotype read attempted',
        'state_QC_pass_count': STATES,
    }
    summary = OUTPUT / 'PI4K_Sanger_SNP_retrieval_plan.json'
    summary.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    print('Saved:', summary, flush=True)
    print('Saved:', site_path, flush=True)
    print('SUCCESS: coordinate-only planning completed; genotype reads: NONE.', flush=True)


if __name__ == '__main__':
    try:
        run()
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}', file=sys.stderr, flush=True)
        raise
