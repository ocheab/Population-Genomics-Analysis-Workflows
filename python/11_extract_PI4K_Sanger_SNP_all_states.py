#!/usr/bin/env python3
"""Pf8 public Sanger SNP-only, bounded, resumable PI4K sample-chunk extractor.

Requires malariagen-data, numpy, bcftools and bgzip installed in active environment.
Outputs a derived, genotype-only VCF (not the original Pf8 VCF), plus NPZ matrices.
"""
import csv
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import malariagen_data
import aiohttp
from fsspec.exceptions import FSTimeoutError

STATES = {"Oyo": 71, "Ogun": 38, "Lagos": 129, "Osun": 15, "Ondo": 29}
REGION = "Pf3D7_05_v3:409403-415758"
CHROM, COORDS = REGION.split(":")
START, END = map(int, COORDS.split("-"))
URL = "https://pf8-release.cog.sanger.ac.uk/snp-only"
EXPECTED_SITES, EXPECTED_PASS = 2471, 1824
ROOT = Path(__file__).resolve().parents[1]
META = ROOT / 'results' / 'Southwest5'
OUT = ROOT / 'variants' / 'Southwest5' / 'Sanger_SNP_complete'
CACHE = OUT / 'verified_chunks'
LOG = OUT / 'progress.jsonl'
MAX_FETCH_ATTEMPTS = 6
FETCH_RETRY_DELAYS_SEC = (5, 10, 20, 30, 45)


def note(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def append_log(event, **kwargs):
    with LOG.open('a', encoding='utf-8') as f:
        f.write(json.dumps({'timestamp': datetime.now(timezone.utc).isoformat(), 'event': event, **kwargs}, sort_keys=True) + '\n')


def fetch_native_genotype_chunk(gt, variant_chunk, sample_chunk, expected_shape):
    """Retry transient remote I/O errors without altering or invalidating saved chunks.

    Each attempt reads exactly the same source Zarr chunk. The caller writes the
    selected sample subset to its verified cache only after a complete successful
    read and shape validation; failed reads are never cached.
    """
    transient_errno = {
        errno.ECONNRESET, errno.ETIMEDOUT, errno.ECONNABORTED,
        errno.EPIPE, errno.ENETUNREACH, errno.EHOSTUNREACH,
    }
    for attempt in range(1, MAX_FETCH_ATTEMPTS + 1):
        try:
            chunk = np.asarray(
                gt.blocks[variant_chunk, sample_chunk, 0].compute(),
                dtype=np.int8,
            )
            if chunk.shape != expected_shape:
                raise RuntimeError(
                    f'Unexpected native chunk shape {chunk.shape}; '
                    f'expected {expected_shape}. Refusing to cache it.'
                )
            if attempt > 1:
                note(f'Source chunk {sample_chunk}: remote read succeeded '
                     f'on attempt {attempt}/{MAX_FETCH_ATTEMPTS}.')
            return chunk
        except (FSTimeoutError, TimeoutError, ConnectionError,
                aiohttp.ClientError, OSError) as exc:
            # Retry network interruptions, NOT permission/404 errors or
            # unexpected server responses that require intervention.
            if isinstance(exc, OSError) and not isinstance(
                exc, (FSTimeoutError, TimeoutError, ConnectionError,
                      aiohttp.ClientError)
            ):
                if exc.errno not in transient_errno:
                    raise
            if isinstance(exc, aiohttp.ClientResponseError) and exc.status not in (408, 429, 500, 502, 503, 504):
                raise
            kind = type(exc).__name__
            if attempt == MAX_FETCH_ATTEMPTS:
                append_log('chunk_retry_exhausted', sample_chunk=sample_chunk,
                           attempts=attempt, error_type=kind, error=str(exc)[:250])
                note(f'Source chunk {sample_chunk}: retry limit reached. '
                     'Previously verified chunk caches remain intact.')
                raise
            delay = FETCH_RETRY_DELAYS_SEC[attempt - 1]
            append_log('chunk_retry', sample_chunk=sample_chunk,
                       attempt=attempt, error_type=kind, backoff_sec=delay)
            note(f'Source chunk {sample_chunk}: {kind} during network read '
                 f'(attempt {attempt}/{MAX_FETCH_ATTEMPTS}); '
                 f'retrying in {delay}s. No partial genotypes saved.')
            time.sleep(delay)
    raise AssertionError('Unreachable retry loop exit')


def read_ids(path):
    ids = [s.strip() for s in path.read_text(encoding='utf-8-sig').splitlines() if s.strip()]
    if len(ids) != len(set(ids)):
        raise RuntimeError(f"Duplicate sample IDs in {path}")
    return ids


def site_sha(positions):
    return hashlib.sha256(np.asarray(positions, dtype='<i4').tobytes()).hexdigest()


def write_npz_atomic(path, **values):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'wb') as f:
        np.savez_compressed(f, **values)
    os.replace(tmp, path)


def load_verified_chunk(path, ids, positions, sc):
    if not path.is_file():
        return None
    try:
        with np.load(path, allow_pickle=False) as z:
            got_ids = z['sample_ids'].astype(str).tolist()
            got_pos = z['positions'].astype('i4')
            g = z['genotypes'].astype('i1')
            nsc = int(z['source_sample_chunk'])
        if (got_ids == ids and nsc == sc and np.array_equal(got_pos, positions)
                and g.shape == (len(positions), len(ids), 2)):
            return g
        note(f'Cache mismatch at sample chunk {sc}; will refetch.')
    except Exception as e:
        note(f'Cache cannot be verified at sample chunk {sc} ({type(e).__name__}); will refetch.')
    return None


def pick_string(z, names):
    for name in names:
        if name in z:
            arr = np.asarray(z[name])
            if arr.ndim == 1 and arr.dtype.kind in ('U','S'):
                return arr.astype(str).tolist()
    return None


def pilot_cache_try(pilot_file, expected_ids, expected_pos):
    """Reuse prior Step 15 cache ONLY if it independently provides verifiable IDs and coordinates."""
    if not pilot_file.exists():
        return None
    try:
        with np.load(pilot_file, allow_pickle=False) as z:
            keys = z.files
            ids = pick_string(z, ['sample_ids','samples','sample_id','selected_sample_ids'])
            pos = None
            for key in ['positions','variant_positions','variant_position','site_positions']:
                if key in z:
                    pos = np.asarray(z[key], dtype='i4')
                    break
            g = None
            for key in ['genotypes','call_genotype','gt','GT']:
                if key in z:
                    g = np.asarray(z[key], dtype='i1')
                    break
            if g is None:
                arrays = [np.asarray(z[k]) for k in keys if np.asarray(z[k]).ndim == 3]
                if len(arrays) == 1:
                    g = arrays[0].astype('i1')
        if ids == expected_ids and pos is not None and np.array_equal(pos, expected_pos) and g is not None and g.shape == (len(pos),len(ids),2):
            note('Verified and reusing Step 15 pilot cache (same positions and sample IDs).')
            return g
        note('Step 15 cache lacks fully verifiable IDs and site coordinates; safely re-reading one source chunk.')
    except Exception as e:
        note(f'Step 15 cache not reusable ({type(e).__name__}: {e}); fetching source chunk.')
    return None


def allele_str(v):
    if v is None:
        return ''
    if isinstance(v, (bytes,np.bytes_)):
        return bytes(v).decode('utf-8').strip()
    text = str(v).strip()
    if text.lower() in ('none','nan','null') or text == '.':
        return ''
    return text


def make_vcf(positions, alleles, passes, ids, genotypes, output):
    """Derived VCF with only basic REF/ALT, boolean Pf8 filter and GT; not source VCF."""
    if shutil.which('bgzip') is None or shutil.which('bcftools') is None:
        raise RuntimeError('bgzip or bcftools is not on PATH; NPZ data remain valid.')
    tmp_tsv = output.with_name(output.name + '.building.vcf')
    tmp_gz = output.with_name(output.name + '.building.gz')
    try:
        with tmp_tsv.open('w', encoding='utf-8', newline='\n') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##source=MalariaGEN_Pf8_Sanger_public_snp_only_derived\n')
            f.write('##reference=PlasmoDB-54-Pfalciparum3D7-Genome\n')
            f.write(f'##contig=<ID={CHROM}>\n')
            f.write('##FILTER=<ID=PF8_NONPASS,Description="Pf8 variant_filter_pass is false; original filter codes not retained">\n')
            f.write('##INFO=<ID=PF8_FILTER_PASS,Number=1,Type=Integer,Description="Pf8 Zarr variant_filter_pass Boolean">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Pf8 call_genotype two-slot representation, not a biological ploidy inference">\n')
            f.write('##pi4k_note=Derived minimal VCF; genotype quality, depth, variant QUAL and original filter labels not included; SNPs only\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n')
            for vi, p in enumerate(positions):
                row = [allele_str(x) for x in alleles[vi]]
                if not row or not row[0]:
                    raise ValueError(f'Missing REF allele at position {p}; refusing potentially incorrect VCF')
                while len(row)>1 and not row[-1]:
                    row.pop()
                if not all(row):
                    raise ValueError(f'Internal empty ALT allele at {p}; refusing to remap genotypes')
                alt = ','.join(row[1:]) if len(row)>1 else '.'
                g = genotypes[vi]
                if g.size and (g.max(initial=0)>=len(row) or g.min(initial=0)<-1):
                    raise ValueError(f'Unrepresentable allele index at {p}: {g.min()}, {g.max()}, n_alleles={len(row)}')
                calls = [('/'.join('.' if int(x) < 0 else str(int(x)) for x in xrow)) for xrow in g]
                filt = 'PASS' if bool(passes[vi]) else 'PF8_NONPASS'
                f.write(f'{CHROM}\t{int(p)}\t.\t{row[0]}\t{alt}\t.\t{filt}\tPF8_FILTER_PASS={int(bool(passes[vi]))}\tGT\t'+'\t'.join(calls)+'\n')
        with tmp_gz.open('wb') as dest:
            subprocess.run(['bgzip','-c',str(tmp_tsv)], stdout=dest, check=True)
        subprocess.run(['bcftools','index','--csi','--force',str(tmp_gz)], check=True)
        os.replace(tmp_gz, output)
        os.replace(str(tmp_gz)+'.csi',str(output)+'.csi')
        subprocess.run(['bcftools','quickcheck','-v',str(output)], check=True)
        actual_ids = subprocess.check_output(['bcftools','query','-l',str(output)], text=True).splitlines()
        if actual_ids != ids:
            raise RuntimeError(f'VCF sample ID mismatch in {output}')
        note(f'Created derived VCF: {output.name} ({output.stat().st_size:,} bytes).')
    finally:
        if tmp_tsv.exists(): tmp_tsv.unlink()
        if tmp_gz.exists(): tmp_gz.unlink()
        tmp_index = Path(str(tmp_gz)+'.csi')
        if tmp_index.exists(): tmp_index.unlink()


def run():
    CACHE.mkdir(parents=True, exist_ok=True)
    note('=== STEP 16: PI4K Pf8 public Sanger SNP extraction ===')
    note('This script reads only selected genomic chunks and 282 samples; source has NO indels.')
    note(f'Network read retry policy: up to {MAX_FETCH_ATTEMPTS} attempts per uncached chunk.')
    expected_ids_by_state = {}
    for state, expected_n in STATES.items():
        path = META / 'clean_sample_lists' / f'{state}_sample_IDs.txt'
        if not path.exists():
            path = META / 'sample_lists' / f'{state}_sample_IDs.txt'
        if not path.exists():
            raise FileNotFoundError(f'Missing state sample ID list: {state}')
        ids = read_ids(path)
        if len(ids) != expected_n: raise RuntimeError(f'{state} expected {expected_n} IDs but found {len(ids)}')
        expected_ids_by_state[state] = ids
        note(f'  {state}: {len(ids)} sample IDs')
    all_ids = [sid for ids in expected_ids_by_state.values() for sid in ids]
    if len(all_ids) != 282 or len(set(all_ids)) != 282:
        raise RuntimeError('Expected exactly 282 distinct sample IDs across five states')
    p = malariagen_data.Pf8(URL)
    ds = p.variant_calls(chunks='native')
    src_ids = [str(x) for x in np.asarray(ds['sample_id'].values)]
    src_id_lookup = {sid: i for i, sid in enumerate(src_ids)}
    missing = [x for x in all_ids if x not in src_id_lookup]
    if missing: raise RuntimeError(f'{len(missing)} sample IDs missing from Pf8 genotype callset: {missing[:5]}')

    gt = ds['call_genotype'].data
    chrom_a, pos_a, pass_a, alleles_a = [ds[k].data for k in ['variant_chrom','variant_position','variant_filter_pass','variant_allele']]
    if gt.ndim != 3 or gt.shape[2] != 2 or len(gt.chunks[2]) != 1:
        raise RuntimeError(f'Unexpected genotype dimensions: shape={gt.shape} chunks={gt.chunks}')
    if len(gt.chunks[0]) != len(chrom_a.chunks[0]) or len(gt.chunks[0]) != len(pos_a.chunks[0]):
        raise RuntimeError('Unexpected variant chunk alignment')
    vchunks = []
    variant_offsets = np.cumsum([0]+list(gt.chunks[0]))
    seen_chrom = False
    note('Locating PI4K variant chunk using small position/chromosome metadata reads...')
    for iv in range(len(gt.chunks[0])):
        chrs = np.asarray(chrom_a.blocks[iv].compute()).astype(str)
        target_chr_mask = chrs == CHROM
        if not np.any(target_chr_mask):
            if seen_chrom: break
            continue
        seen_chrom = True
        ps = np.asarray(pos_a.blocks[iv].compute(), dtype=np.int32)
        loc = np.flatnonzero(target_chr_mask & (ps >= START) & (ps <= END))
        if len(loc):
            fp = np.asarray(pass_a.blocks[iv].compute(), dtype=bool)[loc]
            # allele metadata has one chunk per allele; slice region after lazy selection
            aval = np.asarray(alleles_a[variant_offsets[iv] + loc, :].compute())
            vchunks.append((iv,loc,ps[loc],fp,aval))
            note(f'  Found {len(loc)} sites in variant chunk {iv}.')
        if np.any(target_chr_mask & (ps > END)):
            break
    if not vchunks: raise RuntimeError('PI4K locus could not be found in public SNP callset')
    if len(vchunks) != 1:
        raise RuntimeError(f'Expected one native variant chunk from planner, found {len(vchunks)}; stop before extraction')
    vi, vlocal, positions, site_pass, alleles = vchunks[0]
    if len(positions) != EXPECTED_SITES or site_pass.sum() != EXPECTED_PASS:
        raise RuntimeError(f'Unexpected PI4K positions or PASS count: {len(positions)}/{site_pass.sum()}')
    if not np.all(np.diff(positions) > 0):
        raise RuntimeError('Variant positions not strictly increasing; stop rather than write invalid VCF')
    note(f'PI4K site inventory verified: {len(positions)} SNPs ({int(site_pass.sum())} PASS).')

    source_idx = np.asarray([src_id_lookup[x] for x in all_ids], dtype=np.int64)
    sample_offsets = np.cumsum([0]+list(gt.chunks[1]))
    selected_chunks = sorted(set(np.searchsorted(sample_offsets[1:], source_idx, side='right').tolist()))
    if len(selected_chunks) != 13:
        raise RuntimeError(f'Expected 13 sample chunks from planner; found {len(selected_chunks)}')
    result = np.full((len(positions), 282, 2), -1, dtype=np.int8)
    pilot_path = ROOT/'variants'/'Southwest5'/'Sanger_SNP_chunk_cache'/'PI4K_native_sample_chunk_238.npz'
    report_path = META/'PI4K_Sanger_SNP_pilot_report.json'
    if report_path.exists():
        rep = json.loads(report_path.read_text(encoding='utf-8'))
        if rep.get('status') != 'ONE_CHUNK_PILOT_SUCCESS' or rep.get('region') != REGION:
            note('Pilot report cannot be validated against region; will not rely on Step 15 cache')
            pilot_path = Path('/nonexistent')
    for n, sc in enumerate(selected_chunks, 1):
        sample_mask = (source_idx >= sample_offsets[sc]) & (source_idx < sample_offsets[sc+1])
        local_dest_idx = np.flatnonzero(sample_mask)
        relative_src_idx = (source_idx[local_dest_idx] - sample_offsets[sc]).astype(int)
        ids = [all_ids[i] for i in local_dest_idx]
        cache_path = CACHE/f'chunk_{sc:03d}.npz'
        g = load_verified_chunk(cache_path, ids, positions, sc)
        if g is not None:
            note(f'Chunk {n}/13 (source {sc}): verified local cache ({len(ids)} samples).')
        elif sc == 238:
            g = pilot_cache_try(pilot_path, ids, positions)
            if g is not None:
                append_log('reused_pilot', sample_chunk=sc, n_samples=len(ids))
        if g is None:
            note(f'Chunk {n}/13 (source {sc}): reading one native block; {len(ids)} selected samples...')
            t = time.monotonic()
            whole_block = fetch_native_genotype_chunk(
                gt, vi, sc, expected_shape=(gt.chunks[0][vi], gt.chunks[1][sc], 2)
            )
            g = np.array(whole_block[vlocal[:, None], relative_src_idx[None, :], :], copy=True)
            del whole_block
            if g.shape != (len(positions),len(ids),2):
                raise RuntimeError(f'Unexpected selected chunk shape: {g.shape}')
            append_log('fetched', sample_chunk=sc, n_samples=len(ids), elapsed_sec=round(time.monotonic()-t,2))
            note(f'Chunk {n}/13: retrieved and selected {g.shape} genotypes.')
        result[:,local_dest_idx,:] = g
        write_npz_atomic(cache_path, sample_ids=np.asarray(ids,dtype='U'), positions=positions,
                         genotypes=g, source_sample_chunk=np.int32(sc))
    note('All 282 selected samples were assembled in original state-list order.')

    idsarr = np.asarray(all_ids,dtype='U')
    states = np.asarray([s for s,v in expected_ids_by_state.items() for _ in v],dtype='U')
    all_npz = OUT/'PI4K_Southwest5_282_Sanger_SNP_all_sites.npz'
    write_npz_atomic(all_npz, positions=positions, sample_ids=idsarr, states=states,
                     genotypes=result, filter_pass=site_pass, alleles=alleles.astype('U'),
                     source=np.asarray(URL), region=np.asarray(REGION))
    pass_npz = OUT/'PI4K_Southwest5_282_Sanger_SNP_PASS_sites.npz'
    write_npz_atomic(pass_npz, positions=positions[site_pass], sample_ids=idsarr,
                     states=states, genotypes=result[site_pass], alleles=alleles[site_pass].astype('U'),
                     source=np.asarray(URL), region=np.asarray(REGION))
    state_rows=[]
    outvcfs=[]
    for state, ids in expected_ids_by_state.items():
        subset_ix = np.flatnonzero(states == state)
        one = result[:,subset_ix,:]
        write_npz_atomic(OUT/f'PI4K_{state}_Sanger_SNP_all_sites.npz',
                         positions=positions, sample_ids=np.asarray(ids,dtype='U'), genotypes=one,
                         filter_pass=site_pass, alleles=alleles.astype('U'))
        missing_alleles = int((one<0).sum())
        total_alleles = int(one.size)
        state_rows.append(dict(State=state, Samples=len(ids), SNP_sites_all=len(positions),
                               SNP_sites_Pf8_PASS=int(site_pass.sum()), Missing_allele_proportion=f'{missing_alleles/total_alleles:.6f}',
                               Fully_called_sample_site_proportion=f'{np.all(one>=0,axis=2).mean():.6f}'))
    with (OUT/'PI4K_state_data_QC_preliminary.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=state_rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(state_rows)
    note('Saved 282-sample full and PASS-only NPZ matrices, five state NPZ files, and basic missingness report.')
    try:
        vcf_all = OUT/'PI4K_Southwest5_282_Sanger_SNP_derived.vcf.gz'
        make_vcf(positions, alleles, site_pass, all_ids, result, vcf_all)
        outvcfs.append(str(vcf_all))
        for state, ids in expected_ids_by_state.items():
            idx=np.flatnonzero(states==state)
            path=OUT/f'PI4K_{state}_Sanger_SNP_derived.vcf.gz'
            make_vcf(positions, alleles, site_pass, ids, result[:,idx,:], path)
            outvcfs.append(str(path))
    except Exception as e:
        note(f'VCF export blocked safely: {type(e).__name__}: {e}. NPZ matrices are already saved; no fake VCF created.')
        append_log('vcf_error', error=str(e))
    summary={
        'status':'ALL_282_SNP_GENOTYPES_EXTRACTED',
        'source':URL,'gene':'PF3D7_0509800','region':REGION,
        'sample_total':len(all_ids),'states':STATES,'snp_sites_all':len(positions),
        'snp_sites_pf8_pass':int(site_pass.sum()),'native_variant_chunks':len(vchunks),
        'native_sample_chunks':len(selected_chunks),'native_genotype_chunks_read_or_cached':len(selected_chunks),
        'native_chunk_uncompressed_footprint_mib':round(sum(gt.chunks[0][vi]*gt.chunks[1][sc]*gt.dtype.itemsize*2 for sc in selected_chunks)/2**20,2),
        'actual_network_bytes':'not measured','site_position_sha256':site_sha(positions),
        'all_sites_npz':str(all_npz),'pass_sites_npz':str(pass_npz),'derived_vcfs':outvcfs,
        'important_limitations': ['SNP-only (no indels)','Pf8 site-level PASS is not per-sample genotype QC','two-slot genotype codes do not imply biological diploidy','derived VCF lacks original per-genotype quality and depth fields']
    }
    (OUT/'PI4K_Sanger_SNP_full_extraction_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    append_log('success', samples=len(all_ids), sites=len(positions), pass_sites=int(site_pass.sum()), n_vcf=len(outvcfs))
    note('SUCCESS: 282 samples and all PI4K SNP sites extracted. See summary:')
    note(str(OUT/'PI4K_Sanger_SNP_full_extraction_summary.json'))
    if not outvcfs: note('NOTE: No VCF files exported; check progress.jsonl for VCF error and retain NPZ files.')

if __name__=='__main__':
    try:
        run()
    except Exception as e:
        print(f'ERROR: {type(e).__name__}: {e}', file=sys.stderr, flush=True)
        raise
