#!/usr/bin/env python3
"""PI4K Southwest-five Step 21: offline same-study-year 2020 cohort.

Uses the already-audited Step 20 membership manifest and Step 17 local VCFs.
NO remote requests. NO edits to source datasets. Does not infer genotype ploidy.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import os

import pandas as pd

STATES = ('Oyo', 'Ogun', 'Lagos', 'Osun', 'Ondo')
FULL_COUNTS = {'Oyo': 71, 'Ogun': 38, 'Lagos': 129, 'Osun': 15, 'Ondo': 29}
YEAR_COUNTS = {'Oyo': 16, 'Ogun': 38, 'Lagos': 21, 'Osun': 15, 'Ondo': 29}
STUDY = '1306-PF-NG-NGWA-SM'
YEAR = '2020'


def die(msg: str) -> None:
    raise ValueError(msg)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def run_cmd(args: list[str]) -> str:
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(args)}\n{proc.stderr[-1800:]}")
    return proc.stdout.strip()


def ids_from_vcf(vcf: Path) -> list[str]:
    return [x.strip() for x in run_cmd(['bcftools', 'query', '-l', str(vcf)]).splitlines() if x.strip()]


def count_from_index(vcf: Path) -> int:
    return int(run_cmd(['bcftools', 'index', '-n', str(vcf)]))


def subset_vcf(source: Path, ids: list[str], id_path: Path, output: Path, expect_sites: int) -> dict:
    if not source.is_file():
        raise FileNotFoundError(f'VCF source missing: {source}')
    if not ids:
        die(f'Refusing empty-sample subsetting for {output}')
    if output.exists() or Path(str(output)+'.csi').exists():
        # Verify before replacing stale or partial results.
        if output.exists() and Path(str(output)+'.csi').is_file():
            try:
                if (set(ids_from_vcf(output)) == set(ids) and
                        len(ids_from_vcf(output)) == len(ids) and
                        count_from_index(output) == expect_sites):
                    return {'file':str(output), 'samples':len(ids), 'sites':expect_sites,
                            'bytes':output.stat().st_size, 'sha256':sha256(output), 'action':'verified_existing'}
            except (RuntimeError, ValueError):
                pass
        die(f'Existing output is incomplete or mismatched: {output}. Move it aside manually before rerun.')
    tmp = output.with_name(output.name.replace('.vcf.gz', '.partial.vcf.gz'))
    tmpidx = Path(str(tmp)+'.csi')
    if tmp.exists() or tmpidx.exists():
        die(f'Incomplete temporary output from previous attempt: {tmp}. Review before removing.')
    try:
        run_cmd(['bcftools','view','--samples-file', str(id_path), '-Oz', '-o', str(tmp), str(source)])
        run_cmd(['bcftools','index','-f','-c',str(tmp)])
        found = ids_from_vcf(tmp)
        if len(found) != len(ids) or set(found) != set(ids):
            die(f'Sample ID mismatch for {tmp}: expected {len(ids)}, got {len(found)}')
        sites = count_from_index(tmp)
        if sites != expect_sites:
            die(f'Unexpected site count for {tmp}: {sites} vs {expect_sites}')
        os.replace(tmp, output)
        os.replace(tmpidx, Path(str(output)+'.csi'))
    except Exception:
        # Leave temporary file for troubleshooting; original intact.
        raise
    return {'file':str(output), 'samples':len(ids), 'sites':sites,
            'bytes':output.stat().st_size, 'sha256':sha256(output), 'action':'created'}


def main(root: Path, manifest_only: bool) -> None:
    manifest = root/'results/Southwest5/Step20_Provenance_Cohorts/PI4K_cohort_membership_audit_282.tsv'
    if not manifest.is_file():
        raise FileNotFoundError(f'Missing Step 20 file: {manifest}')
    df = pd.read_csv(manifest, sep='\t', dtype=str, encoding='utf-8-sig', keep_default_na=False)
    required = {'sample_id','state','study','year','Fws','pass_missing_slot_fraction'}
    if not required.issubset(df.columns):
        die(f'Missing columns in Step 20 manifest: {sorted(required-set(df.columns))}')
    df = df.copy()
    for c in ('sample_id','state','study','year'):
        df[c] = df[c].astype(str).str.strip()
    if len(df)!=282 or df.sample_id.nunique()!=282 or df.sample_id.eq('').any():
        die('Step 20 file no longer has exactly 282 distinct nonblank sample IDs')
    if df.state.value_counts().to_dict() != FULL_COUNTS:
        die(f'Unexpected full state counts: {df.state.value_counts().to_dict()}')
    df['Fws'] = pd.to_numeric(df['Fws'], errors='raise')
    df['pass_missing_slot_fraction'] = pd.to_numeric(df['pass_missing_slot_fraction'],errors='raise')
    if not df['Fws'].between(0,1).all() or not df['pass_missing_slot_fraction'].between(0,1).all():
        die('Unexpected Fws or missingness value outside [0,1]')
    all_ids = set(df.sample_id)

    df2020 = df.loc[(df.study == STUDY) & (df.year == YEAR)].copy()
    observed = df2020.state.value_counts().to_dict()
    if observed != YEAR_COUNTS or len(df2020) != 119:
        die(f'Expected 119 samples in same study/year, found {len(df2020)} and counts {observed}')
    if df2020.sample_id.nunique()!=119:
        die('Repeated 2020 sample IDs found')

    out = root/'results/Southwest5/Step21_2020_Controlled_Cohort'
    vcf_out = root/'variants/Southwest5/Sanger_SNP_complete/Step21_2020'
    out.mkdir(parents=True, exist_ok=True)
    all_sorted = df2020.sort_values(['state','sample_id'])
    ids_combined = sorted(df2020.sample_id.tolist())
    list_paths: dict[str,Path] = {}
    combo = out/'PI4K_2020_same_study_119_sample_IDs.txt'
    combo.write_text('\n'.join(ids_combined)+'\n',encoding='utf-8')
    list_paths['COMBINED'] = combo
    for state in STATES:
        ids = sorted(df2020.loc[df2020.state==state, 'sample_id'].tolist())
        path = out/f'PI4K_2020_{state}_sample_IDs.txt'
        path.write_text('\n'.join(ids)+'\n',encoding='utf-8')
        list_paths[state] = path
    all_sorted.to_csv(out/'PI4K_2020_same_study_sample_manifest.tsv',sep='\t',index=False,float_format='%.6f')

    rows = []
    for state in (*STATES, 'TOTAL'):
        s = df2020 if state=='TOTAL' else df2020.loc[df2020.state==state]
        rows.append({'state':state, 'n_same_study_2020':len(s),
                     'n_pass_missing_le10':int((s.pass_missing_slot_fraction<=.10).sum()),
                     'n_fws_ge095':int((s.Fws>=.95).sum()),
                     'n_joint_fws095_pass_missing_le10':int(((s.Fws>=.95)&(s.pass_missing_slot_fraction<=.10)).sum()),
                     'median_Fws':float(s.Fws.median()),
                     'median_PASS_site_missingness':float(s.pass_missing_slot_fraction.median())})
    summary_tbl = pd.DataFrame(rows)
    summary_tbl.to_csv(out/'PI4K_2020_cohort_counts_and_QC.tsv',sep='\t',index=False,float_format='%.6f')
    # Informative, but not used as automatic exclusion.
    for nm, flag in [('pass_missing_le10',df2020.pass_missing_slot_fraction<=.10),
                     ('fws_ge095',df2020.Fws>=.95),
                     ('joint_fws095_missing_le10',(df2020.Fws>=.95)&(df2020.pass_missing_slot_fraction<=.10))]:
        ids = sorted(df2020.loc[flag,'sample_id'].tolist())
        (out/f'PI4K_2020_exploratory_{nm}_sample_IDs.txt').write_text('\n'.join(ids)+'\n',encoding='utf-8')

    print('STEP 21: SAME STUDY/YEAR 2020 COHORT -- OFFLINE')
    print('Source cohort: 282 samples; controlled study-year cohort: 119 samples')
    print(summary_tbl.to_string(index=False))
    source_vcf_dir = root/'variants/Southwest5/Sanger_SNP_complete'
    allvcf = source_vcf_dir/'PI4K_Southwest5_282_Sanger_SNP_derived.vcf.gz'
    passvcf = source_vcf_dir/'PI4K_Southwest5_282_Sanger_SNP_PASS_derived.vcf.gz'
    created = []
    if not manifest_only:
        if shutil.which('bcftools') is None:
            die('bcftools not found in PATH (it was previously installed in WSL)')
        if not allvcf.is_file() or not passvcf.is_file():
            die('Required local Step 17 combined VCF files not found')
        for src, expected_n in [(allvcf,2471),(passvcf,1824)]:
            sample_ids = ids_from_vcf(src)
            if len(sample_ids)!=282 or set(sample_ids)!=all_ids:
                die(f'VCF vs 282-sample manifest mismatch: {src}')
            nsites = count_from_index(src)
            if nsites!=expected_n:
                die(f'Original VCF has unexpected site count {nsites}, expected {expected_n}: {src}')
        vcf_out.mkdir(parents=True, exist_ok=True)
        vcf_tasks=[(allvcf,'PI4K_2020_119_all_sites.vcf.gz','COMBINED',2471),
                   (passvcf,'PI4K_2020_119_PASS_sites.vcf.gz','COMBINED',1824)]
        for state in STATES:
            vcf_tasks.append((passvcf,f'PI4K_2020_{state}_PASS_sites.vcf.gz',state,1824))
        for src, name, group, expected_sites in vcf_tasks:
            ids = ids_combined if group=='COMBINED' else sorted(df2020.loc[df2020.state==group,'sample_id'])
            record = subset_vcf(src,ids,list_paths[group],vcf_out/name,expected_sites)
            created.append(record)
            print(f"{group}: {record['samples']} samples, {record['sites']} sites, {record['action']} -- {name}")

    result = {
        'status':'COHORT_CREATED_AND_LOCAL_VCFS_VERIFIED' if not manifest_only else 'MANIFEST_CREATED_ONLY',
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'source_study':STUDY, 'source_year':int(YEAR),
        'original_n':282,'same_study_year_n':119,'n_excluded_from_controlled_cohort_due_to_study_year_only':163,
        'state_counts':YEAR_COUNTS,
        'cohort_scenario_counts_by_state':summary_tbl.to_dict(orient='records'),
        'source_step20_manifest':str(manifest),'source_step20_manifest_sha256':sha256(manifest),
        'sample_list_sha256':{k:sha256(p) for k,p in list_paths.items()},
        'vcf_outputs':created,
        'caveats':[
            'Study/year restriction controls these two provenance variables but does not prove other sampling confounders absent.',
            'This is a same-study-year cohort, NOT pairs of biologically matched participants.',
            'The full 282-sample cohort should be retained for descriptive and sensitivity analyses.',
            'The 119 cohort is proposed as the primary geographically comparable comparison set, pending assessment of site comparability and quality.',
            'Fws>=0.95 and sample PASS-site missingness<=10% are exploratory sensitivity flags, not applied exclusions.',
            'VCFs preserve the SNP-only derived 2-slot genotype encoding, not biological diploidy.',
            'Presence of 1824 globally PASS SNP sites does not mean 1824 sites are polymorphic in this cohort.',
            'No new genomic data downloaded or genotypes modified at source.'
        ]
    }
    path=out/'PI4K_2020_cohort_build_report.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'SUCCESS: Cohort manifests created; report: {path}')
    print('NO remote requests, NO original genotype files modified.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    p.add_argument('--manifest-only',action='store_true',help='Create validated sample lists and audit; do not run bcftools')
    a=p.parse_args()
    try:
        main(a.root.resolve(),a.manifest_only)
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(1)
