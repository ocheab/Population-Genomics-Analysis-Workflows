#!/usr/bin/env python3
"""Step 19: audit-ONLY integration of Pf8 Fws, provenance and local PI4K QC.

Requires pandas and numpy (installed with the PI4K malariagen-data venv).
No network access, no NPZ changes, no exclusions, no inferential statistics.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import sys
import numpy as np
import pandas as pd

ROOT = Path('/mnt/c/PI4K_Nigeria_Study')
STATES = [('Oyo',71),('Ogun',38),('Lagos',129),('Osun',15),('Ondo',29)]
COUNTS = dict(STATES)
OUT = ROOT / 'results/Southwest5/Step19_Fws_Provenance'
FILES = {
    'metadata': ROOT / 'results/Southwest5/Southwest5_QC_pass_samples.tsv',
    'fws': ROOT / 'metadata/Pf8_fws.tsv',
    'sample_qc': ROOT / 'results/Southwest5/Step18_QC/PI4K_QC_by_sample.tsv',
    'site_qc': ROOT / 'results/Southwest5/Step18_QC/PI4K_QC_by_site.tsv',
}


def read_table(path):
    if not path.is_file():
        raise FileNotFoundError(f'Required input not found: {path}')
    df = pd.read_csv(path, sep='\t', encoding='utf-8-sig', dtype=str, keep_default_na=False)
    df.columns = [str(c).strip().lstrip('\ufeff') for c in df.columns]
    return df


def require_columns(df, columns, file_label):
    absent = sorted(set(columns) - set(df.columns))
    if absent:
        raise ValueError(f'{file_label}: missing required columns {absent}; observed {df.columns.tolist()}')


def check_unique(df, cols, label):
    bad = df.duplicated(cols, keep=False)
    if bad.any():
        example = df.loc[bad, cols].head(5).to_dict('records')
        raise ValueError(f'{label} has duplicated keys {cols}; examples: {example}')


def number_col(frame, name):
    return pd.to_numeric(frame[name].replace('', np.nan), errors='coerce')


def pctl(x, q):
    x = pd.to_numeric(x, errors='coerce').dropna()
    return float(x.quantile(q)) if len(x) else None


def md5_or_sha(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def run():
    print('STEP 19: Integrate Fws, provenance and PI4K sample QC (offline; no exclusions)', flush=True)
    md = read_table(FILES['metadata'])
    fw = read_table(FILES['fws'])
    sq = read_table(FILES['sample_qc'])
    site = read_table(FILES['site_qc'])
    require_columns(md, ['Sample','Study','Country','Admin level 1','Year','All samples same case', '% callable', 'QC pass'], 'metadata')
    require_columns(fw, ['Sample','Fws'], 'Pf8_fws.tsv')
    require_columns(sq, ['sample_id','state','site_set','n_sites','fully_observed_calls','fully_missing_calls','missing_slot_fraction','fully_observed_call_fraction','discordant_two_slot_calls','discordant_two_slot_fraction_among_full_calls'], 'PI4K_QC_by_sample.tsv')
    require_columns(site, ['chrom','pos_1based','pf8_filter_pass','two_or_more_allele_codes_observed','fully_observed_call_fraction'], 'PI4K_QC_by_site.tsv')
    check_unique(md, ['Sample'], 'selected metadata')
    check_unique(fw, ['Sample'], 'Pf8 Fws')
    check_unique(sq, ['sample_id','site_set'], 'PI4K sample QC')
    check_unique(site, ['chrom','pos_1based'], 'PI4K site QC')

    md = md[(md['Country']=='Nigeria') & md['Admin level 1'].isin(COUNTS)].copy()
    md = md[md['QC pass'].str.lower().eq('true')].copy()
    if len(md) != 282:
        raise ValueError(f'Expected 282 selected metadata rows, got {len(md)}')
    md['state'] = md['Admin level 1']
    md['sample_id'] = md['Sample'].str.strip()
    if (md['sample_id']=='').any():
        raise ValueError('Metadata contains blank sample IDs')
    md_by_state = md['state'].value_counts().to_dict()
    if md_by_state != COUNTS:
        raise ValueError(f'Metadata state sample counts mismatch: {md_by_state} vs {COUNTS}')

    list_rows=[]
    for state, n in STATES:
        p = ROOT / 'results/Southwest5/sample_lists' / f'{state}_sample_IDs.txt'
        if not p.exists(): raise FileNotFoundError(p)
        ids = [line.strip().lstrip('\ufeff') for line in p.read_text(encoding='utf-8-sig').splitlines() if line.strip()]
        if len(ids) != n or len(set(ids)) != n:
            raise ValueError(f'{state} expected {n} unique IDs in {p}, got {len(ids)} rows, {len(set(ids))} unique')
        list_rows.extend((sample, state) for sample in ids)
    checklist = pd.DataFrame(list_rows, columns=['sample_id','state'])
    check_unique(checklist, ['sample_id'], 'sample ID lists')
    if set(checklist.sample_id) != set(md.sample_id):
        raise ValueError(f'Manifest mismatch: {len(set(checklist.sample_id)-set(md.sample_id))} IDs missing from metadata; {len(set(md.sample_id)-set(checklist.sample_id))} metadata IDs not in list')
    compare = checklist.merge(md[['sample_id','state']], on='sample_id', suffixes=('_listed','_meta'), validate='1:1')
    if not (compare.state_listed==compare.state_meta).all():
        raise ValueError('Sample lists disagree with state metadata')

    for scope, nsite in [('ALL',2471),('PF8_PASS',1824)]:
        f = sq.loc[sq['site_set']==scope]
        if len(f)!=282 or not (number_col(f,'n_sites')==nsite).all():
            raise ValueError(f'Unexpected {scope} QC row/site counts: {len(f)} rows')
        if set(f['sample_id']) != set(checklist['sample_id']):
            raise ValueError(f'{scope}: QC sample IDs do not match manifest')
        if not (f.set_index('sample_id').loc[compare.sample_id,'state'].to_numpy() == compare.state_listed.to_numpy()).all():
            raise ValueError(f'{scope}: QC state assignments differ from manifest')
    if len(site)!=2471:
        raise ValueError(f'Expected 2471 site QC rows, got {len(site)}')
    pass_site = site['pf8_filter_pass'].str.lower().isin(['1','true','yes'])
    if int(pass_site.sum()) != 1824:
        raise ValueError(f'Expected 1824 PASS sites; got {int(pass_site.sum())}')

    fw['sample_id'] = fw['Sample'].str.strip()
    fw['Fws_numeric'] = number_col(fw, 'Fws')
    meta_columns = ['sample_id','state','Study','Year','All samples same case','% callable','ENA','Sample type','Sample was in Pf7']
    meta_columns = [c for c in meta_columns if c in md.columns]
    integrated = checklist.merge(md[meta_columns], on=['sample_id','state'], validate='1:1')
    integrated = integrated.merge(fw[['sample_id','Fws_numeric']], on='sample_id', how='left', validate='1:1')
    integrated = integrated.rename(columns={'Fws_numeric':'Fws', 'Study':'study','Year':'collection_year_raw','All samples same case':'case_group_id','% callable':'genome_percent_callable','ENA':'ena_accession','Sample type':'sample_type','Sample was in Pf7':'in_pf7'})
    for scope, prefix in [('ALL','all'),('PF8_PASS','pass')]:
        f=sq.loc[sq['site_set']==scope].copy().drop(columns=['site_set','state'])
        f=f.rename(columns={c:f'{prefix}_{c}' for c in f.columns if c!='sample_id'})
        integrated=integrated.merge(f,on='sample_id', validate='1:1')
    integrated['collection_year'] = number_col(integrated,'collection_year_raw')
    integrated['collection_year'] = integrated['collection_year'].apply(lambda x: int(x) if pd.notna(x) and x==int(x) else pd.NA).astype('Int64')
    integrated['genome_percent_callable'] = number_col(integrated,'genome_percent_callable')
    for col in integrated.columns:
        if col.startswith(('all_','pass_')) and col!='sample_id':
            integrated[col]=number_col(integrated,col)

    # Descriptive sensitivity flags only. Never drop rows here.
    integrated['fws_available'] = integrated['Fws'].notna()
    integrated['fws_ge_0_95_exploratory'] = integrated['Fws']>=0.95
    integrated['fws_ge_0_90_exploratory'] = integrated['Fws']>=0.90
    integrated['pass_missing_le_0_05_exploratory'] = integrated['pass_missing_slot_fraction']<=0.05
    integrated['pass_missing_le_0_10_exploratory'] = integrated['pass_missing_slot_fraction']<=0.10
    integrated['candidate_fws95_missing10_sensitivity_only'] = integrated['fws_ge_0_95_exploratory'] & integrated['pass_missing_le_0_10_exploratory']
    integrated['case_group_id'] = integrated['case_group_id'].astype(str).str.strip()
    integrated['valid_case_group'] = ~integrated['case_group_id'].isin(['', 'nan','None','NA'])

    if integrated['Fws'].notna().any():
        outliers=integrated['Fws'].notna() & ~integrated['Fws'].between(0,1)
        if outliers.any(): print(f'CAUTION: {outliers.sum()} Fws values outside [0,1], inspect original data.',flush=True)
    state_stats=[]
    for state, _ in STATES:
        d=integrated.loc[integrated.state==state]
        fv=d.Fws.dropna()
        grp=d.loc[d.valid_case_group].groupby('case_group_id').size()
        state_stats.append({
            'state':state,'n_qc_pass':len(d),
            'n_fws_available':len(fv),'n_fws_missing':int(d.Fws.isna().sum()),
            'fws_median':pctl(fv,.5),'fws_q1':pctl(fv,.25),'fws_q3':pctl(fv,.75),
            'n_fws_ge_0_95':int((fv>=.95).sum()),'n_fws_lt_0_95':int((fv<.95).sum()),
            'n_fws_ge_0_90':int((fv>=.90).sum()),
            'pass_missing_median':pctl(d.pass_missing_slot_fraction,.5),
            'pass_missing_q1':pctl(d.pass_missing_slot_fraction,.25),
            'pass_missing_q3':pctl(d.pass_missing_slot_fraction,.75),
            'n_pass_missing_gt_0_05':int((d.pass_missing_slot_fraction>.05).sum()),
            'n_pass_missing_gt_0_10':int((d.pass_missing_slot_fraction>.10).sum()),
            'n_pass_missing_gt_0_20':int((d.pass_missing_slot_fraction>.20).sum()),
            'n_fws95_and_missing10_sensitivity':int(d.candidate_fws95_missing10_sensitivity_only.sum()),
            'distinct_studies':int(d.study.nunique()),
            'distinct_years_known':int(d.collection_year.nunique()),
            'n_year_missing':int(d.collection_year.isna().sum()),
            'n_case_group_missing':int((~d.valid_case_group).sum()),
            'n_repeated_case_group_ids':int((grp>1).sum()),
            'n_samples_in_repeated_case_groups':int(grp[grp>1].sum()),
            'genome_percent_callable_median':pctl(d.genome_percent_callable,.5),
        })
    state_df=pd.DataFrame(state_stats)

    counts=integrated.loc[integrated.valid_case_group].groupby('case_group_id',dropna=False).agg(n_samples=('sample_id','size'),states=('state',lambda s:';'.join(sorted(set(s)))),sample_ids=('sample_id',lambda s:';'.join(sorted(s))),studies=('study',lambda s:';'.join(sorted(set(s))))).reset_index()
    repeated_cases=counts.loc[counts.n_samples>1].sort_values(['n_samples','case_group_id'],ascending=[False,True])
    study_year=(integrated.groupby(['state','study','collection_year'],dropna=False).agg(n_samples=('sample_id','size'),n_fws_available=('fws_available','sum')).reset_index().sort_values(['state','study','collection_year']))
    site2=site.copy()
    for col in ['fully_observed_call_fraction','two_or_more_allele_codes_observed']:
        site2[col]=number_col(site2,col)
    pass_sites=site2.loc[pass_site]
    n_pass_variable=int((pass_sites['two_or_more_allele_codes_observed']==1).sum())
    summary={
        'status':'INTEGRATED_QC_COMPLETE_NO_EXCLUSIONS',
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'selected_samples':len(integrated),'state_counts':COUNTS,
        'fws_available':int(integrated.fws_available.sum()),
        'fws_missing':int(integrated.Fws.isna().sum()),
        'fws_ge_0_95':int(integrated.fws_ge_0_95_exploratory.sum()),
        'fws_lt_0_95_among_available':int(((integrated.Fws<0.95)&integrated.fws_available).sum()),
        'pass_missing_gt_0_10':int((integrated.pass_missing_slot_fraction>.1).sum()),
        'fws95_and_pass_missing_le_0_10_sensitivity_only':int(integrated.candidate_fws95_missing10_sensitivity_only.sum()),
        'repeated_case_group_ids':int(len(repeated_cases)),
        'unique_studies':int(integrated.study.nunique()),
        'known_collection_years':sorted(int(y) for y in integrated.collection_year.dropna().unique()),
        'sites_all':len(site),'sites_pf8_pass':int(pass_site.sum()),
        'pass_sites_two_or_more_allele_codes':n_pass_variable,
        'files_sha256':{k:md5_or_sha(v) for k,v in FILES.items()},
        'interpretation_caveats':[
            'No exclusions made. All thresholds are descriptive and/or potential sensitivity cohorts, not automatically justified filters.',
            'Fws >= 0.95 is an exploratory low-within-host-diversity threshold, not proof of monoclonality.',
            'Repeat case-group IDs flag possible non-independent samples; review before cross-validation.',
            'Year/study may confound state classification. Splitting by sample alone can leak provenance.',
            'Two genotype slots do not constitute two independent haploid isolates.',
            'PASS is a Pf8 variant-site flag, not a per-sample genotype quality metric.',
            'SNP-only callset lacks indels and original per-genotype quality/depth.',
            '197 polymorphic-code candidate sites are not yet validated independent-isolate polymorphisms.',
        ]
    }
    OUT.mkdir(parents=True, exist_ok=True)
    integrated.to_csv(OUT/'PI4K_integrated_sample_QC.tsv',sep='\t',index=False,na_rep='')
    state_df.to_csv(OUT/'PI4K_Fws_QC_by_state.tsv',sep='\t',index=False,na_rep='')
    study_year.to_csv(OUT/'PI4K_provenance_by_state_study_year.tsv',sep='\t',index=False,na_rep='')
    repeated_cases.to_csv(OUT/'PI4K_repeated_case_groups.tsv',sep='\t',index=False,na_rep='')
    (OUT/'PI4K_integrated_QC_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('SUCCESS: Integrated data for 282 samples; no genotype downloads and no filtering.',flush=True)
    print(f'Fws available: {summary["fws_available"]}; missing: {summary["fws_missing"]}',flush=True)
    print(f'Fws >= 0.95: {summary["fws_ge_0_95"]} (exploratory only)',flush=True)
    print(f'PASS-site sample missingness >10%: {summary["pass_missing_gt_0_10"]}',flush=True)
    print(f'Repeated case-group identifiers: {summary["repeated_case_group_ids"]}',flush=True)
    print(f'Outputs: {OUT}',flush=True)
    print('\nSTATE AUDIT')
    print(state_df[['state','n_qc_pass','n_fws_available','n_fws_ge_0_95','pass_missing_median','n_pass_missing_gt_0_10','distinct_studies','distinct_years_known','n_repeated_case_group_ids']].to_string(index=False))

if __name__=='__main__':
    try:
        run()
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(1)
