#!/usr/bin/env python3
"""PI4K Southwest5 Step 20: offline provenance overlap and cohort audit.

READ-ONLY inputs, no genotype reads, no filtering of source files, no inference.
Requires pandas and numpy (already present in pi4k_env).
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

STATES = ('Oyo', 'Ogun', 'Lagos', 'Osun', 'Ondo')
EXPECTED = {'Oyo': 71, 'Ogun': 38, 'Lagos': 129, 'Osun': 15, 'Ondo': 29}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f'Missing Step 19 input: {path}')
    df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, encoding='utf-8-sig')
    df.columns = [str(c).lstrip('\ufeff').strip() for c in df.columns]
    return df


def truth(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.strip().str.lower()
    unknown = ~values.isin(['true','false','1','0'])
    if unknown.any():
        raise ValueError(f'Unrecognised boolean values {sorted(values[unknown].unique())[:5]}')
    return values.isin(['true','1'])


def clean_year(y: str) -> str:
    s = str(y).strip()
    if not s or s.lower() in ('nan','none','na','<na>'):
        return 'UNKNOWN'
    try:
        f = float(s)
        if f.is_integer() and 1900 <= f <= 2100:
            return str(int(f))
    except ValueError:
        pass
    raise ValueError(f'Invalid collection year: {s!r}')


def stratum_table(df: pd.DataFrame, key_cols: list[str], value='n_samples') -> pd.DataFrame:
    tab = (df.groupby(key_cols + ['state'], dropna=False).size()
           .unstack('state', fill_value=0).reindex(columns=STATES, fill_value=0))
    tab['Total'] = tab.sum(axis=1)
    return tab.reset_index()


def run(root: Path) -> None:
    src = root / 'results/Southwest5/Step19_Fws_Provenance/PI4K_integrated_sample_QC.tsv'
    summary_src = root / 'results/Southwest5/Step19_Fws_Provenance/PI4K_integrated_QC_summary.json'
    df = load(src)
    needed = {'sample_id','state','study','collection_year','case_group_id','Fws',
              'pass_missing_slot_fraction'}
    absent = sorted(needed - set(df.columns))
    if absent:
        raise ValueError(f'Step 19 integrated table missing columns: {absent}')
    if len(df) != 282 or df['sample_id'].duplicated().any() or (df['sample_id'].str.strip()=='').any():
        raise ValueError('Expected exactly 282 distinct nonblank sample IDs')
    if df['state'].value_counts().to_dict() != EXPECTED:
        raise ValueError(f'State counts changed: {df.state.value_counts().to_dict()}')

    df['study'] = df['study'].str.strip().replace('', 'UNKNOWN')
    df['year'] = df['collection_year'].map(clean_year)
    df['Fws'] = pd.to_numeric(df['Fws'].replace('', np.nan), errors='raise')
    df['pass_missing_slot_fraction'] = pd.to_numeric(df['pass_missing_slot_fraction'].replace('', np.nan), errors='raise')
    if df['Fws'].isna().any() or df['pass_missing_slot_fraction'].isna().any():
        raise ValueError('Fws or PASS-site missingness unexpectedly unavailable')
    if not df['Fws'].between(0,1).all() or not df['pass_missing_slot_fraction'].between(0,1).all():
        raise ValueError('Fws or missingness outside [0,1]')
    if 'fws_ge_0_95_exploratory' in df:
        declared = truth(df['fws_ge_0_95_exploratory'])
        if not (declared == df['Fws'].ge(.95)).all():
            raise ValueError('Step 19 Fws flag disagrees with numeric Fws')
    if 'candidate_fws95_missing10_sensitivity_only' in df:
        declared = truth(df['candidate_fws95_missing10_sensitivity_only'])
        expectedflag = df['Fws'].ge(.95) & df['pass_missing_slot_fraction'].le(.10)
        if not (declared == expectedflag).all():
            raise ValueError('Step 19 joint flag disagrees with numeric measurements')

    df['cohort_full_282'] = True
    df['cohort_missing_le10'] = df['pass_missing_slot_fraction'].le(.10)
    df['cohort_missing_le05'] = df['pass_missing_slot_fraction'].le(.05)
    df['cohort_fws_ge095'] = df['Fws'].ge(.95)
    df['cohort_fws_ge090'] = df['Fws'].ge(.90)
    df['cohort_joint_fws095_missing10'] = df['cohort_fws_ge095'] & df['cohort_missing_le10']
    df['cohort_joint_fws090_missing10'] = df['cohort_fws_ge090'] & df['cohort_missing_le10']

    observed = {'n_full':int(df.cohort_full_282.sum()),
                'n_missing_le10':int(df.cohort_missing_le10.sum()),
                'n_fws095':int(df.cohort_fws_ge095.sum()),
                'n_joint_fws095_missing10':int(df.cohort_joint_fws095_missing10.sum())}
    expected_totals = {'n_full':282, 'n_missing_le10':248, 'n_fws095':196, 'n_joint_fws095_missing10':166}
    if observed != expected_totals:
        raise ValueError(f'Cohort counts differ from Step 19 summary: {observed}, expected {expected_totals}')
    if summary_src.is_file():
        prior = json.loads(summary_src.read_text(encoding='utf-8'))
        checks = {'selected_samples':282,'fws_ge_0_95':196,
                  'fws95_and_pass_missing_le_0_10_sensitivity_only':166}
        for k,v in checks.items():
            if prior.get(k) != v:
                raise ValueError(f'Previous summary {k} changed: {prior.get(k)!r} vs {v}')

    out = root / 'results/Southwest5/Step20_Provenance_Cohorts'
    out.mkdir(parents=True, exist_ok=True)

    cols = ['sample_id','state','study','year','case_group_id','Fws','pass_missing_slot_fraction',
            'cohort_full_282','cohort_missing_le10','cohort_missing_le05',
            'cohort_fws_ge095','cohort_fws_ge090','cohort_joint_fws095_missing10',
            'cohort_joint_fws090_missing10']
    manifest = df[cols].sort_values(['state','sample_id']).copy()
    manifest.to_csv(out/'PI4K_cohort_membership_audit_282.tsv', sep='\t', index=False)
    cohort_cols = [c for c in cols if c.startswith('cohort_')]
    cohort_data = []
    for state in (*STATES, 'TOTAL'):
        s = df if state == 'TOTAL' else df[df['state']==state]
        row = {'state':state, 'n_all':len(s), 'n_fws_lt095':int((s.Fws<.95).sum()),
               'n_pass_missing_gt10':int((s.pass_missing_slot_fraction>.10).sum()),
               'median_Fws':float(s.Fws.median()),
               'median_PASS_missingness':float(s.pass_missing_slot_fraction.median()),
               'studies':int(s.study.nunique()), 'years':int(s.year.nunique())}
        row.update({c:int(s[c].sum()) for c in cohort_cols})
        cohort_data.append(row)
    cohort = pd.DataFrame(cohort_data)
    cohort.to_csv(out/'PI4K_cohort_counts_by_state.tsv', sep='\t', index=False, float_format='%.6f')

    study_tab = stratum_table(df, ['study'])
    year_tab = stratum_table(df, ['year'])
    study_year_tab = stratum_table(df, ['study','year'])
    study_tab.to_csv(out/'PI4K_state_by_study.tsv', sep='\t', index=False)
    year_tab.to_csv(out/'PI4K_state_by_year.tsv', sep='\t', index=False)
    study_year_tab.to_csv(out/'PI4K_state_by_study_year.tsv', sep='\t', index=False)

    overlap_rows = []
    for a,b in itertools.combinations(STATES,2):
        a_df = df[df.state==a]
        b_df = df[df.state==b]
        a_years,b_years=set(a_df.year)-{'UNKNOWN'},set(b_df.year)-{'UNKNOWN'}
        a_studies,b_studies=set(a_df.study)-{'UNKNOWN'},set(b_df.study)-{'UNKNOWN'}
        shared_years=sorted(a_years & b_years)
        shared_studies=sorted(a_studies & b_studies)
        ay = a_df.groupby(['study','year']).size().to_dict()
        by = b_df.groupby(['study','year']).size().to_dict()
        matched = sorted((k,ay[k],by[k]) for k in ay.keys() & by.keys()
                         if k[0]!='UNKNOWN' and k[1]!='UNKNOWN')
        overlap_rows.append({
            'state_a':a,'state_b':b,'samples_a':len(a_df),'samples_b':len(b_df),
            'shared_collection_years':';'.join(shared_years),
            'number_shared_years':len(shared_years),
            'shared_studies':';'.join(shared_studies),
            'number_shared_studies':len(shared_studies),
            'shared_study_year_strata':';'.join(f'{st}@{yr}' for (st,yr),na,nb in matched),
            'n_shared_study_year_strata':len(matched),
            'samples_a_in_shared_study_year_strata':sum(na for key,na,nb in matched),
            'samples_b_in_shared_study_year_strata':sum(nb for key,na,nb in matched),
            'max_one_to_one_matched_pairs':sum(min(na,nb) for key,na,nb in matched),
            'flag_no_shared_study_year':len(matched)==0,
        })
    overlaps=pd.DataFrame(overlap_rows)
    overlaps.to_csv(out/'PI4K_pairwise_provenance_overlap.tsv', sep='\t', index=False)
    # Study-year overlaps by cohort, not just all-sample distribution.
    sens = stratum_table(df[df.cohort_joint_fws095_missing10], ['study','year'])
    sens.to_csv(out/'PI4K_sensitivity_by_study_year.tsv', sep='\t', index=False)
    # Verify case group keys remain unique. No record deletion.
    case = df['case_group_id'].fillna('').astype(str).str.strip()
    valid = ~case.str.lower().isin(['','none','nan','na'])
    repeats = case[valid].value_counts()
    repeat_ids = int((repeats>1).sum())
    if repeat_ids:
        print(f'CAUTION: {repeat_ids} repeated case-group keys found; review independence.', flush=True)

    n_no_shared = int(overlaps['flag_no_shared_study_year'].sum())
    warnings = [
        'No genotype data retrieved or modified; cohort membership is descriptive only.',
        'Do not treat combined Fws>=0.95 and missingness<=0.10 as a proven primary cohort.',
        'No matching study-year strata does not prove absence of geography effects; it means direct adjustment is limited.',
        'Collection time/study provenance may be inseparable from geography in some contrasts.',
        'Two genotype slots are not two independent haploid parasite chromosomes.',
        'SNP-only data exclude indels and lack original per-genotype DP and GQ values.',
        'The 197 preliminary variable-code sites have not been independently validated for state population-genetics testing.',
        'This audit does not implement cross-validation, test state predictability or make inferential claims.'
    ]
    summary = {
        'status':'PROVENANCE_AND_COHORT_AUDIT_COMPLETE_NO_EXCLUSIONS',
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'source_integrated_QC_sha256':sha256(src),
        'original_samples':282,'state_counts':EXPECTED,
        'cohort_totals':observed,'source_repeated_case_group_count':repeat_ids,
        'total_state_pairs':len(overlaps),'state_pairs_without_shared_study_year':n_no_shared,
        'state_pairs_with_shared_study_year':len(overlaps)-n_no_shared,
        'study_count':int(df.study.nunique()),
        'collection_years':sorted(set(df.year)-{'UNKNOWN'}),
        'outputs': [p.name for p in sorted(out.glob('*.tsv'))],
        'interpretation_caveats':warnings,
    }
    (out/'PI4K_step20_design_audit.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print('SUCCESS: Step 20 provenance/cohort audit completed OFFLINE; NO data filtered.')
    print('All samples:', len(df), '| PASS missing<=10%:', observed['n_missing_le10'],
          '| Fws>=0.95:', observed['n_fws095'],
          '| joint sensitivity:', observed['n_joint_fws095_missing10'])
    print('\nSTATE COHORT AUDIT')
    print(cohort[['state','n_all','cohort_missing_le10','cohort_fws_ge095',
                  'cohort_joint_fws095_missing10','studies','years']].to_string(index=False))
    print('\nPROVENANCE OVERLAP')
    print(f'Pairs with no shared study-year stratum: {n_no_shared} / {len(overlaps)}')
    print(overlaps[['state_a','state_b','n_shared_study_year_strata',
                    'max_one_to_one_matched_pairs']].to_string(index=False))
    print('\nOutputs:', out)
    print('Share PI4K_state_by_study_year.tsv, PI4K_pairwise_provenance_overlap.tsv and the JSON report.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent,
                        help='Study root (default: directory above scripts/)')
    args = parser.parse_args()
    try:
        run(args.root.resolve())
    except Exception as e:
        print(f'ERROR: {type(e).__name__}: {e}', file=sys.stderr)
        sys.exit(1)
