#!/usr/bin/env python3
"""Step 18: Offline, non-destructive QC audit of extracted MalariaGEN Pf8 PI4K SNPs.

Input: Step 16 verified NPZ. This script makes descriptive QC tables only;
no exclusion, imputation, FST, PCA, biological haplotype inference, or network access.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

STATES = {"Oyo": 71, "Ogun": 38, "Lagos": 129, "Osun": 15, "Ondo": 29}
CHROM = "Pf3D7_05_v3"
START, END = 409403, 415758
ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "variants" / "Southwest5" / "Sanger_SNP_complete"
SRC = SRC_DIR / "PI4K_Southwest5_282_Sanger_SNP_all_sites.npz"
OUT = ROOT / "results" / "Southwest5" / "Step18_QC"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, header: list[str], rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def safe_float(n: int | float, d: int | float) -> float | None:
    return float(n / d) if d else None


def flt(n: int | float, d: int | float) -> str:
    v = safe_float(n, d)
    return "NA" if v is None else f"{v:.6f}"


def allele_code_validity(gt: np.ndarray, alleles: np.ndarray):
    if np.any(gt < -1):
        raise ValueError("Genotype codes below -1 observed")
    if gt.ndim != 3 or gt.shape[-1] != 2:
        raise ValueError("Expected 3D 2-slot genotypes")
    for i in range(gt.shape[0]):
        al = [str(a).strip() for a in alleles[i]]
        def empty(s: str) -> bool:
            return s.lower() in ("", "none", "null", "nan", ".")
        while len(al) > 1 and empty(al[-1]):
            al.pop()
        if not al or not al[0] or any(empty(a) for a in al):
            raise ValueError(f"Missing or noncontiguous allele strings at SNP index {i}")
        if np.any(gt[i] >= len(al)):
            raise ValueError(f"Allele code out of range at SNP index {i}")


def metrics_for_calls(gt: np.ndarray) -> dict[str, int | str]:
    """gt: n_sites x n_samples x 2. Both slots are treated as encoded slots, not independent parasite chromosomes."""
    called = gt >= 0
    full = np.all(called, axis=-1)
    no = np.all(~called, axis=-1)
    partial = ~(full | no)
    discordant = full & (gt[:, :, 0] != gt[:, :, 1])
    alt = np.any(gt > 0, axis=-1)
    ns, nsm = gt.shape[:2]
    return dict(
        sites=ns,
        samples=nsm,
        slots_total=int(gt.size),
        slots_missing=int(np.sum(~called)),
        calls_fully_observed=int(full.sum()),
        calls_partially_observed=int(partial.sum()),
        calls_fully_missing=int(no.sum()),
        calls_discordant_two_slot=int(discordant.sum()),
        calls_with_alt_slot=int(alt.sum()),
        fraction_missing_slots=flt(int(np.sum(~called)),gt.size),
        fraction_fully_observed_calls=flt(int(full.sum()),ns*nsm),
        fraction_discordant_among_fully_observed=flt(int(discordant.sum()),int(full.sum())),
    )


def run():
    if not SRC.is_file():
        raise FileNotFoundError(f"Extracted genotype input absent: {SRC}")
    with np.load(SRC, allow_pickle=False) as z:
        ids = z["sample_ids"].astype(str)
        states = z["states"].astype(str)
        positions = z["positions"].astype(np.int32)
        passmask = z["filter_pass"].astype(bool)
        alleles = z["alleles"].astype(str)
        gt = z["genotypes"].astype(np.int8)
    if (gt.shape != (2471,282,2) or len(ids) != 282 or len(set(ids)) != 282 or
        len(positions) != 2471 or int(passmask.sum()) != 1824 or
        np.any(np.diff(positions) <= 0) or positions[0] < START or positions[-1] > END):
        raise RuntimeError("Input genotype shape/site filter/locus fails integrity requirements")
    if alleles.shape[0] != len(positions) or len(states) != len(ids):
        raise RuntimeError("Input allele/sample state arrays do not align")
    for state,n in STATES.items():
        idx = np.flatnonzero(states == state)
        if len(idx) != n:
            raise RuntimeError(f"State {state} count mismatch: {len(idx)} != {n}")
        idfile = ROOT / "results" / "Southwest5" / "clean_sample_lists" / f"{state}_sample_IDs.txt"
        if not idfile.exists():
            idfile = ROOT / "results" / "Southwest5" / "sample_lists" / f"{state}_sample_IDs.txt"
        if not idfile.exists():
            raise FileNotFoundError(f"Missing verified sample list: {idfile}")
        listed = [v.strip() for v in idfile.read_text(encoding="utf-8-sig").splitlines() if v.strip()]
        if ids[idx].tolist() != listed:
            raise RuntimeError(f"ID order does not match verified list for {state}")
    allele_code_validity(gt, alleles)
    OUT.mkdir(parents=True, exist_ok=True)

    # Sample-level descriptive QC, all-sites and PASS-sites, no exclusions.
    sample_header = ["sample_id","state","site_set","n_sites","fully_observed_calls","partially_observed_calls",
                     "fully_missing_calls","missing_slots","total_slots","missing_slot_fraction",
                     "fully_observed_call_fraction","discordant_two_slot_calls",
                     "discordant_two_slot_fraction_among_full_calls","calls_with_alt_slot"]
    sample_rows=[]
    for tag,mask in (("ALL",np.ones(len(positions),bool)),("PF8_PASS",passmask)):
        for j,(sid,st) in enumerate(zip(ids,states)):
            m=metrics_for_calls(gt[mask,j:j+1,:])
            sample_rows.append([sid,st,tag,m["sites"],m["calls_fully_observed"],m["calls_partially_observed"],
                                m["calls_fully_missing"],m["slots_missing"],m["slots_total"],
                                m["fraction_missing_slots"],m["fraction_fully_observed_calls"],
                                m["calls_discordant_two_slot"],m["fraction_discordant_among_fully_observed"],
                                m["calls_with_alt_slot"]])
    write_tsv(OUT/"PI4K_QC_by_sample.tsv",sample_header,sample_rows)

    # Per state: site counts, encoded missingness, no genotype exclusions.
    state_header=["state","site_set","samples","sites","total_sample_site_calls","fully_observed_calls",
                  "partially_observed_calls","fully_missing_calls","missing_slots","total_slots",
                  "missing_slot_fraction","fully_observed_call_fraction","discordant_two_slot_calls",
                  "discordant_two_slot_fraction_among_full_calls","sites_with_any_alt_slot",
                  "sites_with_two_or_more_alleles_observed"]
    state_rows=[]
    for state in STATES:
        sm = states==state
        for tag,mask in (("ALL",np.ones(len(positions),bool)),("PF8_PASS",passmask)):
            sub=gt[mask][:,sm,:]
            m=metrics_for_calls(sub)
            n_alt=int(np.any(sub>0,axis=(1,2)).sum())
            n_var=0
            for site in sub:
                vals=site[site>=0]
                if vals.size and np.unique(vals).size>=2:
                    n_var+=1
            state_rows.append([state,tag,int(sm.sum()),m["sites"],m["sites"]*m["samples"],
                               m["calls_fully_observed"],m["calls_partially_observed"],m["calls_fully_missing"],
                               m["slots_missing"],m["slots_total"],m["fraction_missing_slots"],
                               m["fraction_fully_observed_calls"],m["calls_discordant_two_slot"],
                               m["fraction_discordant_among_fully_observed"],n_alt,n_var])
    write_tsv(OUT/"PI4K_QC_by_state.tsv",state_header,state_rows)

    # Per site, whole five-state cohort.
    site_header=["chrom","pos_1based","pf8_filter_pass","n_samples","fully_observed_calls",
                 "partially_observed_calls","fully_missing_calls","fully_observed_call_fraction",
                 "discordant_two_slot_calls","any_alt_slot_in_cohort","distinct_allele_codes_observed",
                 "two_or_more_allele_codes_observed","exploratory_95pct_full_call_threshold"]
    site_rows=[]
    pass_variants=0
    pass_with_alt=0
    for i,p in enumerate(positions):
        calls=gt[i,:,:]
        called=calls>=0
        fully=np.all(called,axis=1)
        missing=np.all(~called,axis=1)
        partial=~(fully|missing)
        disc=fully & (calls[:,0]!=calls[:,1])
        obs=calls[called]
        uniq=int(np.unique(obs).size)
        alt=bool(np.any(obs>0))
        if passmask[i] and uniq>=2: pass_variants +=1
        if passmask[i] and alt: pass_with_alt +=1
        site_rows.append([CHROM,int(p),int(passmask[i]),len(ids),int(fully.sum()),int(partial.sum()),
                          int(missing.sum()),flt(int(fully.sum()),len(ids)),int(disc.sum()),
                          int(alt),uniq,int(uniq>=2),int(fully.mean()>=0.95)])
    write_tsv(OUT/"PI4K_QC_by_site.tsv",site_header,site_rows)

    # Source integrity and machine-readable summary.
    source_digest=digest(SRC)
    summary={
        "status":"OFFLINE_QC_COMPLETE",
        "generated_utc":datetime.now(timezone.utc).isoformat(),
        "source_file":str(SRC),
        "source_sha256":source_digest,
        "gene":"PF3D7_0509800",
        "region":"Pf3D7_05_v3:409403-415758",
        "source":"MalariaGEN Pf8 public SNP-only Zarr; local verified NPZ",
        "samples":len(ids),
        "states":STATES,
        "sites_all":len(positions),
        "sites_pf8_pass":int(passmask.sum()),
        "pf8_pass_sites_with_two_or_more_allele_codes_observed":pass_variants,
        "pf8_pass_sites_with_at_least_one_alt_slot":pass_with_alt,
        "overall_all_site_qc":metrics_for_calls(gt),
        "overall_pass_site_qc":metrics_for_calls(gt[passmask]),
        "qc_tables":["PI4K_QC_by_sample.tsv","PI4K_QC_by_state.tsv","PI4K_QC_by_site.tsv"],
        "limitations":[
          "No clinical interpretation or geographic classification performed.",
          "Pf8 PASS is a site-level flag, not a per-sample call quality guarantee.",
          "Two genotype slots are source encoding, not evidence of biological diploidy.",
          "Discordant two-slot calls should not be equated to a definitive mixed-infection diagnosis.",
          "Alleles in two slots must not be counted as independent haploid parasite observations.",
          "Variable-in-cohort means >=2 distinct allele codes are observed across the encoded slots; does not establish independent between-isolate variation.",
          "No original DP/GQ/QUAL values are included in these derived genotypes.",
          "SNP-only dataset excludes indels; sites are globally discovered, not necessarily polymorphic in the study cohort.",
          "A 95% full-call column is a descriptive scenario only; no threshold applied to exclude sites."
        ]
    }
    (OUT/"PI4K_QC_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print("SUCCESS: Offline QC audit completed; no genomic downloads and no genotype filtering.")
    print("Samples: 282; PI4K SNPs: 2471; Pf8 PASS positions: 1824")
    print("PASS positions with >=2 observed allele codes:",pass_variants)
    print("PASS positions with an ALT-coded genotype slot:",pass_with_alt)
    print("PASS-site missing slot fraction:",summary["overall_pass_site_qc"]["fraction_missing_slots"])
    print("PASS-site discordant two-slot fraction among fully observed calls:",summary["overall_pass_site_qc"]["fraction_discordant_among_fully_observed"])
    print("Reports:",OUT)


if __name__=="__main__":
    run()
