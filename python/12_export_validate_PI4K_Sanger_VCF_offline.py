#!/usr/bin/env python3
"""OFFLINE Step 17: repair/validate PI4K derived VCFs using existing verified SNP NPZ.
No network, no genomic data download, no changes to NPZ source files.
"""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import numpy as np

STATES = {"Oyo": 71, "Ogun": 38, "Lagos": 129, "Osun": 15, "Ondo": 29}
CHROM = 'Pf3D7_05_v3'
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'variants' / 'Southwest5' / 'Sanger_SNP_complete'
META = ROOT / 'results' / 'Southwest5'
FILE = OUT / 'PI4K_Southwest5_282_Sanger_SNP_all_sites.npz'
VCF = OUT / 'PI4K_Southwest5_282_Sanger_SNP_derived.vcf.gz'


def run(args, **kwargs):
    return subprocess.run(args, check=True, **kwargs)


def allele_str(x):
    if isinstance(x, (bytes, np.bytes_)):
        return bytes(x).decode('utf-8').strip()
    s = str(x).strip()
    return '' if s.lower() in ('none', 'nan', 'null', '.') else s


def valid_alleles(row):
    r = [allele_str(x) for x in row]
    while len(r) > 1 and not r[-1]:
        r.pop()
    if not r or not all(r):
        raise ValueError('Invalid REF or ALT list in source NPZ')
    return r


def check_vcf(path, ids, positions, passes, alleles, gt):
    if not path.exists():
        return False, 'file does not exist'
    try:
        if not Path(str(path) + '.csi').exists():
            run(['bcftools', 'index', '-c', '-f', str(path)], stdout=subprocess.DEVNULL)
        run(['bcftools', 'view', '-h', str(path)], stdout=subprocess.DEVNULL)
        with subprocess.Popen(['bcftools', 'query', '-l', str(path)], stdout=subprocess.PIPE, text=True) as p:
            found = [s.strip() for s in p.stdout if s.strip()]
            status = p.wait()
            if status != 0 or found != ids:
                raise ValueError('VCF sample ID list or ordering differs from NPZ')
        with subprocess.Popen(['bcftools', 'view', '-H', str(path)], stdout=subprocess.PIPE, text=True) as p:
            n = 0
            for line in p.stdout:
                if n >= len(positions):
                    raise ValueError('Too many variant records')
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9 + len(ids):
                    raise ValueError(f'Variant {n}: wrong number of columns')
                want_al = valid_alleles(alleles[n])
                if (fields[0] != CHROM or int(fields[1]) != int(positions[n]) or
                    fields[3] != want_al[0] or fields[4] != (','.join(want_al[1:]) if len(want_al)>1 else '.') or
                    fields[6] != ('PASS' if bool(passes[n]) else 'PF8_NONPASS')):
                    raise ValueError(f'Variant {n}: locus/alleles/filter differs from NPZ')
                for i, samplefield in enumerate(fields[9:]):
                    wanted = '/'.join('.' if int(a)<0 else str(int(a)) for a in gt[n,i,:])
                    if samplefield.split(':')[0].replace('|','/') != wanted:
                        raise ValueError(f'Variant {n}: genotype mismatch for {ids[i]}')
                n += 1
            status = p.wait()
            if status != 0:
                raise ValueError('bcftools could not finish reading all variant records')
            if n != len(positions):
                raise ValueError(f'VCF has {n} records; expected {len(positions)}')
        return True, f'{len(positions)} variants, {len(ids)} samples, genotype-by-genotype verified'
    except Exception as e:
        return False, f'{type(e).__name__}: {e}'


def write_vcf(path, ids, positions, passes, alleles, gt):
    print('Rebuilding combined VCF locally from NPZ; no network access', flush=True)
    with tempfile.TemporaryDirectory(prefix='pi4k_vcf_', dir=str(OUT)) as temp:
        plain = Path(temp)/'build.vcf'
        gz = Path(temp)/'build.vcf.gz'
        with plain.open('w', encoding='utf-8', newline='\n') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##source=MalariaGEN_Pf8_Sanger_public_snp_only_derived\n')
            f.write('##reference=PlasmoDB-54-Pfalciparum3D7-Genome\n')
            f.write(f'##contig=<ID={CHROM}>\n')
            f.write('##FILTER=<ID=PF8_NONPASS,Description="Pf8 variant_filter_pass is false; original filter codes not retained">\n')
            f.write('##INFO=<ID=PF8_FILTER_PASS,Number=1,Type=Integer,Description="Pf8 Zarr variant_filter_pass Boolean">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Original Pf8 two-slot call_genotype encoding; biological ploidy is not inferred">\n')
            f.write('##pi4k_note=SNP-only derived VCF without original QUAL/GQ/DP, with two-slot genotype encoding\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n')
            for n,p in enumerate(positions):
                al = valid_alleles(alleles[n])
                if np.min(gt[n],initial=0) < -1 or np.max(gt[n],initial=0) >= len(al):
                    raise ValueError(f'Out of range allele codes at site {n}')
                calls = ['/'.join('.' if int(a)<0 else str(int(a)) for a in pair) for pair in gt[n]]
                f.write(f'{CHROM}\t{int(p)}\t.\t{al[0]}\t{(",".join(al[1:]) if len(al)>1 else ".")}\t.\t{("PASS" if bool(passes[n]) else "PF8_NONPASS")}\tPF8_FILTER_PASS={int(passes[n])}\tGT\t'+'\t'.join(calls)+'\n')
        with gz.open('wb') as target:
            run(['bgzip', '-c', str(plain)], stdout=target)
        run(['bcftools', 'index', '-c', '-f', str(gz)])
        ok, reason = check_vcf(gz, ids, positions, passes, alleles, gt)
        if not ok:
            raise RuntimeError(f'Newly generated VCF failed verification: {reason}')
        os.replace(gz, path)
        os.replace(Path(str(gz)+'.csi'), Path(str(path)+'.csi'))
    return reason


def main():
    for command in ['bcftools', 'bgzip']:
        if shutil.which(command) is None:
            raise RuntimeError(f'{command} not installed or not on PATH')
    if not FILE.is_file():
        raise FileNotFoundError(FILE)
    with np.load(FILE, allow_pickle=False) as d:
        positions = np.asarray(d['positions'], dtype=np.int32)
        ids = d['sample_ids'].astype(str).tolist()
        states = d['states'].astype(str).tolist()
        genotypes = np.asarray(d['genotypes'], dtype=np.int8)
        passes = np.asarray(d['filter_pass'], dtype=bool)
        alleles = d['alleles'].astype(str)
    if (len(ids) != 282 or len(set(ids)) != 282 or len(positions) != 2471 or
        int(passes.sum()) != 1824 or genotypes.shape != (2471,282,2) or
        alleles.shape[0] != 2471 or not np.all(np.diff(positions)>0) or
        int(positions[0]) < 409403 or int(positions[-1]) > 415758):
        raise RuntimeError('Source NPZ failed expected shape, ID, locus or PASS checks')
    for state,n in STATES.items():
        index = [i for i,s in enumerate(states) if s == state]
        if len(index) != n:
            raise RuntimeError(f'Unexpected number of samples from {state}')
        idfile = META/'clean_sample_lists'/f'{state}_sample_IDs.txt'
        if not idfile.exists():
            idfile = META/'sample_lists'/f'{state}_sample_IDs.txt'
        want = [s.strip() for s in idfile.read_text(encoding='utf-8-sig').splitlines() if s.strip()]
        if want != [ids[i] for i in index]:
            raise RuntimeError(f'NPZ sample order does not match verified {state} manifest')
    print('Verified: original NPZ has 282 unique samples and 2471 sites (1824 PASS).',flush=True)
    ok, reason = check_vcf(VCF,ids,positions,passes,alleles,genotypes)
    if ok:
        print(f'Existing combined VCF is intact: {reason}',flush=True)
    else:
        print(f'Existing combined VCF not usable: {reason}',flush=True)
        reason = write_vcf(VCF,ids,positions,passes,alleles,genotypes)
        print(f'Combined VCF rebuilt and verified: {reason}',flush=True)
    report=[]
    # Use bcftools solely on the small local VCF; never reopen the source Zarr / VCF URL.
    for state,n in STATES.items():
        out=OUT/f'PI4K_{state}_Sanger_SNP_derived.vcf.gz'
        if out.exists():
            state_ix=[i for i,s in enumerate(states) if s==state]
            ok,why=check_vcf(out,[ids[i] for i in state_ix],positions,passes,alleles,genotypes[:,state_ix,:])
            if ok:
                print(f'{state}: existing VCF validated ({n} samples)',flush=True)
                report.append((state,n,2471,'VALIDATED_EXISTING',out.name))
                continue
            else:
                print(f'{state}: replacing invalid existing file: {why}',flush=True)
        # all IDs are unique and ordered exactly as combined VCF, so selecting by ID cannot reorder unexpectedly
        with tempfile.TemporaryDirectory(prefix=f'pi4k_{state}_',dir=str(OUT)) as tdir:
            tmp=Path(tdir)/'subset.vcf.gz'
            idx=[i for i,s in enumerate(states) if s == state]
            idfile=Path(tdir)/'samples.txt'
            idfile.write_text('\n'.join(ids[i] for i in idx)+'\n',encoding='utf-8')
            run(['bcftools','view','-S',str(idfile),'-Oz','-o',str(tmp),str(VCF)])
            run(['bcftools','index','-c','-f',str(tmp)])
            ok,why=check_vcf(tmp,[ids[i] for i in idx],positions,passes,alleles,genotypes[:,idx,:])
            if not ok:
                raise RuntimeError(f'{state} derived VCF validation failed: {why}')
            os.replace(tmp,out)
            os.replace(Path(str(tmp)+'.csi'),Path(str(out)+'.csi'))
        print(f'{state}: exported {n} samples, all 2471 sites; verified against NPZ',flush=True)
        report.append((state,n,2471,'EXPORTED_VALIDATED',out.name))
    # Filter by original Pf8 site-level PASS flag; 1824 sites retained.
    allpass=OUT/'PI4K_Southwest5_282_Sanger_SNP_PASS_derived.vcf.gz'
    with tempfile.TemporaryDirectory(prefix='pi4k_pass_',dir=str(OUT)) as tdir:
        tmp=Path(tdir)/'pass.vcf.gz'
        run(['bcftools','view','-f','PASS','-Oz','-o',str(tmp),str(VCF)])
        run(['bcftools','index','-c','-f',str(tmp)])
        ok,why=check_vcf(tmp,ids,positions[passes],passes[passes],alleles[passes],genotypes[passes,:,:])
        if not ok:
            raise RuntimeError(f'PASS-only combined file validation failed: {why}')
        os.replace(tmp,allpass)
        os.replace(Path(str(tmp)+'.csi'),Path(str(allpass)+'.csi'))
    print('PASS-only combined VCF exported: 1824 sites, 282 samples',flush=True)
    with (OUT/'PI4K_Sanger_SNP_VCF_export_validation.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['State','Samples','Sites','Result','Filename']);w.writerows(report)
    payload={'status':'VCF_EXPORT_VALIDATED_OFFLINE', 'source':'saved Pf8 SNP-only NPZ','sample_total':282,'states':STATES,'sites_all':2471,'sites_Pf8_PASS':1824,'files':[VCF.name,allpass.name]+[x[-1] for x in report], 'caveats':['Derived genotype-only VCF; not the original Pf8 VCF','Original GQ, DP, QUAL and filter labels not recovered','Two-slot genotype encoding does not establish biological diploidy','No indels in source SNP-only dataset','All 2471 sites are globally ascertained, not necessarily segregating in these samples']}
    (OUT/'PI4K_Sanger_SNP_VCF_export_validation.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print('SUCCESS: 5 state VCFs + combined all-sites VCF + combined PASS VCF, fully cross-checked against saved NPZ.',flush=True)
    print('Saved validation report:',OUT/'PI4K_Sanger_SNP_VCF_export_validation.json',flush=True)

if __name__=='__main__':
    try:
        main()
    except Exception as e:
        print('ERROR:',type(e).__name__,str(e), file=sys.stderr)
        raise
