#!/usr/bin/env python3
"""Lightweight public-repo preflight. Not a substitute for human review."""
from pathlib import Path
import re, sys
repo = Path(__file__).resolve().parent.parent
excluded_dirs = {'results','variants','metadata','data','private','manuscript','thesis','reference','logs','output','figures','.git','__pycache__'}
forbidden_exts = {'.tsv','.csv','.vcf','.bcf','.npz','.rds','.tif','.tiff','.png','.jpg','.jpeg','.pdf','.zip','.bam','.bai','.cram','.fastq','.fasta','.gff','.parquet'}
scan = [p for p in repo.rglob('*') if p.is_file() and '.git' not in p.parts]
issues=[]
for f in scan:
    relative=f.relative_to(repo)
    if any(part.lower() in excluded_dirs for part in relative.parts[:-1]): issues.append(('BLOCKED_DIRECTORY',str(relative)))
    if f.suffix.lower() in forbidden_exts: issues.append(('BLOCKED_EXTENSION',str(relative)))
    if f.suffix.lower() in {'.r','.py','.md','.sh','.txt'}:
        s=f.read_text(errors='replace')
        for needle,label in [
            (r'ghp_[A-Za-z0-9]{20,}', 'POSSIBLE_GITHUB_TOKEN'),
            (r'github_pat_[A-Za-z0-9_]{20,}', 'POSSIBLE_GITHUB_TOKEN'),
            (r'-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----','POSSIBLE_PRIVATE_KEY'),
        ]:
            if re.search(needle,s): issues.append((label,str(relative)))
print('Repository:',repo)
print('Files scanned:',len(scan))
print('Python sources:',len(list((repo/'python').glob('*.py'))))
print('R sources:',len(list((repo/'R').glob('*.R'))))
for flag,p in issues: print(flag,p)
if issues:
    print('FAIL: Review/remove forbidden files or secrets before publishing.')
    sys.exit(1)
print('PASS: No blocked data-like files or basic secret patterns detected.')
print('NOTE: Scripts are exact project sources and may disclose study identity, counts, coordinates or findings; inspect manually.')
