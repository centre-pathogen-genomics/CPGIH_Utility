#!/usr/bin/env python

# Build the per-species prokaryotic genome size table used by ont_genomesqc.sh
#
# USAGE: build_ncbi_genome_size_db.py --workdir /local/dir [--output data/prok_species_genome_size_stats.tsv]
#
# Steps
#   1. Pull GenBank assembly metadata for Bacteria (taxid 2) and Archaea (taxid 2157) with the
#      NCBI Datasets CLI (--assembly-source genbank --exclude-atypical; MAGs included). A taxon pull
#      is skipped if already completed, so the script can be re-run after a failure.
#   2. Map accession -> species_taxid from NCBI assembly_summary_genbank.txt (column 7)
#   3. Download the current NCBI taxdump (for taxonkit species names)
#   4. Summarise ungapped genome length per species_taxid for all assemblies (_all) and
#      'Complete Genome' assemblies only (_complete)
#
# Quality filtering
#   --exclude-atypical removes assemblies NCBI flags as too small/large, contaminated, partial etc.
#   CheckM filtering is controlled by --checkm_filter:
#     none          (default) no CheckM filter - CheckM values are missing for ~2/3 of GenBank
#                   assemblies and the CheckM filter falsely fails reduced-genome taxa (e.g. Mycoplasma)
#     exclude_fail  drop assemblies with CheckM values failing the thresholds; keep those with no CheckM
#     pass_only     keep only assemblies with CheckM values passing the thresholds
#
# Outputs
#   --output                 compact table for ont_genomesqc.sh
#                            (species, species_taxid, n_all, median_all, n_complete, median_complete)
#   <workdir>/prok_species_genome_size_stats_full.tsv
#                            full table (mean/median/SD, CheckM status counts, name source)
#   <workdir>/build.log      progress log
#
# Requires: ncbi-datasets-cli (datasets, dataformat), taxonkit, curl, pandas
# Large intermediate files (~2 GB) are written to --workdir - keep this off OneDrive

import argparse
import io
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd

TAXA = {'bacteria': 2, 'archaea': 2157}
FIELDS = ['accession', 'organism-name', 'organism-tax-id', 'assminfo-level', 'assmstats-total-ungapped-len',
          'checkm-completeness', 'checkm-contamination']
META_COLS = ['accession', 'organism_name', 'organism_taxid', 'assembly_level', 'ungapped_len',
             'checkm_completeness', 'checkm_contamination']
ASM_SUMMARY_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'
TAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
DEFAULT_OUTPUT = Path(__file__).resolve().parent.parent / 'data' / 'prok_species_genome_size_stats.tsv'

log = logging.getLogger('genome_size_db')


def run_pipe(cmd, desc):
    log.info(f'{desc}: {cmd}')
    res = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], capture_output=True, text=True)
    if res.returncode != 0:
        log.error(f'{desc} failed (exit {res.returncode})\n{res.stderr.strip()}')
        sys.exit(1)


def check_tools():
    missing = [t for t in ['datasets', 'dataformat', 'taxonkit', 'curl'] if shutil.which(t) is None]
    if missing:
        log.error(f'Required tools not found on PATH: {", ".join(missing)}')
        sys.exit(1)
    for cmd in (['datasets', '--version'], ['taxonkit', 'version']):
        log.info(subprocess.run(cmd, capture_output=True, text=True).stdout.strip())


def pull_metadata(workdir, api_key):
    # one file per taxon; a completed pull is marked with a .done file and not re-pulled
    key = f' --api-key {api_key}' if api_key else ''
    out_files = []
    for name, taxid in TAXA.items():
        out = workdir / f'meta_{name}.tsv'
        done = out.with_suffix('.done')
        out_files.append(out)
        if done.exists():
            log.info(f'{out.name} already complete - skipping')
            continue
        tmp = out.with_suffix('.tmp')
        cmd = (f'datasets summary genome taxon {taxid} --assembly-source genbank --exclude-atypical '
               f'--as-json-lines{key} | '
               f'dataformat tsv genome --fields {",".join(FIELDS)} --elide-header > {tmp}')
        run_pipe(cmd, f'Pulling {name} (taxid {taxid})')
        tmp.rename(out)
        done.touch()
        n = sum(1 for _ in open(out))
        log.info(f'{out.name}: {n} records')
    return out_files


def get_acc2species(workdir):
    # stream the ~2 GB summary and keep only accession + species_taxid
    out = workdir / 'acc2species_taxid.tsv'
    if out.exists() and out.with_suffix('.done').exists():
        log.info(f'{out.name} already complete - skipping')
        return out
    tmp = out.with_suffix('.tmp')
    run_pipe(f'curl -fsSL {ASM_SUMMARY_URL} | grep -v "^#" | cut -f 1,7 > {tmp}',
             'Downloading assembly_summary_genbank.txt (accession, species_taxid)')
    tmp.rename(out)
    out.with_suffix('.done').touch()
    return out


def get_taxdump(workdir):
    tdir = workdir / 'taxdump'
    if (tdir / 'names.dmp').exists() and (tdir / '.done').exists():
        log.info('taxdump already present - skipping')
        return tdir
    tdir.mkdir(exist_ok=True)
    tgz = tdir / 'taxdump.tar.gz'
    log.info(f'Downloading {TAXDUMP_URL}')
    urllib.request.urlretrieve(TAXDUMP_URL, tgz)
    with tarfile.open(tgz) as tf:
        for m in ['names.dmp', 'nodes.dmp', 'merged.dmp', 'delnodes.dmp']:
            tf.extract(m, tdir)
    tgz.unlink()
    (tdir / '.done').touch()
    return tdir


def taxonkit_names(taxids, data_dir):
    res = subprocess.run(['taxonkit', 'lineage', '-L', '-n', '--data-dir', str(data_dir)],
                         input='\n'.join(str(t) for t in taxids) + '\n', capture_output=True, text=True)
    if res.returncode != 0:
        log.error(f'taxonkit failed\n{res.stderr}')
        sys.exit(1)
    names = pd.read_csv(io.StringIO(res.stdout), sep='\t', header=None,
                        names=['species_taxid', 'taxonkit_name'], dtype={'taxonkit_name': str})
    names['species_taxid'] = names['species_taxid'].astype('Int64')
    return names.dropna(subset=['taxonkit_name']).drop_duplicates('species_taxid').set_index('species_taxid')


def summarise(df, suffix):
    g = df.groupby('species_taxid')['ungapped_len']
    out = pd.DataFrame({'n': g.size(), 'mean': g.mean(), 'median': g.median(), 'sd': g.std(ddof=1)})
    for c in ['mean', 'median', 'sd']:
        out[c] = out[c].round(0).astype('Int64')
    out.columns = [f'{c}_{suffix}' for c in out.columns]
    return out


def fallback_name(organism_names):
    # most common first-two-words of organism_name; keep the full name for 'sp.' / uncultured
    # organisms so distinct taxids don't collapse to the same name
    def short(n):
        w = n.split()
        if len(w) >= 2 and (w[1] == 'sp.' or w[0] in ('uncultured', 'Candidatus')):
            return n
        return ' '.join(w[:2])
    return organism_names.map(short).value_counts().index[0]


def build_table(meta_files, acc2sp_file, taxdump, args):
    meta = pd.concat([pd.read_csv(f, sep='\t', header=None, names=META_COLS,
                                  dtype={'accession': str, 'organism_name': str, 'assembly_level': str},
                                  low_memory=False) for f in meta_files], ignore_index=True)
    for c in ['ungapped_len', 'checkm_completeness', 'checkm_contamination']:
        meta[c] = pd.to_numeric(meta[c], errors='coerce')
    n_raw = len(meta)
    meta = meta.drop_duplicates('accession')

    acc2sp = pd.read_csv(acc2sp_file, sep='\t', header=None, names=['accession', 'species_taxid'],
                         dtype={'accession': str})
    acc2sp['species_taxid'] = pd.to_numeric(acc2sp['species_taxid'], errors='coerce').astype('Int64')
    meta = meta.merge(acc2sp.drop_duplicates('accession'), on='accession', how='left')

    has_checkm = meta['checkm_completeness'].notna() & meta['checkm_contamination'].notna()
    checkm_ok = (meta['checkm_completeness'] >= args.min_completeness) & \
                (meta['checkm_contamination'] <= args.max_contamination)
    meta['checkm_status'] = 'missing'
    meta.loc[has_checkm & checkm_ok, 'checkm_status'] = 'pass'
    meta.loc[has_checkm & ~checkm_ok, 'checkm_status'] = 'fail'

    if args.checkm_filter == 'none':
        passing = pd.Series(True, index=meta.index)
    elif args.checkm_filter == 'exclude_fail':
        passing = meta['checkm_status'] != 'fail'
    else:
        passing = meta['checkm_status'] == 'pass'
    missing_sp = meta['species_taxid'].isna()

    log.info('QC counts')
    log.info(f'  rows read:                       {n_raw}')
    log.info(f'  total assemblies (dedup):        {len(meta)}')
    log.info(f'  CheckM pass/fail/missing:        ' + '/'.join(
        str((meta['checkm_status'] == s).sum()) for s in ['pass', 'fail', 'missing']) +
        f' (completeness >= {args.min_completeness}, contamination <= {args.max_contamination})')
    log.info(f'  missing species_taxid:           {missing_sp.sum()}')
    log.info(f'  missing ungapped length:         {meta["ungapped_len"].isna().sum()}')
    log.info(f'  passing CheckM filter ({args.checkm_filter}): {passing.sum()}')

    keep = meta[passing & ~missing_sp & meta['ungapped_len'].notna()].copy()
    log.info(f'  assemblies used:                 {len(keep)}')

    stats = summarise(keep, 'all')
    counts = keep.groupby(['species_taxid', 'checkm_status']).size().unstack(fill_value=0)
    for s in ['pass', 'fail', 'missing']:
        stats[f'n_checkm_{s}'] = counts[s] if s in counts else 0
    stats = stats.join(summarise(keep[keep['assembly_level'] == 'Complete Genome'], 'complete'), how='left')
    stats['n_complete'] = stats['n_complete'].fillna(0).astype(int)

    names = taxonkit_names(stats.index.tolist(), taxdump)
    stats['species'] = names['taxonkit_name'].reindex(stats.index)
    stats['name_source'] = 'taxonkit'
    no_name = stats['species'].isna()
    if no_name.any():
        fb = keep[keep['species_taxid'].isin(stats.index[no_name])].groupby('species_taxid')['organism_name'] \
            .agg(fallback_name)
        stats.loc[no_name, 'species'] = fb.reindex(stats.index[no_name])
        stats.loc[no_name, 'name_source'] = 'organism_name_fallback'
    log.info(f'  species:                         {len(stats)} ({no_name.sum()} named by organism_name fallback)')

    stats = stats.reset_index().sort_values(['species', 'species_taxid'])
    full_cols = ['species', 'species_taxid', 'name_source',
                 'n_all', 'n_checkm_pass', 'n_checkm_fail', 'n_checkm_missing',
                 'mean_all', 'median_all', 'sd_all',
                 'n_complete', 'mean_complete', 'median_complete', 'sd_complete']
    return stats[full_cols]


def main():
    parser = argparse.ArgumentParser(description='Build per-species prokaryotic genome size table from NCBI GenBank metadata')
    parser.add_argument('--workdir', required=True, help='Directory for intermediate files (~2 GB; keep off OneDrive)')
    parser.add_argument('--output', default=str(DEFAULT_OUTPUT), help=f'Compact table for ont_genomesqc.sh (default {DEFAULT_OUTPUT})')
    parser.add_argument('--checkm_filter', choices=['none', 'exclude_fail', 'pass_only'], default='none',
                        help='CheckM filtering mode (default none)')
    parser.add_argument('--min_completeness', type=float, default=90, help='CheckM completeness threshold (default 90)')
    parser.add_argument('--max_contamination', type=float, default=5, help='CheckM contamination threshold (default 5)')
    parser.add_argument('--api_key', default=os.environ.get('NCBI_API_KEY'), help='NCBI API key (default $NCBI_API_KEY)')
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        handlers=[logging.FileHandler(workdir / 'build.log'), logging.StreamHandler()])
    log.info(f'Starting build; workdir={workdir}; checkm_filter={args.checkm_filter}')

    check_tools()
    meta_files = pull_metadata(workdir, args.api_key)
    acc2sp = get_acc2species(workdir)
    taxdump = get_taxdump(workdir)
    stats = build_table(meta_files, acc2sp, taxdump, args)

    full_out = workdir / 'prok_species_genome_size_stats_full.tsv'
    stats.to_csv(full_out, sep='\t', index=False, na_rep='NA')
    log.info(f'Full table written to {full_out}')

    # compact table: one row per species name (Kraken lookups are by name);
    # where names collide keep the taxid with the most assemblies
    compact = stats.sort_values('n_all', ascending=False).drop_duplicates('species')
    n_dup = len(stats) - len(compact)
    if n_dup:
        log.info(f'  {n_dup} rows dropped from compact table due to duplicated species names (kept largest n_all)')
    compact = compact.sort_values('species')[['species', 'species_taxid', 'n_all', 'median_all', 'n_complete', 'median_complete']]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    compact.to_csv(args.output, sep='\t', index=False, na_rep='NA')
    log.info(f'Compact table ({len(compact)} species) written to {args.output}')


if __name__ == '__main__':
    main()
