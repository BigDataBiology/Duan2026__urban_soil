import polars as pl
from pathlib import Path
from collections import Counter


def create_emapper_summary(base):
    base = Path(base)
    tag = base.parts[-1]
    gff = pl.read_csv(base / f'{tag}.emapper.genepred.gff',
                    has_header=False,
                    separator='\t',
                    comment_prefix='#',
                    new_columns=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    gff = gff.rename({'seqid': 'contig'})
    seen = Counter()
    full_geneid = []
    for geneid in gff['contig'].to_list():
        seen[geneid] += 1
        full_geneid.append(f'{geneid}_{seen[geneid]}')
    gff = gff.with_columns(seqid=pl.Series(full_geneid))

    emapper = pl.read_csv(base / f'{tag}.emapper.annotations', has_header=True, separator='\t', comment_prefix='##')
    emapper = emapper.rename({'#query': 'query'})
    assert not set(emapper['query'].to_list()) - set(gff['seqid'].to_list())
    EMAPPER_COLS = ['query', 'COG_category', 'Preferred_name', 'KEGG_ko', 'KEGG_Module']

    data = gff['seqid', 'contig', 'start', 'end', 'strand'].join(
            emapper[EMAPPER_COLS], left_on='seqid', right_on='query')
    oname = base / f'{tag}.emapper_summary.tsv'
    data.write_csv(oname, separator='\t')
    return oname

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Create summary file for emapper annotations')
    parser.add_argument('base', help='Base directory containing emapper output files')
    args = parser.parse_args()
    outname = create_emapper_summary(args.base)
    print(f'Created summary file: {outname}')

