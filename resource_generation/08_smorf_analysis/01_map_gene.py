def map_gene(infile1,infile2,outfile):
    from fasta import fasta_iter
    contig_dict = {}
    for h,seq in fasta_iter(infile1):
        contigpart = h.replace('mag_','').split('_')
        contig = f'{contigpart[0]}_{contigpart[1]}_{contigpart[2]}'
        if contig not in contig_dict.keys():
            contig_dict[contig] = []
        info = h.replace('mag_','').split(' # ID')[0]
        contig_dict[contig].append(info)

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                copies = line.strip().split('\t')[4].split(',')
                for copy in copies:
                    contigpart = copy.replace('mag_','').split('_')
                    contig = f'{contigpart[0]}_{contigpart[1]}_{contigpart[2]}'
                    infopart = copy.split(' # ')
                    start = int(infopart[1])
                    end = int(infopart[2])
                    strand = int(infopart[3])
                    if strand not in (1, -1):
                        raise ValueError(f'Invalid smORF strand: {copy}')

                    genes = contig_dict[contig]
                    upstream = []
                    downstream = []
                    for gene in genes:
                        g_infopart = gene.split(' # ')
                        g_start = int(g_infopart[1])
                        g_end = int(g_infopart[2])
                        if g_end < start:
                            if strand == 1:
                                upstream.append(gene)
                            else:
                                downstream.append(gene)
                        elif g_start > end:
                            if strand == 1:
                                downstream.append(gene)
                            else:
                                upstream.append(gene)
                    upstream.sort(key=lambda g: (
                        start - int(g.split(' # ')[2]) if strand == 1
                        else int(g.split(' # ')[1]) - end, g))
                    downstream.sort(key=lambda g: (
                        int(g.split(' # ')[1]) - end if strand == 1
                        else start - int(g.split(' # ')[2]), g))
                    name = line.strip().split('\t')[0]
                    upstreams = ','.join(upstream[:10])
                    downstreams = ','.join(downstream[:10])
                    if upstreams == '':
                        upstreams='NA'
                    if downstreams == '':
                        downstreams='NA'
                    out.write(f'{name}\t{upstreams}\t{downstreams}\n')

infile1 = 'all_genes.faa.gz'
infile2 = 'smorf_mapped_contig.tsv'
outfile = 'smorf_genes_10.tsv'
map_gene(infile1,infile2,outfile)
