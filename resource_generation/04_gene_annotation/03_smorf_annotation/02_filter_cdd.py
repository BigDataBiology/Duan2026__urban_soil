def add_length(infile1,infile2,outfile):
    cdd_dict = {}
    with open(infile1,'rt') as f1:
        for line in f1:
            linelist = line.strip().split('\t')
            cdd_dict[linelist[0]] = linelist[4]

    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f2:
            for line in f2:
                line = line.strip()
                linelist = line.split('\t')
                pssm = linelist[1].split('|')[2]
                if cdd_dict[pssm] != '0':
                    out.write(f'{line}\t{cdd_dict[pssm]}\n')

def filter_cov(infile,outfile):
    import pandas as pd
    result = pd.read_csv(infile,sep='\t',header=None,names=['smorf','cdd','query_length','score','align_length','identity','evalue','target_length'])
    result['tcov'] = result['align_length']/result['target_length']
    result = result[result['tcov'] >0.8]
    result.to_csv(outfile,sep='\t',index=None)

def filter_family(infile3,outfile2,outfile3):
    cluster = set()
    with open(infile3,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            cluster.add(linelist[2])
    with open(outfile3,'wt') as out:
        with open(outfile2,'rt') as f:
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in cluster:
                    out.write(line)

def cal(infile,outfile):
    import pandas as pd
    result = pd.read_csv(infile,sep='\t',header=None,names=['smorf','cdd','query_length','score','align_length','identity','evalue','target_length','tcov'])
    df = result.smorf.value_counts()
    bg_count = df.shape[0]
    bg_multi_count = df[df>1].shape[0]
    bg_perc = bg_count/2206901
    bg_multi_perc = bg_multi_count/2206901

    bg_cdd_count = result.cdd.value_counts()
    bg_cdd_count = bg_cdd_count.to_frame(name='bg_count').reset_index().rename(columns={'index':'cdd'})

    bg_cdd_count['bg_percentage'] = bg_cdd_count['bg_count']/bg_cdd_count['bg_count'].sum()
    bg_cdd_count.to_csv(outfile,sep='\t',index=None)

def map_info(infile,outfile):
    import pandas as pd
    result = pd.read_csv(infile,sep='\t',header=None,names=['PSSM_Id','GMSC_count','GMSC_fraction'],skiprows=1)
    result['PSSM_Id'] = result['PSSM_Id'].apply(lambda x:x.split('|')[2])
    result['PSSM_Id'] = result['PSSM_Id'].astype('int')
    cdd = pd.read_csv(r'D:\soil\smorf\cddid_all.tbl',sep='\t',header=None,names=['PSSM_Id','accession','short_name','description','PSSM_Length'])
    result = result.merge(cdd,'left',on='PSSM_Id')
    result.to_csv(outfile,sep='\t',index=None)

def pfam2clan(infile1,infile2,outfile):
    pfam_dict = {}
    with open(infile1,'rt') as f1:
        for line in f1:
            linelist = line.strip().split('\t')
            pfam = linelist[0].replace('PF','')
            pfam_dict[pfam] = f'{linelist[1]}\t{linelist[2]}\t{linelist[4]}' 
    with open(outfile,'wt') as out:
        out.write(f'PSSM_Id\tGMSC_count\tGMSC_fraction\taccession\tshort_name\tdescription\tPSSM_Length\tclan_id\tclan\tshort_description\n')
        with open(infile2,'rt') as f2:
            for line in f2:
                if line.startswith('PSSM_Id'):
                    continue
                else:
                    line = line.replace('\n','')
                    linelist = line.split('\t')
                    pf = linelist[3].replace('pfam','')
                    if pf in pfam_dict.keys():
                        out.write(f'{line}\t{pfam_dict[pf]}\n')
                    else:
                        out.write(f'{line}\n')

def format(infile,outfile):
    pfam_dict = {}
    with open(infile,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('#=GF ID'):
                id = line.strip().replace('#=GF ID   ','')
                pfam_dict[id] = ['',[]]
            if line.startswith('#=GF DE'):
                de = line.strip().replace('#=GF DE   ','')
                pfam_dict[id][0] = de
            if line.startswith('#=GF MB'):
                pfam = line.strip().replace('#=GF MB   PF','').replace(';','')
                pfam_dict[id][1].append(pfam)
    with open(outfile,'wt') as out:
        for key,value in pfam_dict.items():
            for item in value[1]:
                out.write(f'{item}\t{key}\t{value[0]}\n')

def clanc(infile1,infile2,outfile):
    import pandas as pd
    df = pd.read_csv(infile1,sep='\t')
    df = df[df['accession'].str.startswith('pfam')]
    clan_c = pd.read_csv(infile2,sep='\t',header=None,dtype=str,names=['accession','clan','clan_description'])
    df['accession'] = df['accession'].str.replace('pfam','')
    df=clan_c.merge(df,'right',on=['accession','clan'])
    df.to_csv(outfile,sep='\t',index=None)

# group ribosomal unknown function protein
def modify(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('accession'):
                    out.write(f'{line}\tgroup\n')
                else:
                    linelist = line.split('\t')
                    if 'ribosomal' in linelist[-1].casefold():
                        out.write(f'{line}\tRibosomal protein\n')
                    else:
                        if linelist[1]!='':
                            if linelist[2]!='':
                                out.write(f'{line}\t{linelist[2]}\n')
                            else:
                                out.write(f'{line}\t{linelist[1]}\n')
                        else:
                            if 'uncharacterised' in linelist[-1].casefold() or 'unknown' in linelist[-1].casefold():
                                out.write(f'{line}\tDomain of unknown function\n')
                            else:
                                group = linelist[-1].split(',')[0]
                                out.write(f'{line}\t{group}\n')          

infile1 = 'cddid_all.tbl'
infile2 = 'all.mapped.smorfs_dedup.tsv'
infile3 = 'cluster_count.tsv'
infile4 = 'Pfam-A.clans.tsv'
infile5 = 'pfam-c.txt'
infile6 = 'pfam-c.tsv'
outfile1 = 'all.mapped.smorfs_dedup_cdd.tsv'
outfile2 = 'all.mapped.smorfs_dedup_cdd_tcov.tsv'
outfile3 = 'all.mapped.smorfs_dedup_cdd_tcov_family.tsv'
outfile4 = 'all.mapped.smorfs_dedup_cdd_tcov_family_bg.tsv'
outfile5 = 'all.mapped.smorfs_dedup_cdd_tcov_family_bg_info.tsv'
outfile6 = 'all.mapped.smorfs_dedup_cdd_tcov_family_bg_info_A.tsv'
outfile7 = 'all.mapped.smorfs_dedup_cdd_tcov_family_bg_info_A_C.tsv'
outfile8 = 'all.mapped.smorfs_dedup_cdd_tcov_family_bg_info_A_C_group.tsv'

add_length(infile1,infile2,outfile1)
filter_cov(outfile1,outfile2)
filter_family(infile3,outfile2,outfile3)
cal(outfile3,outfile4)
map_info(outfile4,outfile5)
pfam2clan(infile4,outfile5,outfile6)
format(infile5,infile6)
clanc(outfile6,infile6,outfile7)
modify(outfile7,outfile8)