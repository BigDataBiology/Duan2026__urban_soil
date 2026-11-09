import requests
from fasta import fasta_iter
def get_sequences(sample,outfile):
    from glob import glob
    from fasta import fasta_iter
    fs = glob('/home1/duanyq/soil/04_binning/cpsnj04/barrnap/*.faa')
    with open(outfile,'wt') as out:
        for f in fs:
            for h,seq in fasta_iter(f):
                mag = f.split('/')[-1].split('.')[0]
                if h.startswith('16S_rRNA::'):
                    h = f'{sample}_{mag}_{h}'
                    out.write(f'>{h}\n{seq}\n')

#sample = 'cpsnj04'
#outfile = f'~/soil/04_binning/barrnap/{sample}_16S.fna'
#get_sequences(sample,outfile)

def filter_sgb(infile1,infile2,outfile):
    from fasta import fasta_iter
    hq = set()
    with open(infile1,'rt') as f:
        for line in f:
            hq.add(line.strip())
    with open(outfile,'wt') as out:
        for h,seq in fasta_iter(infile2):
            mag = h.split('#')[0]
            if mag in hq:
                out.write(f'>{h}\n{seq}\n')
       
infile1 = 'sgb_hq.txt'
infile2 = 'mag_16S.fna'
outfile = 'mag_hq_16S.fna'

import requests
from fasta import fasta_iter

def lookup_1seq(h, seq):
    from time import sleep
    sleep(0.2)  # to avoid hitting the server too hard
    URL = 'https://microbeatlas.org/ma_gw/index.py'
    payload = {
        "action": "mapseq",
        "seq": f'>{h}\\n{seq}',
        "seq_original": f'>{h}\n{seq}',
        "command": f'@map:mapseq(">{h}\\n{seq}", "-1")',
    }
    payload

    r = requests.post(
            URL,
            data=payload,
            headers = {
                'User-Agent': 'Python/requests (custom script. Hi from Brisbane, this is Luis. How is everyone doing in Zurich?)',
            }
        )
    if r.status_code != 200:
        print(r.status_code)
        raise OSError(f"Error {r.status_code} from Microbe Atlas server")
    return r.status_code, r.json()

def reorganize_results(results):
    """
    Reorganize the results from the Microbe Atlas server into a DataFrame.
    Parameters:
    results (list): List of tuples containing sequence, headers, and response from the server.

    Returns:
    pd.DataFrame: DataFrame with columns ['MAG', 'Seq', 'OTU'].
    """
    import pandas as pd
    reorg = []
    for seq, bin, header,(st, r) in results:
        assert st == 200, f"Error {st} from Microbe Atlas server"
        otus = r['otutable']['otus']
        if len(otus) == 0:
            otu = 'NA'
        elif len(otus) == 1:
            [otu] = otus
        else:
            raise ValueError(f"Unexpected number of OTUs: {len(otus)} for sequence {seq}")
        reorg.append((bin,header,seq, otu))
    data = pd.DataFrame(reorg, columns=['MAG', 'header', 'Seq', 'OTU'])
    data.sort_values(by='MAG', inplace=True)
    return data.reset_index(drop=True)


def save_results(final):
    """
    Save the final DataFrame to a CSV file.
    """
    import pandas as pd
    final.to_csv('microbe-atlas.csv', index=False)

results = []
for h,seq in fasta_iter('mag_hq_16S.fna'):
    binid = h.split('#')[0]
    contig = h.split(':')[2]
    results.append((seq, binid, h, lookup_1seq(h, seq)))
    final = reorganize_results(results)
    save_results(final)

def cal(infile1):
    otu90 = 0
    otu96 = 0
    otu97 = 0
    otu98 = 0
    otu99 = 0
    mag_dict = set()
    with open(infile1,'rt') as f:
        for line in f:
            if line.startswith('MAG'):
                continue
            else:
                mag,header,seq,otu = line.strip().split(',')
                otulist = otu.split(';')
                if mag not in mag_dict:
                    if otulist[-1].startswith('90'):
                        otu90 += 1
                    if otulist[-1].startswith('96'):
                        otu96 += 1  
                    if otulist[-1].startswith('97'):
                        otu97 += 1  
                    if otulist[-1].startswith('98'):
                        otu98 += 1
                    if otulist[-1].startswith('99'):
                        otu99 += 1
                mag_dict.add(mag)
    print(f'OTU90: {otu90}')
    print(f'OTU96: {otu96}') 
    print(f'OTU97: {otu97}')    
    print(f'OTU98: {otu98}')
    print(f'OTU99: {otu99}')        

infile1 = 'microbe-atlas.csv'
cal(infile1)

def map_tax(infile1,infile2,outfile):
    gtdbdict = {}
    with open(infile1,'rt') as f:
        for line in f:
            linelist = line.strip().split(',')
            gtdbdict[linelist[7]] = line.strip()
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                if line.startswith('MAG'):
                    continue
                else:
                    linelist = line.strip().split(',')
                    out.write(f'{gtdbdict[linelist[0]]},{linelist[1]},{linelist[2]},{linelist[3]}\n')


infile1 = 'sgb_gtdb.csv'
infile2 = r'D:\soil\barrnap\microbe-atlas.csv'
outfile = r'sgb_16s_gtdb.csv'
map_tax(infile1,infile2,outfile)