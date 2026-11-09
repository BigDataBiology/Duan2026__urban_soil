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
for h,seq in fasta_iter('test.fna'):
    binid = h.split('#')[0]
    contig = h.split(':')[2]
    results.append((seq, binid, h, lookup_1seq(h, seq)))
    final = reorganize_results(results)
    save_results(final)