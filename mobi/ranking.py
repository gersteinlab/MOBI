import os
import pandas as pd

def get_rank_file(infile, outfile, rank):
    df = pd.read_csv(infile, sep="\t", header=None)
        
    if rank.startswith("RankSPP"):
        df = df.sort_values(3, ascending=False)
    elif rank.startswith("RankCrowdness"):
        df = df.sort_values([4, 3], ascending=[True, False])
    elif rank.startswith("RankHotness"):
        df = df.sort_values([4, 3], ascending=[False, False])
    elif rank.startswith("RankLinear"):
        weight = float(rank.split("RankLinear")[1])

        df = df.sort_values(3, ascending=False)
        df[6] = range(1, len(df)+1) # spp rank

        df = df.sort_values([4, 3], ascending=[True, False])
        df[7] = range(1, len(df)+1) # c-score rank

        df[8] = df[6] + weight*df[7]
        df = df.sort_values([8, 3], ascending=[True, False]).iloc[:,:6] # combine rank

    df[6] = os.path.basename(infile).replace(".bed", "")
    df[7] = range(1, len(df)+1)
    df.to_csv(outfile, sep="\t", header=False, index=False)