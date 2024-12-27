import os
import pandas as pd
import tempfile
import pybedtools


def get_fasta(rank_file, genome_fasta, out_dir, allow_n_seq):
    
    df = pd.read_csv(rank_file, sep="\t", header=None)
    df = df[df[7] <= allow_n_seq].copy() # keep top n for motif inference
    df_g = df.groupby(6).max()
    tf2max = df_g[7].to_dict()

    for i,j in tf2max.items(): # if not enough sequences, keep top half
        if j < allow_n_seq:
            tf2max[i] = j//2
    
    df[8] = df[6].map(tf2max)
    df = df[df[7]<=df[8]].copy()
    
    motifs = df.iloc[:,-3].unique()

    with tempfile.TemporaryDirectory() as tmp_dir:
        for mm in motifs:
            df[df.iloc[:,-3] == mm].to_csv(os.path.join(tmp_dir, mm), sep="\t", header=False, index=False)
            bed = pybedtools.BedTool(os.path.join(tmp_dir, mm))
            bed.sequence(fi=genome_fasta).save_seqs(os.path.join(out_dir, mm+".fasta"))