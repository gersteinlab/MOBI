import re
import numpy as np
import pandas as pd
import numpy.ma as ma


def tomtom_cmd(out_dir, out_file, infer_meme, known_meme):
    cmd = "module load MEME; tomtom -no-ssc -oc %s -verbosity 1 -min-overlap 1 -dist pearson -evalue -thresh 10 %s %s; cp %s/tomtom.tsv %s\n" % (out_dir, infer_meme, known_meme, out_dir, out_file)
    return(cmd)


def get_inference_motif_order(motif_file):

    with open(motif_file, "r") as f:
        lines = f.readlines()

    ind = 0
    seq2ind = {}
    for line in lines:
        if line.startswith("MOTIF"):
            seq = line.split(" ")[1].rstrip()
            if seq not in seq2ind.keys():
                seq2ind[seq] = ind
            ind += 1
    return seq2ind


def get_inference_motif_array(mfile):
    motifs = []
    motif_names = []
    with open(mfile, "r") as f:
        lines = [line for line in f.readlines() if line.strip()]
        lines += ["ENDLINE"]
    for kk in range(len(lines)):
        if lines[kk].startswith("letter-probability"):
            for hh in range(kk+1, len(lines)):
                if not (lines[hh].lstrip())[0].isdigit():
                    break
            motif = np.array(lines)[kk+1:hh]

            motif_array = []
            for ii in motif:
                ii_strip = ii.rstrip().lstrip()
                motif_array.append(re.split("\s+", ii_strip)[:4])
            if motif.size != 0:
                motif = np.array(motif_array).astype(float)
                motifs.append(motif)
            else:
                motifs.append(np.nan)
        if lines[kk].startswith("MOTIF"):
            name = lines[kk].split(" ")[1].rstrip()
            motif_names.append(name)
    n2a = {}
    for ii,jj in zip(motifs, motif_names):
        if not np.all(np.isnan(ii)):
            n2a[jj] = ii
    return(n2a)


def motif_entropy(motif):
    log_term = ma.log(motif).filled(0)
    e_array = motif * log_term
    result = np.mean(-np.sum(e_array, axis=1))
    return(result)


def tomtom_result(meme_file, tomtom_file):
    n2i = get_inference_motif_order(meme_file) # name - order
    n2a = get_inference_motif_array(meme_file) # name - pwm array
    
    tomtom_df = pd.read_csv(tomtom_file, sep="\t", comment='#') # read tomtom result
    tomtom_df = tomtom_df.sort_values('p-value').drop_duplicates("Query_ID") # for one predict motif, only use the best matched known motif

    tomtom_df['topN'] = tomtom_df['Query_ID'].map(n2i) # order
    tomtom_df['-logp'] = -np.log(tomtom_df['p-value']) # -logP
    tomtom_df['entropy'] = tomtom_df['Query_ID'].map(n2a).map(motif_entropy) # motif entropy
    return(tomtom_df)