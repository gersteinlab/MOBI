"""Functions to prepare motifs into uniform format
"""

import numpy as np
from numpy import ma
import pandas as pd
from . import utils

def cisbp_meta_filter_in_vitro(file, outfile):
    """
    Filter motifs from cisbp based on meta data file

    Parameters
    ----------
    file : str
        TF_Information.txt download from cisbp
    outfile : str or None

    Returns
    -------
    DataFrame
    """

    cisbp_df = pd.read_csv(file, sep="\t")
    cisbp_df = cisbp_df[cisbp_df["TF_Status"] == "D"]  # direct
    cisbp_df = cisbp_df[cisbp_df["Motif_Type"].isin(
        ["SELEX", "PBM", "High-throughput SELEX SAGE", "B1H"])]  # exp type
    cisbp_df = cisbp_df.loc[:, ["TF_Name", "Motif_ID"]]
    cisbp_df["TF_Name"] = cisbp_df["TF_Name"].str.replace(
        "(", "").str.replace(")", "").str.replace(":", "")  # remove ():
    cisbp_df = cisbp_df.reset_index(drop=True).copy()

    if outfile:
        cisbp_df.to_csv(outfile, sep="\t", index=False)

    return cisbp_df


def stack_motif(sep_motif_dir, stack_motif_dir, split_str="--"):
    """Merge pwm files of multiple motifs of the same TF into a single file.
    Duplicated exact same pwm will be included only once.

    Parameters
    ----------------
    sep_motif_dir: str
        Every file contain only one motif. Files names are in format:
            [TF name][split_str][motif name].meme
    stack_motif_dir: str
        results in this dir would be [TF name].meme
    split_str: str
        see sep_motif_dir for file name format

    Returns
    ----------------
    None
    """

    all_motif = [i for i in os.listdir(sep_motif_dir) if i.endswith(
        ".meme") and not i.startswith(".")]
    u_motif = np.unique([i.split(split_str)[0] for i in all_motif])

    utils.os.mkdir_empty(stack_motif_dir, "Stack motif dir existed")

    # Stack motif of same TF
    # delete the redundant same motifs

    for i in u_motif:

        result = []
        for motif in all_motif:
            if motif.startswith(i+split_str):
                with open(os.path.join(sep_motif_dir, motif), "r") as f:
                    lines = f.readlines()
                    for k in range(len(lines)):
                        if lines[k].startswith("MOTIF"):
                            m_list = [m for m in lines[k:] if m.startswith(" ")]
                            final_str = ""
                            for m_str in m_list:
                                final_str += " "
                                final_str += "\t".join(m_str.split())
                                final_str += "\n"
                            result.append(final_str)

        u_ind = np.unique(result, return_index=True)[1]
        with open(os.path.join(stack_motif_dir, i+".meme"), "w") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")

            for motif_i in range(len(u_ind)):
                f.write("MOTIF %d\n" % motif_i)
                f.write("letter-probability matrix: alength= 4\n")
                f.write(result[u_ind[motif_i]])

    files = [i for i in os.listdir(stack_motif_dir) if i.endswith(".meme") and not i.startswith(".")]

    for i in files:
        with open(stack_motif_dir+i, "r") as f:
            content = f.read()
        name = i.replace(".meme", "")
        new_content = re.sub(r"MOTIF\s", "MOTIF %s_" % name, content)
        with open(stack_motif_dir+i, "w") as f:
            f.write(new_content)
    return


def motif_entropy(motif):
    """calculate motif entropy.

    Parameters
    ----------
    motif : np.array
        motif pwm

    Returns
    -------
    float
        entropy
    """
    log_term = ma.log(motif).filled(0)
    e_array = motif * log_term
    return(-e_array.sum())
