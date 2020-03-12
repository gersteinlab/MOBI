"""
Functions related to motifs. Filter, preprocess of motifs from CIS-bp, run FIMO etc.
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from src.utils import os_mkdir_empty

__all__ = ["cisbp_meta_filter_in_vitro",
           "stack_motif",
           "run_FIMO",
           "reformat_FIMO_result_to_bed",
           "filter_FIMO_result"]

def cisbp_meta_filter_in_vitro(file, outfile):
    """
    Filter motifs from cisbp based on meta data file

    Parameters
    ----------------
    file: str
        TF_Information.txt download from cisbp
    outfile: str or None
    
    Returns
    ----------------
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
    """
    Merge pwm files of multiple motifs of the same TF into a single file.
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

    if not os_mkdir_empty(stack_motif_dir, ""):
        raise OSError("mkdir fail. Stack motif dir existed.")

    # Stack motif of same TF
    # delete the redundant same motifs

    for i in u_motif:

        result = []
        for motif in all_motif:
            if motif.startswith(i+split_str):
                with open("%s/%s" % (sep_motif_dir, motif), "r") as f:
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
        with open("%s/%s.meme" % (stack_motif_dir, i), "w") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")

            for motif_i in range(len(u_ind)):
                f.write("MOTIF %d\n" % motif_i)
                f.write("letter-probability matrix: alength= 4\n")
                f.write(result[u_ind[motif_i]])

    files = [i for i in os.listdir(stack_motif_dir) if i.endswith(
        ".meme") and not i.startswith(".")]

    for i in files:
        with open(stack_motif_dir+i, "r") as f:
            content = f.read()
        name = i.replace(".meme", "")
        new_content = re.sub(r"MOTIF\s", "MOTIF %s_" % name, content)
        with open(stack_motif_dir+i, "w") as f:
            f.write(new_content)
    return

def run_FIMO(motif_dir, genome_file, fimo_output_dir, joblist_file, bgfreq=False, subset=None):
    """
    Get the command line joblist file of runing FIMO
    
    Parameters
    ----------------
    motif_dir: str
        dir of motif pwm files, current naming: [TF name].meme
    genome_file: str
    fimo_output_dir: str
    joblist_file: str
    bgfreq: boolean
        whether to run fimo with "--bgfile --motif--"
    subset: list or np.array
        list of TF names, current naming: [TF name],
        only TF in this subset will be run

    Returns
    ----------------
    None
    """

    if not os_mkdir_empty(fimo_output_dir, msg="skip fimo"):
        return

    if isinstance(subset, (list, np.ndarray)):
        motif_files = [motif_dir+i for i in os.listdir(motif_dir) if i.endswith(".meme") and not i.startswith(".") and (i.replace(".meme", "") in subset)]
    else:
        motif_files = [motif_dir+i for i in os.listdir(motif_dir) if i.endswith(".meme") and not i.startswith(".")]

    with open(joblist_file, "w") as f:
        for motif in motif_files:
            if bgfreq:
                f.write("fimo --bgfile --motif-- --thresh 0.001 --max-stored-scores 10000000 --o %s/%s/ %s %s\n" %
                        (fimo_output_dir, os.path.basename(motif), motif, genome_file))
            else:
                f.write("fimo --thresh 0.001 --max-stored-scores 10000000 --o %s/%s/ %s %s\n" %
                        (fimo_output_dir, os.path.basename(motif), motif, genome_file))
    return


def reformat_FIMO_result_to_bed(fimo_output_dir, output_dir):
    """
    Reformat the raw FIMO output to bedfile
    
    Parameters
    ----------------
    fimor_output_dir: str
        the direct output directory from fimo, current structure: fimo_output_dir/[TF name].meme/fimo.txt
    output_dir: str
        current structure: output_fir/[TF name].bed

    Returns
    ----------------
    """

    if not os_mkdir_empty(output_dir, msg="skip clean fimo"):
        return

    motif_files = [i.replace(".meme", "") for i in os.listdir(
        fimo_output_dir) if i.endswith(".meme") and not i.startswith(".")]
    for motif in motif_files:
        os.system("sed '1d' %s/%s.meme/fimo.txt |\
                   awk '{print $0\"\\t\"$1\"\\t\"$2}' |\
                   cut -f3- > %s/%s.bed" %
                  (fimo_output_dir, motif, output_dir, motif))
    return

def filter_FIMO_result(input_dir, output_dir, fimo_p=0.0001, fimo_q=None):
    """
    Filter to get significant motif hit

    Parameters
    ----------------
    input_dir: str
    output_dir: str
    fimo_p: float
        use p-value as filter, not compatible with fimo_q
    fimo_q: float
        use q-value as filter, not compatible with fimo_p

    Returns
    ----------------
    None
    """

    if not os_mkdir_empty(output_dir, msg="skip filter fimo"):
        return

    motif_files = [i for i in os.listdir(input_dir) if i.endswith(
        ".bed") and not i.startswith(".")]
    for motif in motif_files:
        if fimo_p:
            os.system("awk '{if ($6<%f){print $0}}' %s/%s > %s/%s" %
                      (fimo_p, input_dir, motif, output_dir, motif))
        elif fimo_q:
            os.system("awk '{if ($7<%f){print $0}}' %s/%s > %s/%s" %
                      (fimo_q, input_dir, motif, output_dir, motif))
    return