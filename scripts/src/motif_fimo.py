import os
import numpy as np
from . import utils

def run_FIMO(motif_dir, genome_file, fimo_output_dir, joblist_file, bgfreq=False, subset=None):
    """Run FIMO
    Generate slurm joblist file

    Parameters
    ----------
    motif_dir : str
        dir of motif pwm files, current naming: [TF name].meme
    genome_file : str
    fimo_output_dir : str
    joblist_file : str
    bgfreq : boolean
        whether to run fimo with "--bgfile --motif--"
    subset : list or np.array
        list of TF names, current naming: [TF name],
        only TF in this subset will be run

    Returns
    -------
    None
    """

    utils.os.mkdir_empty(fimo_output_dir, msg="skip fimo")

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
    """Reformat the raw FIMO output to bedfile

    Parameters
    ----------
    fimor_output_dir : str
        the direct output directory from fimo, current structure: fimo_output_dir/[TF name].meme/fimo.txt
    output_dir : str
        current structure: output_fir/[TF name].bed

    Returns
    -------
    """

    utils.os.mkdir_empty(output_dir, msg="skip reformat_fimo")

    motif_files = [i.replace(".meme", "") for i in os.listdir(fimo_output_dir) if i.endswith(".meme") and not i.startswith(".")]
    for motif in motif_files:
        os.system("sed '1d' %s/%s.meme/fimo.txt |\
                   awk '{print $0\"\\t\"$1\"\\t\"$2}' |\
                   cut -f3- > %s/%s.bed" %
                  (fimo_output_dir, motif, output_dir, motif))
    return


def filter_FIMO_result(input_dir, output_dir, fimo_p=0.0001, fimo_q=None):
    """Filter to get significant motif hit

    Parameters
    ----------
    input_dir : str
    output_dir : str
    fimo_p : float
        use p-value as filter, not compatible with fimo_q
    fimo_q : float
        use q-value as filter, not compatible with fimo_p

    Returns
    -------
    None
    """

    utils.os.mkdir_empty(output_dir, msg="filter fimo")

    motif_files = [i for i in os.listdir(input_dir) if i.endswith(".bed") and not i.startswith(".")]
    for motif in motif_files:
        if fimo_p:
            os.system("awk '{if ($6<%f){print $0}}' %s/%s > %s/%s" %
                      (fimo_p, input_dir, motif, output_dir, motif))
        elif fimo_q:
            os.system("awk '{if ($7<%f){print $0}}' %s/%s > %s/%s" %
                      (fimo_q, input_dir, motif, output_dir, motif))
    return
