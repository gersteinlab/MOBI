import os
import shutil
import uuid
import numpy as np
import pandas as pd
from src.utils import *
from src.narrowPeak import *

__all__ = [
    "get_all_summit_file",
    "get_crowdness_sliding_window",
    "get_motif_hit",
    "get_rank_file",
    "get_motif_enrichment",
    "summarize_motif_enrichment",
    "get_bed_from_rank_file",
    "get_fasta_from_bed",
    "run_inference_joblist",
    "run_tomtom_joblist",
    "inference_output_to_minimal_meme",
    "get_tomtom_hits",
    "get_tomtom_summary_data"]


def get_all_summit_file(TF_files, all_summit_file):
    """
    Given a list of TF bs, fetch the summits of bs and save into a single file
    output: chr, summit start, summit end

    Parameters
    ----------------
    TF_files: list of str
        ENCODE bed.gz files
        $1 chr, $2 peak start, $10 peak summit
    all_summit_file: str

    Returns
    ----------------
    None
    """

    if os.path.exists(all_summit_file):
        print("skip get all_summit_file")
        return

    tmp_dir = os_mkdir_tmp()

    for tf in TF_files:
        get_fix_width_region(
            infile=tf, 
            outfile=os.path.join(tmp_dir, os.path.basename(tf)),
            width=-1,
            summit_col=10,
            keep_col=None)

    os.system("cat %s/* | cut -f-3 | sort -k1,1 -k2,2n > %s" % (tmp_dir, all_summit_file))
    
    shutil.rmtree(tmp_dir)
    return


def get_crowdness_sliding_window(TF_files, result_dir, all_summit_file=None, crowdness_n_sliding=250):
    """
    Get crowdness score by sliding window method
    input: ENCODE bed.gz files
    input: $1 chr, $2 peak start, $10 peak summit, $7 SPP
    output: $1-3 coor of crowdness_n_crowdness bp around each peak summit, $4 spp, $5 crowdness

    For each summit of a TF, the crowdness is
    how many other summits (from all TF) lies within
    crowdness_n_sliding bp window
    line 3: if start < 0, set to 0 (<0 not possible in genome coordinate)
    line 5: count summits
    
    Parameters
    ----------------
    TF_files: list of str
    result_dir: str
    all_summit_file: str or None
        if not given, generate a new all_summit_file
    crowdness_n_sliding: int
        actual length when used to calculate the crowdness
        this will also be the binding site length of the output

    Returns
    ----------------
    None
    """

    if not os_mkdir_empty(result_dir, "skip get_crowdness"):
        return

    tmp_dir1 = os_mkdir_tmp()
    tmp_dir2 = os_mkdir_tmp()

    for tf in TF_files:
        # trim to the window size
        get_fix_width_region(
            infile=tf, 
            outfile=os.path.join(tmp_dir1, os.path.basename(tf)), 
            width=crowdness_n_sliding,
            summit_col=10,
            keep_col=7)
        # remove abnormal coordinates
        rm_column_negative(
            infile=os.path.join(tmp_dir1, os.path.basename(tf)), 
            outfile=os.path.join(tmp_dir2, os.path.basename(tf)), 
            col=2)
        # intersect with the all_summit files to count summits in the window
        os.system("\
            sort -k1,1 -k2,2n %s|\
            bedtools intersect -c -a stdin -b %s|\
            awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5}' \
            > %s" % (
                os.path.join(tmp_dir2, os.path.basename(tf)),
                all_summit_file,
                os.path.join(result_dir, os.path.basename(tf))))
    
    shutil.rmtree(tmp_dir1)
    shutil.rmtree(tmp_dir2)

    return


def get_motif_hit(TF_files, motif_fimo_files, result_dir):
    """
    Count motif occurance in each binding site
    If motif_fimo_files is given, then this is done via intersecting the FIMO result
    If not, assign pseudo count, which is all 0

    output: 
        add col: motif_count
        add col: TF name

    Parameters
    ----------------
    TF_files: list of str
        bed file
    motif_fimo_files: list of str, or None
        corresponde to TF_files
    result_dir: str

    Returns
    ----------------
    None
    """

    if not os_mkdir_empty(result_dir, "skip get_motif_hit_fimo_intersect"):
        return

    if motif_fimo_files:
        for tf, motif_fimo in zip(TF_files, motif_fimo_files):
            #os.system("bedtools intersect -nonamecheck -c -a %s -b %s |\
            os.system("bedtools intersect -c -a %s -b %s |\
                       sed -e \"s/$/\t%s/g\" \
                       > %s/%s"
                      % (tf, motif_fimo,
                          os.path.basename(tf).replace(".bed", ""),
                          result_dir, os.path.basename(tf)))
    else:
        for tf in TF_files:
            with open(tf, "r") as f1:
                with open(result_dir+os.path.basename(tf), "w") as f2:
                    for line in f1:
                        new_line = line.replace("\n", "\t0\t%s\n" % os.path.basename(tf).replace(".bed", ""))
                        f2.write(new_line)
    return


def get_rank_file(rank_method, TF_files, result_file, bs_half_length=200, rank_linear_alpha=None, rank_exclude_prop=None):
    """
    Rank all binding sites according to the ranking criteria
    All TF bs into a single file
    input: $4 spp, $5 crowdness, $6 motif_count(arbitrary), $7 TF name
    output: $1 chr, $2 summit, $3-6 same as $4-7 in input, $8 rank

    Parameters
    ----------------
    rank_method: {"RankSPP", "RankCrowdness", "RankHotness", "RankLinear_*", "RankExcludeArbitrary_*"}
    TF_files: list of str
    result_file: str
    bs_half_length: int
    rank_linear_alpha: int or None
    rank_exclude_prop: float or None

    Returns
    ----------------
    None
    """

    if os.path.exists(result_file):
        print("skip get_rank_file")
        return

    if rank_method.startswith("RankSPP"):
        for tf in TF_files:
            os.system("awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                       sort -k3,3nr  | \
                       awk '{print $0\"\\t\"NR}' \
                       >> %s"
                      % (bs_half_length, tf,
                          result_file))

    elif rank_method.startswith("RankCrowdness"):
        for tf in TF_files:
            os.system("awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                       sort -k4,4n -k3,3nr  | \
                       awk '{print $0\"\\t\"NR}' \
                       >> %s"
                      % (bs_half_length, tf,
                          result_file))

    elif rank_method.startswith("RankHotness"):
        for tf in TF_files:
            os.system("awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                       sort -k4,4nr -k3,3nr  | \
                       awk '{print $0\"\\t\"NR}' \
                       >> %s"
                      % (bs_half_length, tf,
                          result_file))

    elif rank_method.startswith("RankLinear"):
        for tf in TF_files:
            os.system("awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                       sort -k3,3nr  | awk '{print $0\"\\t\"NR}' | \
                       sort -k4,4n -k3,3nr | awk '{print $0\"\\t\"NR}' | \
                       awk '{print $0\"\\t\"$7+$8*%g}' | \
                       sort -k9,9n | \
                       awk '{print $0\"\\t\"NR}' | \
                       cut -f-6,10- \
                       >> %s"
                      % (bs_half_length, tf,
                          rank_linear_alpha,
                          result_file))

    elif rank_method.startswith("RankExcludeArbitrary"):

        tmp_file = os_mkdir_tmp()

        for tf in TF_files:
            os.system("cat %s >> %s" % (tf, tmp_file))  # all peaks from all TF
        df = pd.read_csv(tmp_file, sep="\t", header=None)
        # index = total_n*proportion
        ind = len(df[4]) - int(np.round(len(df[4]) * rank_exclude_prop))
        # crowdness_rank[index] is the crowdness threshold to keep the peak
        thresh = np.sort(df[4].values)[ind]
        os.remove(tmp_file)

        # remove all peaks that are hotter than the threshold
        for tf in TF_files:
            os.system("awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                       sort -k3,3nr | \
                       awk '{if($4<=%d){print $0}}' | \
                       awk '{print $0\"\\t\"NR}' \
                       >> %s"
                      % (bs_half_length, tf,
                          thresh,
                          result_file))

    else:
        raise Exception("Unknown rank method")
    return


def get_motif_enrichment(rank_file, n_seq, result_dir, enrichment_method="EnrichmentFraction"):
    """
    Get motif enrichment
    Fraction of binding sites that have motif count >=1 
    HalfInsuff: if not enough sequence, use the top half only

    Parameters
    ----------------
    rank_file: str
    n_seq: int
    result_dir: str
    enrichment_method: {"EnrichmentFraction", "EnrichmentFractionHalfInsuff"}

    Returns
    ----------------
    None
    """
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    if os.path.exists("%s/%s_%d.txt" % (result_dir, enrichment_method, n_seq)):
        print("skip %s_%d.txt" % (enrichment_method, n_seq))
        return

    tmp_file = "/home/jg2447/log/TMP_%d" % uuid.uuid4()

    # keep all peaks with ranking < n_seq
    # then group by TF
    # then calculate motif enrichment

    enrichment = pd.DataFrame()

    os.system("awk '{if($7<=%d){print $0}}' %s > %s" %
              (n_seq, rank_file, tmp_file))

    df = pd.read_csv(tmp_file, sep="\t", header=None)
    os.remove(tmp_file)

    if enrichment_method == "EnrichmentFractionHalfInsuff":
        tfs = df[5].unique()
        new_df = pd.DataFrame()
        for tf in tfs:
            tf_df = df[df[5] == tf].copy()
            if len(tf_df) < n_seq:
                tf_df = tf_df.iloc[:len(tf_df)//2, :].copy()
            new_df = pd.concat([new_df, tf_df])
        df = new_df.reset_index(drop=True).copy()
    elif enrichment_method == "EnrichmentFraction":
        pass
    elif enrichment_method == "EnrichmentAvgHitPerBS":
        raise Exception("TBD, not implemented yet")
    else:
        raise Exception("Unknown enrichment method")

    df[4] = df[4].replace(0, np.nan)
    df = df.groupby(5).count()
    n_seq_result = (df[4] / df[0])

    n_seq_result.to_csv("%s/%s_%d.txt" %
                        (result_dir, enrichment_method, n_seq), sep="\t", header=False, float_format="%.3f")   
    return


def summarize_motif_enrichment(enrichment_dir):
    """
    Compress and summarize enrichment results

    Parameters
    ----------------
    enrichment_dir: str
        the main dir of enrichment, structure is main/SubFolderByRank/fileByNseq.txt

    Returns
    ----------------
    None
    """

    enrichment_rank_dirs = [i for i in os.listdir(enrichment_dir) if i.startswith("Rank")]

    for rank_dir in enrichment_rank_dirs:

        edir = enrichment_dir + rank_dir + "/"
        files = [i for i in os.listdir(edir) if i.startswith("EnrichmentFraction")]

        df = pd.read_csv(edir+files[0], sep="\t", header=None, index_col=0)
        del df[1]

        for f in files:
            df2 = pd.read_csv(edir+f, sep="\t", header=None, index_col=0)
            df = pd.merge(df, df2, left_index=True, right_index=True, how='outer')

        df.columns = [i.replace(".txt", "").replace("EnrichmentFraction_", "") for i in files]
        df.columns = df.columns.astype(int)
        df = df.reindex(sorted(df.columns), axis=1)
        df.to_csv(enrichment_dir+rank_dir+".txt", sep="\t", header=True, index=True, float_format="%.3f")
        
    for rank_dir in enrichment_rank_dirs:
        edir = enrichment_dir + rank_dir + "/"   
        shutil.rmtree(edir)

    return


def get_bed_from_rank_file(rank_file, inference_n_seq, inference_nbp, result_dir):
    """
    split rank file by TF name into individual BED file
    keep top 2*inference_n_seq sequences(if not enough, keep all)
    each region truncated to +- inference_nbp around summit

    Parameters
    ----------------
    rank_file: str
    inference_n_seq: int
    inference_nbp: int
    result_dir: str

    Returns
    ----------------
    None
    """

    os_mkdir_empty(result_dir, "skip split_rank_file")

    allow_n_seq = 2 * inference_n_seq
    os.system("awk '{if($7<=%d){print $1\"\\t\"$2-%d\"\\t\"$2+%d >\"%s/\"$6\".bed\"}}' %s" \
              % (allow_n_seq, inference_nbp, inference_nbp, result_dir, rank_file))

    return


def get_fasta_from_bed(bed_dir, genome_file, result_dir):
    """
    convert BED to fasta
    only keep top half sequences

    Parameters
    ----------------
    bed_dir: str
    genome_file: str
    result_dir: str

    Returns
    ----------------
    None
    """

    if os.path.exists(result_dir):
        print("skip get_fasta")
        return
    os.makedirs(result_dir)
    
    files = [i.split(".bed")[0] for i in os.listdir(bed_dir) if i.endswith(".bed")]
        
    for fi in files:
        with open("%s/%s.bed" % (bed_dir, fi), "r") as f:
            flen = len(f.readlines())
        os.system("head -n %d %s/%s.bed |\
                   bedtools getfasta -fi %s -bed stdin -fo %s/%s.fasta" \
                   % (flen//2, bed_dir, fi,\
                   genome_file, result_dir, fi))
    
    return


def run_inference_joblist(inference_tool, fasta_dir, hpc_joblist_file, inference_motif_dir, keep_motif_only=True):
    """
    Run inference tools via slurm
    This function only generate sbatch joblist

    Parameters
    ----------------
    inference_tool: {"DREME", "MEME", "HOMER"}
    fasta_dir: str
    hpc_joblist_file: str
    inference_motif_dir: str
    keep_motif_only: bool
        if true, only keep the meme file
        if false, keep all the output from the inference tools

    Returns
    ----------------
    None
    """

    os_mkdir_empty(inference_motif_dir, "run_inference output dir exist")
      
    tfs = [i.split(".fasta")[0] for i in os.listdir(fasta_dir) if i.endswith(".fasta")]

    with open(hpc_joblist_file, "w") as f:
        for tf in tfs:
            if keep_motif_only:
                outdir = inference_motif_dir.rstrip("/")
                if inference_tool == "DREME":
                    f.write("dreme -oc %s/%s/ -k 8 -m 10 -e 1000 -p %s/%s.fasta; mv %s/%s/dreme.txt %s_%s.meme; rm -rf %s/%s/\n" \
                            % (inference_motif_dir, tf, fasta_dir, tf, outdir, tf, inference_motif_dir, tf))
                elif inference_tool == "HOMER":
                    f.write("findMotifs.pl %s/%s.fasta fasta %s/%s/ -nocheck -nogo -len 8 -noknown -basic -S 10; Rscript /home/jg2447/slayman/motif_inference/scripts/src/MoVRs_Motif2meme.R %s/%s/homerMotifs.motifs8 %s/%s/homerMotifs.meme; ; mv %s/%s/homerMotifs.meme %s_%s.meme; rm -rf %s/%s/\n" \
                            % (fasta_dir, tf, inference_motif_dir, tf, inference_motif_dir, tf, inference_motif_dir, tf, outdir, tf, inference_motif_dir, tf))
                elif inference_tool == "MEME":
                    f.write("meme -oc %s/%s/ -dna -nmotifs 10 -w 8 -maxsize 150000 -nostatus %s/%s.fasta; mv %s/%s/meme.txt %s_%s.meme; rm -rf %s/%s/\n" \
                            % (inference_motif_dir, tf, fasta_dir, tf, outdir, tf, inference_motif_dir, tf))
            else:    
                if inference_tool == "DREME":
                    f.write("dreme -oc %s/%s/ -k 8 -m 10 -e 1000 -p %s/%s.fasta\n" \
                            % (inference_motif_dir, tf, fasta_dir, tf))
                elif inference_tool == "HOMER":
                    f.write("findMotifs.pl %s/%s.fasta fasta %s/%s/ -nocheck -nogo -len 8 -noknown -basic -S 10; Rscript /home/jg2447/slayman/motif_inference/scripts/src/MoVRs_Motif2meme.R %s/%s/homerMotifs.motifs8 %s/%s/homerMotifs.meme\n" \
                            % (fasta_dir, tf, inference_motif_dir, tf))
                elif inference_tool == "MEME":
                    f.write("meme -oc %s/%s/ -dna -nmotifs 10 -w 8 -maxsize 150000 -nostatus %s/%s.fasta\n" \
                            % (inference_motif_dir, tf, fasta_dir, tf))


    return


def run_tomtom_joblist(inference_tool, fasta_dir, hpc_joblist_file, inference_motif_dir, known_motif_dir, inference_tomtom_dir, keep_motif_only=True):
    """
    Run tomtom via slurm
    This function only generate sbatch joblist

    Parameters
    ----------------
    inference_tool: {"DREME", "MEME", "HOMER"}
    fasta_dir: str
    hpc_joblist_file: str
    inference_motif_dir: str
    known_motif_dir: str
    inference_tomtom_dir: str

    Returns
    ----------------
    None
    """
    os_mkdir_empty(inference_tomtom_dir, "run_tomtom output dir exist")
    
    tfs = [i.split(".fasta")[0] for i in os.listdir(fasta_dir) if i.endswith(".fasta")]
    
    with open(hpc_joblist_file, "w") as f:
        for tf in tfs:
            if keep_motif_only:
                outdir = inference_motif_dir.rstrip("/")
                outdir2 = inference_tomtom_dir.rstrip("/")
                f.write("tomtom -no-ssc -oc %s/%s/ -verbosity 1 -min-overlap 1 -dist pearson -evalue -thresh 10 %s_%s.meme %s/%s.meme; mv %s/%s/tomtom.txt %s_%s.tomtom; rm -rf %s/%s/\n" \
                        % (inference_tomtom_dir, tf, outdir, tf, known_motif_dir, tf, inference_tomtom_dir, tf, outdir2, tf, inference_tomtom_dir, tf))
            else:
                f.write("tomtom -no-ssc -oc %s/%s/ -verbosity 1 -min-overlap 1 -dist pearson -evalue -thresh 10 %s/%s/dreme.txt %s/%s.meme\n" \
                        % (inference_tomtom_dir, tf, inference_motif_dir, tf, known_motif_dir, tf))
    return


def inference_output_to_minimal_meme(infile, outfile, bg_freq=[0.25,0.25,0.25,0.25], tool="DREME"):
    """
    Reformat the output from inference tools to the same minimal meme format

    Parameters
    ----------------
    infile: str
    outfile: str
    bg_freq: list of 4 float
        backgroud frequency, need to sum up to 1, A match T and C match G
    tool: {"DREME", "MEME", "HOMER"}

    Returns
    ----------------
    None
    """


    with open(infile, "r") as f:
        c = f.read()
    motif_file = c.split("\n")
    motif_file = [i.lstrip() for i in motif_file]

    name_id = 1
    seq_list = []
    matrix_info_line_list = []
    pwm_list = []
    for i in range(len(motif_file)):
        if motif_file[i].startswith("MOTIF"):
            if (tool == "DREME") or (tool == "MEME"):
                id_line = motif_file[i].split()
                seq_list.append(id_line[1])
            elif tool == "HOMER":
                id_line = motif_file[i].split()
                seq = id_line[2].split("-")[1]
                seq_list.append(seq)
            continue
        if motif_file[i].startswith("letter-probability matrix"):
            matrix_info_line_list.append(motif_file[i])
            j = i + 1
            pwm = []
            while motif_file[j].startswith("0") or motif_file[j].startswith("1"):
                pwm.append(re.sub(r"\s+", "\t", motif_file[j]))
                j += 1
            pwm_list.append(pwm)
            continue

    name_list = [tool+"-"+str(i) for i in range(1,len(seq_list)+1)]

    with open(outfile, "w") as out_f:
        out_f.write("""MEME version 4

ALPHABET "DNA" DNA-LIKE
A "Adenine" CC0000 ~ T "Thymine" 008000
C "Cytosine" 0000CC ~ G "Guanine" FFB300
N "Any base" = ACGT
X = ACGT
. = ACGT
V "Not T" = ACG
H "Not G" = ACT
D "Not C" = AGT
B "Not A" = CGT
M "Amino" = AC
R "Purine" = AG
W "Weak" = AT
S "Strong" = CG
Y "Pyrimidine" = CT
K "Keto" = GT
U = T
END ALPHABET

strands: + -

Background letter frequencies:
A %f C %f G %f T %f

""" % (bg_freq[0], bg_freq[1], bg_freq[2], bg_freq[3]))
        for i in range(len(seq_list)):
            out_f.write("MOTIF %s %s\n" % (seq_list[i], name_list[i]))
            out_f.write("%s\n" % (matrix_info_line_list[i]))
            out_f.write("%s\n\n" % ("\n".join(pwm_list[i])))
    return


def _get_inference_motif_order(motif_file):
    """
    seq2ind:
    for each motif sequence,
    return the order of that sequence in the inferred motif file (single pwm with multiple motifs)

    Parameters
    ----------------
    motif_file: str

    Returns
    ----------------
    dict
    """

    with open(motif_file, "r") as f:
        lines = f.readlines()
    
    ind = 0
    seq2ind = {}
    for line in lines:
        if line.startswith("MOTIF"):
            seq = line.split(" ")[1]
            seq2ind[seq] = ind
            ind += 1
    return seq2ind


def get_tomtom_hits(tomtom_dir, motif_dir, tool, rank, nn, match_type, n_tomtom_motif, tf_subset=None):
    """
    Parse the tomtom output directory and report the results
    
    Parameters
    ----------------
    tomtom_dir: str
    motif_dir: str
    tool: str
    rank: str
    nn: int
    match_type: {"count", "order", "p-value", "order_p-value", "known_hit"}
    n_tomtom_motif: list of int
        the nth infered motif to be counted, [0] is top 1 and [0,1,2,3,4] is top 5
    
    Returns
    ----------------
    np.array
        length of array = number of TF
    """

    # get available TFs
    tfs = [i for i in os.listdir(tomtom_dir) if i.startswith("%s_%s_%s" % (tool, rank, nn))]
    if rank.startswith("RankLinear"):
        tfs = [i.split("_")[4] for i in tfs]
    else:
        tfs = [i.split("_")[3] for i in tfs]
    tfs = [i.replace(".tomtom", "") for i in tfs]
    tfs = np.unique(tfs)

    # if given subset, only use TFs in the subset
    if tf_subset:
        tfs = np.array([i for i in tfs if i in tf_subset])
    
    result_tf = []
    result = []
    for i in range(len(tfs)):

        # get inferred motifs order
        seq2ind_DREME = _get_inference_motif_order("%s/DREME_%s_%s_%s.meme" % (motif_dir, rank, nn, tfs[i]))
        seq2ind_HOMER = _get_inference_motif_order("%s/HOMER_%s_%s_%s.meme" % (motif_dir, rank, nn, tfs[i]))
        seq2ind_MEME = _get_inference_motif_order("%s/MEME_%s_%s_%s.meme" % (motif_dir, rank, nn, tfs[i]))

        # use the min number of inference motif across all tools' result
        n_min_motif = min(len(seq2ind_DREME), len(seq2ind_HOMER), len(seq2ind_MEME))
        n_motif = np.array(n_tomtom_motif)[np.array(n_tomtom_motif) < n_min_motif]
        # skip if no motifs
        if len(n_motif) == 0:
            continue

        # focus on the one currently been run
        if tool == "DREME":
            seq2ind = seq2ind_DREME
        elif tool == "HOMER":
            seq2ind = seq2ind_HOMER
        elif tool == "MEME":
            seq2ind = seq2ind_MEME      

        # read raw tomtom output
        df = pd.read_csv("%s/%s_%s_%s_%s.tomtom" % (tomtom_dir, tool, rank, nn, tfs[i]), sep="\t")
        try:
            spp = pd.read_csv("%s/%s_RankSPP_100_%s.tomtom" % (tomtom_dir, tool, tfs[i]), sep="\t")
        except Exception:
            continue
        
        # keep P-value < 0.05 if the statistics is not relative to p-value
        if (match_type != "p-value") and (match_type != "order_p-value"):
            df = df[df['p-value'] <= 0.05]

        # for each query, only keep the known motif match with most significant level (i.e. smallest p-value)
        if df.shape[0] != 0:
            # match_df = df[df['Target ID'].str.startswith(tfs[i])].copy() # use this line if target is mix of motifs
            match_df = df.copy() # use this line if target is single motif
            match_df = match_df.sort_values('p-value').drop_duplicates("#Query ID") 
        else:
            match_df = pd.DataFrame([])

        # keep the top n or nth motif(define by range n_tomtom_motif)
        if match_df.shape[0] != 0:
            match_ind = [seq2ind[i] for i in match_df["#Query ID"].values]
            match_df_ind = [i for i in range(len(match_ind)) if match_ind[i] in n_tomtom_motif]
            match_df = match_df.iloc[match_df_ind,:].copy()

        # results base on match_type
        
        # for no hits
        if match_df.shape[0] == 0:
            result.append(0)    
        else:
            match_ind = [seq2ind[i] for i in match_df["#Query ID"].values]
            match_ind_pt = 10 - np.array(match_ind)
            p_value_pt = -np.log(match_df['p-value'].values)

            if match_type == "count":
                result.append(match_df.shape[0])
            elif match_type == "order":
                result.append(np.sum(match_ind_pt))
            elif match_type == "p-value":
                result.append(np.sum(p_value_pt))
            elif match_type == "order_p-value":
                result.append(np.sum(match_ind_pt * p_value_pt))
            elif match_type == "known_hit":
                result.append(len(match_df["Target ID"].unique()))
        
        # record the used TF
        # some TF might be skip b/c no motifs
        result_tf.append(tfs[i])

    return((np.array(result_tf), np.array(result)))


def get_tomtom_summary_data(tomtom_summary_dir, nbp_list, rank_list, tool, top, mtype="count", improve=True, tf_subset=None):
    """
    Retrieve tomtom summary data of different experiments

    Parameters
    ----------------
    tomtom_summary_dir: str
    nbp_list: list of int
    rank_list: list of str
    tool: {"DREME", "MEME", "HOMER"}
    top: {"top5", "top1"}
    improve: boolean
        if true return improvement comparing to the convention, which is RankSPP and nbp=100
        require "RankSPP" in rank_list and 100 in nbp_list
        else return raw count
    tf_subset: list or np.array
        if given, only return results based on TFs in this subset

    Returns
    ----------------
    DataFrame
    """

    nn_data = []
    for nn in nbp_list:
        rank_data = []
        for rank in rank_list:
            
            data = pd.read_csv(
                    "%s/%s_%s_%d_%s_%s.txt" % (
                        tomtom_summary_dir, tool, rank, nn, mtype, top),
                    sep="\t", header=None)

            if tf_subset is not None:    
                data = data[data[0].isin(tf_subset)][1].values
            else:
                data = data[1].values
            
            rank_data.append(data.mean()) # mean of all TF, then all rank
        nn_data.append(rank_data) # then all nn
    
    df = pd.DataFrame(nn_data)
    df.columns = rank_list
    df.index = nbp_list
    
    # get improvement    
    if improve:
        df = (df - df.loc[100, "RankSPP"]) / df.loc[100, "RankSPP"]
    
    return df