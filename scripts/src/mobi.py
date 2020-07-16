
import os
import shutil
import numpy as np
import pandas as pd
import utils
import narrow_peak

def get_all_summit_file(TF_files, all_summit_file):
    """Given a list of TF bs, fetch the summits of bs and save into a single file
    output: chr, summit start, summit end

    Parameters
    ----------
    TF_files : list of str
        ENCODE bed.gz files
        $1 chr, $2 peak start, $10 peak summit
    all_summit_file : str

    Returns
    -------
    None
    """

    tmp_dir = utils.os.mkdir_tmp()

    for tf in TF_files:
        narrow_peak.get_fix_width_region(
            infile=tf,
            outfile=os.path.join(tmp_dir, os.path.basename(tf)),
            width=-1,
            summit_col=10,
            keep_col=None)

    # concat all into one file
    # keep only coordinates
    # sort
    os.system("cat %s/* | cut -f-3 | sort -k1,1 -k2,2n > %s" % (tmp_dir, all_summit_file))

    shutil.rmtree(tmp_dir)
    return


def get_crowdness_sliding_window(TF_files, result_dir, all_summit_file, crowdness_n_sliding=250):
    """Get crowdness score by sliding window method
    input: ENCODE bed.gz files
    input: $1 chr, $2 peak start, $10 peak summit, $7 SPP
    output: $1-3 coor of crowdness_n_crowdness bp around each peak summit, $4 spp, $5 crowdness

    For each summit of a TF, the crowdness is
    how many other summits (from all TF) lies within
    crowdness_n_sliding bp window
    line 3: if start < 0, set to 0 (<0 not possible in genome coordinate)
    line 5: count summits

    Parameters
    ----------
    TF_files : list of str
    all_summit_file : str
    result_dir : str
    crowdness_n_sliding : int
        actual length when used to calculate the crowdness
        this will also be the binding site length of the output

    Returns
    -------
    None
    """

    tmp_dir1 = utils.os_mkdir_tmp()
    tmp_dir2 = utils.os_mkdir_tmp()

    for tf in TF_files:
        # trim to the window size
        narrow_peak.get_fix_width_region(
            infile=tf,
            outfile=os.path.join(tmp_dir1, os.path.basename(tf)),
            width=crowdness_n_sliding,
            summit_col=10,
            keep_col=7)
        # remove abnormal coordinates
        narrow_peak.rm_column_negative(
            infile=os.path.join(tmp_dir1, os.path.basename(tf)),
            outfile=os.path.join(tmp_dir2, os.path.basename(tf)),
            col=2)
        # intersect with the all_summit files to count summits in the window
        os.system("
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
    """Count motif occurance in each binding site
    If motif_fimo_files is given, then this is done via intersecting the FIMO result
    If not, assign pseudo count, which is all 0

    output format:
        add col: motif_count
        add col: TF name

    Parameters
    ----------
    TF_files : list of str
        bed file
    motif_fimo_files : list of str, or None
        corresponde to TF_files
    result_dir : str

    Returns
    -------
    None
    """

    if motif_fimo_files:
        for tf, motif_fimo in zip(TF_files, motif_fimo_files):
            os.system("
                bedtools intersect -c -a %s -b %s |\
                sed -e \"s/$/\t%s/g\" \
                > %s/%s" % (
                    tf, motif_fimo,
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
    """Rank all binding sites according to the ranking criteria
    All TF bs into a single file
    input: $4 spp, $5 crowdness, $6 motif_count(arbitrary), $7 TF name
    output: $1 chr, $2 summit, $3-6 same as $4-7 in input, $8 rank

    Parameters
    ----------
    rank_method : str
        one of {"RankSPP", "RankCrowdness", "RankHotness", "RankLinear_*", "RankExcludeArbitrary_*"}
    TF_files : list of str
    result_file : str
    bs_half_length : int
    rank_linear_alpha : int or None
    rank_exclude_prop : float or None

    Returns
    -------
    None
    """

    if os.path.exists(result_file):
        raise OSError("skip get_rank_file")

    if rank_method.startswith("RankSPP"):
        for tf in TF_files:
            os.system("
                awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                sort -k3,3nr  | \
                awk '{print $0\"\\t\"NR}' \
                >> %s" % (
                    bs_half_length, tf, result_file))

    elif rank_method.startswith("RankCrowdness"):
        for tf in TF_files:
            os.system("
                awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                sort -k4,4n -k3,3nr  | \
                awk '{print $0\"\\t\"NR}' \
                >> %s" % (
                    bs_half_length, tf, result_file))

    elif rank_method.startswith("RankHotness"):
        for tf in TF_files:
            os.system("
                awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                sort -k4,4nr -k3,3nr  | \
                awk '{print $0\"\\t\"NR}' \
                >> %s" % (
                    bs_half_length, tf, result_file))

    elif rank_method.startswith("RankLinear"):
        for tf in TF_files:
            os.system("
                awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                sort -k3,3nr  | awk '{print $0\"\\t\"NR}' | \
                sort -k4,4n -k3,3nr | awk '{print $0\"\\t\"NR}' | \
                awk '{print $0\"\\t\"$7+$8*%g}' | \
                sort -k9,9n | \
                awk '{print $0\"\\t\"NR}' | \
                cut -f-6,10- \
                >> %s" % (
                    bs_half_length, tf, rank_linear_alpha, result_file))

    elif rank_method.startswith("RankExcludeArbitrary"):

        tmp_file = utils.os.mkdir_tmp()

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
            os.system("
                awk '{print $1\"\\t\"$2+%d\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' %s | \
                sort -k3,3nr | \
                awk '{if($4<=%d){print $0}}' | \
                awk '{print $0\"\\t\"NR}' \
                >> %s" % (
                    bs_half_length, tf, thresh, result_file))
    else:
        raise Exception("Unknown rank method")
    return


def get_motif_enrichment(rank_file, n_seq, result_dir, enrichment_method="EnrichmentFraction"):
    """Get motif enrichment
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
    os.makedirs(result_dir, exist_ok=True)

    if os.path.exists("%s/%s_%d.txt" % (result_dir, enrichment_method, n_seq)):
        raise OSError("skip %s_%d.txt" % (enrichment_method, n_seq))

    tmp_dir = utils.os.mkdir_tmp()
    tmp_file = os.path.join(tmp_dir, "tmpfile")

    # keep all peaks with ranking < n_seq
    # then group by TF
    # then calculate motif enrichment

    enrichment = pd.DataFrame()

    os.system("awk '{if($7<=%d){print $0}}' %s > %s" % (n_seq, rank_file, tmp_file))

    df = pd.read_csv(tmp_file, sep="\t", header=None)
    shutil.rmtree(tmp_dir)

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
        raise Exception("TBD, not implemented yet")
    elif enrichment_method == "EnrichmentAvgHitPerBS":
        raise Exception("TBD, not implemented yet")
    else:
        raise Exception("Unknown enrichment method")

    df[4] = df[4].replace(0, np.nan)
    df = df.groupby(5).count()
    n_seq_result = (df[4] / df[0])

    n_seq_result.to_csv(
        "%s/%s_%d.txt" % (result_dir, enrichment_method, n_seq),
        sep="\t",
        header=False,
        float_format="%.3f")
    return


def summarize_motif_enrichment(enrichment_dir):
    """Compress and summarize enrichment results

    Parameters
    ----------
    enrichment_dir: str
        the main dir of enrichment, structure is main/SubFolderByRank/fileByNseq.txt

    Returns
    -------
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
    """Split rank file by TF name into individual BED file
    keep top 2*inference_n_seq sequences(if not enough, keep all)
    each region truncated to +- inference_nbp around summit

    Parameters
    ----------
    rank_file: str
    inference_n_seq: int
    inference_nbp: int
    result_dir: str

    Returns
    -------
    None
    """

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
                    f.write("
                        dreme -oc {infer_out_dir}/{tf}/ -k 8 -m 10 -e 1000 -p {fasta_dir}/{tf}.fasta;\
                        mv {infer_out_dir}/{tf}/dreme.txt {out_dir}_{tf}.meme;\
                        rm -rf {infer_out_dir}/{tf}/\n".format(
                            infer_out_dir=inference_motif_dir,
                            tf=tf,
                            fasta_dir=fasta_dir,
                            outdir=outdir))
                elif inference_tool == "HOMER":
                    f.write("
                        findMotifs.pl {fasta_dir}/{tf}.fasta fasta {infer_out_dir}/{tf}/ -nocheck -nogo -len 8 -noknown -basic -S 10;\
                        Rscript /home/jg2447/slayman/motif_inference/scripts/src/MoVRs_Motif2meme.R {infer_out_dir}/{tf}/homerMotifs.motifs8 {infer_out_dir}/{tf}/homerMotifs.meme;\
                        mv {infer_out_dir}/{tf}/homerMotifs.meme {out_dir}_{tf}.meme;\
                        rm -rf {infer_out_dir}/{tf}/\n".format(
                            infer_out_dir=inference_motif_dir,
                            tf=tf,
                            fasta_dir=fasta_dir,
                            outdir=outdir))
                elif inference_tool == "MEME":
                    f.write("
                        meme -oc {infer_out_dir}/{tf}/ -dna -nmotifs 10 -w 8 -maxsize 150000 -nostatus {fasta_dir}/{tf}.fasta\
                        mv {infer_out_dir}/{tf}/dreme.txt {out_dir}_{tf}.meme;\
                        rm -rf {infer_out_dir}/{tf}/\n".format(
                            infer_out_dir=inference_motif_dir,
                            tf=tf,
                            fasta_dir=fasta_dir,
                            outdir=outdir))
            else:
                if inference_tool == "DREME":
                    f.write("
                        dreme -oc {infer_out_dir}/{tf}/ -k 8 -m 10 -e 1000 -p {fasta_dir}/{tf}.fasta\n".format(
                            infer_out_dir=inference_motif_dir,
                            tf=tf,
                            fasta_dir=fasta_dir))
                elif inference_tool == "HOMER":
                    f.write("
                        findMotifs.pl {fasta_dir}/{tf}.fasta fasta {infer_out_dir}/{tf}/ -nocheck -nogo -len 8 -noknown -basic -S 10;\
                        Rscript /home/jg2447/slayman/motif_inference/scripts/src/MoVRs_Motif2meme.R {infer_out_dir}/{tf}/homerMotifs.motifs8\n".format(
                            infer_out_dir=inference_motif_dir,
                            tf=tf,
                            fasta_dir=fasta_dir))
                elif inference_tool == "MEME":
                    f.write("
                        meme -oc {infer_out_dir}/{tf}/ -dna -nmotifs 10 -w 8 -maxsize 150000 -nostatus {fasta_dir}/{tf}.fasta\n".format(
                            infer_out_dir=inference_motif_dir,
                            tf=tf,
                            fasta_dir=fasta_dir))
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
                f.write("
                    tomtom -no-ssc -oc {tomtom_out_dir}/{tf}/ -verbosity 1 -min-overlap 1 -dist pearson -evalue -thresh 10 {outdir}_{tf}.meme {known_dir}/{tf}.meme;\
                    mv {tomtom_out_dir}/{tf}/tomtom.txt {outdir2}_{tf}.tomtom;\
                    rm -rf {tomtom_out_dir}/{tf}/\n".format(
                        tomtom_out_dir=inference_tomtom_dir,
                        tf=tf,
                        outdir=outdir,
                        known_dir=known_motif_dir,
                        outdir2=outdir2))
            else:
                f.write("
                    {tomtom_out_dir}/{tf}/ -verbosity 1 -min-overlap 1 -dist pearson -evalue -thresh 10 {outdir}_{tf}.meme {known_dir}/{tf}.meme\n".format(
                        tomtom_out_dir=inference_tomtom_dir,
                        tf=tf,
                        outdir=outdir,
                        known_dir=known_motif_dir))
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


def get_inference_motif_order(motif_file):
    """
    seq2ind:
    for each motif sequence,
    return the order of that sequence in the inferred motif file (single pwm with multiple motifs)

    Parameters
    ----------
    motif_file : str

    Returns
    -------
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


def get_inference_motif_array(mfile):
    """return dict {id: motif pwm array}

    Parameters
    ----------
    mfile : str
        path to motif meme file

    Returns
    -------
    dict
        {id: motif pwm array}

    """
    motifs = []
    motif_names = []
    with open(mfile, "r") as f:
        lines = f.readlines()
        for k in range(len(lines)):
            if lines[k].startswith("letter-probability"):
                motif = np.array(lines)[k+1:k+9] # fix length 8
                motif_array = []
                for i in motif:
                    motif_array.append(i.rstrip().split(" "))
                motif = np.array(motif_array).astype(float)
                motifs.append(motif)
            if lines[k].startswith("MOTIF"):
                name = lines[k].split(" ")[1]
                motif_names.append(name)
    n2a = {}
    for i,j in zip(motifs, motif_names):
        n2a[j] = i
    return(n2a)


def get_tomtom_result(meme_file_list, tomtom_file_list, top_nn):
    """calculate accuracy, coverage, entropy of top nn inferred motif, entropy of top nn significant inferred motif.

    Parameters
    ----------
    meme_file_list : list of str
        list of paths to all inferred motifs
    tomtom_file_list : list of str
        list of paths to all inferred motifs tomtom result, order should be matched with meme_file_list
    top_nn : int
        number of top motif that should be keep

    Returns
    -------
    DataFrame
        row: TFs
        columns: accuracy, coverage, top_entropy, top_entropy_list
    """

    accuracy = []
    coverage = []
    top_entropy = []
    top_sig_entropy = []

    for meme_file, tomtom_file in zip(meme_file_list, tomtom_file_list):
        #### calculate accuracy and coverage
        tomtom_df = pd.read_csv(tomtom_file, sep="\t")
        tomtom_df = tomtom_df.sort_values('p-value').drop_duplicates("#Query ID") # for one predict motif, only use the best matched known motif
        orders = [n2i[jj] for jj in tomtom_df["#Query ID"].values] # get inferred orders for each query inferred motif
        tomtom_df_ind = [jj for jj in range(len(orders)) if orders[jj] in list(range(top_nn))] # only keep predict motif that is in the top 5
        tomtom_df = tomtom_df.iloc[tomtom_df_ind,:].copy() # all top 5 inferred motif result
        tomtom_df_sig = tomtom_df[tomtom_df['p-value'] <= 0.05].copy() # significant result only

        accuracy.append(tomtom_df_sig.shape[0])
        coverage.append(len(tomtom_df_sig["Target ID"].unique()))

        #### entropy measurement
        # get top 5 motif pwm array
        n2i = get_inference_motif_order(meme_file)
        n2a = get_inference_motif_array(meme_file)
        i2n = {}
        for j in n2i.keys():
            i2n[n2i[j]] = j
        motif_array = [n2a[i2n[j]].reshape(-1) for j in range(top_nn)]

        top_entropy_list = [motif_entropy(n2a[ii]) for ii in tomtom_df["#Query ID"].values]
        top_sig_entropy_list = [motif_entropy(n2a[ii]) for ii in tomtom_df_sig["#Query ID"].values]

        top_entropy.append(np.mean(top_entropy_list) if top_entropy_list else np.nan)
        top_sig_entropy.append(np.mean(top_sig_entropy_list) if top_sig_entropy_list else np.nan)

        # order_pt = 5 - np.array([n2i[jj] for jj in tomtom_df_sig["#Query ID"].values]) # score for significant matches, accouning order and accuracy
        # pvalue_pt = -np.log(tomtom_df_sig['p-value'].values)
        # order_pt_result = np.sum(order_pt)
        # pvalue_pt_result = np.sum(pvalue_pt)
        # order_p-value_pt_result = np.sum(order_pt * pvalue_pt)

    tfs = [os.path.basename(i).replace(".meme", "") for i in meme_file]
    result = pd.DataFrame([tfs, accuracy, coverage, top_entropy, top_sig_entropy]).T
    result.columns = ["TF", "accuracy", "coverage", "top_entropy", "top_sig_entropy"]
    return(result)


def get_tomtom_summary_data(tomtom_summary_dir, nbp_list, rank_list, tool, suffix, improve=True, tf_subset=None, measurement="accuracy"):
    """Retrieve tomtom summary data of different experiments

    Parameters
    ----------
    tomtom_summary_dir: str
    nbp_list: list of int
    rank_list: list of str
    tool: str
        {"DREME", "MEME", "HOMER"}
    suffix: str
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
                    "%s/%s_%s_%d_%s.txt" % (
                        tomtom_summary_dir, tool, rank, nn, top),
                    sep="\t")

            if tf_subset:
                data = data[data["TF"].isin(tf_subset)][measurement].values
            else:
                data = data[measurement].values

            rank_data.append(data.mean()) # mean of all TF within one rank, append for all rank
        nn_data.append(rank_data) # append for all nn

    df = pd.DataFrame(nn_data)
    df.columns = rank_list
    df.index = nbp_list

    # get improvement
    if improve:
        df = (df - df.loc[100, "RankSPP"]) / df.loc[100, "RankSPP"]
    return df
