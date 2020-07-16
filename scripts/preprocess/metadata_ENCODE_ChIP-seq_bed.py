"""
Filter for ENCODE samples

Input:
	Metadata file from ENCODE
Output:
	Clean metadata for certain bio sample, e.g. worm/fly/GM12878/K562
"""

from src import encode_meta

# add TF to metafiles
df = encode_meta.add_TF(
    encode_file="/home/jg2447/slayman/data/ENCODE_ChIP_seq/metadata/ENCODE_TF-ChIP-seq_human.tsv",
    outfile="/home/jg2447/slayman/data/ENCODE_ChIP_seq/metadata/ENCODE_TF-ChIP-seq_human_TFadded.tsv",
    species="human")

# subset metafiles for each dataset (human all cell line/GM12878/K562/worm/fly and etc)
df = encode_meta.filter_subset(
    encode_file="/home/jg2447/slayman/data/ENCODE_ChIP_seq_metadata/ENCODE_TF-ChIP-seq_human_TFadded.tsv",
    outfile=None,
    Assembly="hg19"
    Biosample_type="cell line",
    File_format="bed narrowPeak",
    Output_type="optimal idr thresholded peaks")
    # if filtering for specific cell line, use Biosample_term_name = "GM12878"

# add motif info
meta = encodeMeta.add_motif(
    encode_file=df,
    motif_dir="/home/jg2447/slayman/data/FIMO/human_stack/FIMOp_0.000100_chr/",
    outfile="/home/jg2447/slayman/motif_inference/result/metadata/human-All.txt")

# download the ChIP-seq data
encodeMeta.download_file(
    encode_file=meta,
    TF_dir="/home/jg2447/slayman/data/ENCODE_ChIP_seq/bed/human_hg19/",
    ext=".bed.gz")
