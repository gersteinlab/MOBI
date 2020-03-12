from src import utils
from src import motif_utils

# run FIMO with slurm
motif_utils.run_FIMO(
    motif_dir="/home/jg2447/slayman/data/motif_cisbp/meme_human_stack/",
    genome_file="/home/jg2447/slayman/data/genome/GRCh37.p13.genome.fa",
    fimo_output_dir="/home/jg2447/scratch60/fimo_output_human/",
    joblist_file="/home/jg2447/scratch60/fimo_joblist_human.txt",
    bgfreq=False,
    subset=None)
utils.dsq_get_sh("/home/jg2447/scratch60/fimo_joblist_human.txt")