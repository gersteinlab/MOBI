import os
import sys
sys.path.append("/home/jg2447/slayman/motif_inference/MOBI/scripts")

from src import motif_fimo

# run FIMO with slurm
motif_fimo.run_FIMO(
    motif_dir="/home/jg2447/slayman/data/motif/cisbp/meme_human_stack/",
    genome_file="/home/jg2447/slayman/data/genome/GRCh37.p13.genome.fa",
    fimo_output_dir="/home/jg2447/scratch60/fimo_output_human/",
    joblist_file="/home/jg2447/scratch60/fimo_joblist_human.txt",
    bgfreq=False,
    subset=None)
os.system("dSQ --job-file /home/jg2447/scratch60/fimo_joblist_human.txt")
