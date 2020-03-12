import os
import shutil
from src import utils
from src import narrowPeak
from src import motif_utils

fimo_output_dir = "/home/jg2447/scratch60/fimo_output_human/"
main_result_dir = "/home/jg2447/slayman/data/FIMO/human_stack/"


human_genome_name = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY", "chrM"]
fly_genome_name = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM"]
worm_genome_name = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"]


# reformat FIMO output to bed format
motif_utils.reformat_FIMO_result_to_bed(
    fimo_output_dir=fimo_output_dir,
    output_dir="%s/FIMO_output/" % main_result_dir)
shutil.rmtree(fimo_output_dir)

# filter to get significant hit
motif_utils.filter_FIMO_result(
    input_dir="%s/FIMO_output/" % main_result_dir, 
    output_dir="%s/FIMOp_0.000100/" % main_result_dir, 
    fimo_p=0.0001,
    fimo_q=None)

# keep only regions in chromosome, discard scaffolds
utils.os_mkdir_empty("%s/FIMOp_0.000100_chr/" % main_result_dir)
for i in os.listdir("%s/FIMOp_0.000100/" % main_result_dir):
    narrowPeak.filter_chr_name(
        "%s/FIMOp_0.000100/%s" % (main_result_dir, i),
        "%s/FIMOp_0.000100_chr/%s" % (main_result_dir, i),
        human_genome_name)
shutil.rmtree("%s/FIMOp_0.000100/" % main_result_dir)

# remove files with no significant hit regions
utils.os_remove_size_zero_file("%s/FIMOp_0.000100_chr/" % main_result_dir)

"""
Also run for worm and fly,
change path name and human_genome_name
"""