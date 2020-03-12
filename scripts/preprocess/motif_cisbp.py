"""
Filter experiment derived motifs from CIS-BP, convert to meme format, stack motifs

Input:
	TF_information.txt file from cisbp
	Raw motif pwm from cisbp
Output:
	Cleaned stacked motifs
"""

import os
import shutil
from src import motif_utils
from src import utils

# get cisbp df
cisbp_df = motif_utils.cisbp_meta_filter_in_vitro(
    "/home/jg2447/slayman/data/motif_raw/CISBP_Homo_sapiens/TF_Information.txt", 
    outfile=None)

# pwm to meme: cp pwm to tmp dir, format convert, remove tmp dir
tmp_dir = utils.os_mkdir_tmp()
for i in range(cisbp_df.shape[0]):
    shutil.copy2(
        "/home/jg2447/slayman/data/motif_raw/CISBP_Homo_sapiens/pwms_all_motifs/%s.txt" % cisbp_df.iloc[i, 1], 
         "%s/%s--%s.pwm" % (tmp_dir, cisbp_df.iloc[i, 0], cisbp_df.iloc[i, 1]))
os.system("src/pwm2meme.sh %s /home/jg2447/slayman/data/motif_cisbp/meme_human" % tmp_dir)
shutil.rmtree(tmp_dir)

# stack motif
motif_utils.stack_motif(
    sep_motif_dir="/home/jg2447/slayman/data/motif_cisbp/meme_human/",
    stack_motif_dir="/home/jg2447/slayman/data/motif_cisbp/meme_human_stack/")


"""
Also run for worm and fly
Need to change the corresponding folder name
"""