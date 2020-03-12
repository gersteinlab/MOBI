"""
Get motifs from Kheradpour and Kellis (2014)
http://compbio.mit.edu/encode-motifs/


Required:
    matrix2meme from MEME
Input:
    Download "motifs.txt (1.1M): all the known and discovered motifs." from the website above.
    A single file containning all motifs.
Output:
	Cleaned stacked motifs
"""

import os
import shutil
from src import motif_utils
from src import utils

with open("/home/jg2447/slayman/data/motif_raw/motifs_manolis.txt", "r") as f:
    content = f.read()
    motif_list = content.split(">")
    
pwm_dir = utils.os_mkdir_tmp()
for motif in motif_list:
    if motif == "":
        continue
    name = motif.split("\n")[0].split(" ")[0]
    name = name.replace("_disc", "--disc")
    with open("%s/%s.pwm" % (pwm_dir, name), "w") as f:
        f.write(motif)
os.system("for file in %s/*.pwm; do sed 1d $file | cut -d " " -f2- | matrix2meme -dna > ~/slayman/data/motif_manolis/meme_human/$(basename $file .pwm).meme; done" % pwm_dir)
shutil.rmtree(pwm_dir)

motif_utils.stack_motif(
    sep_motif_dir="/home/jg2447/slayman/data/motif_manolis/meme_human/",
    stack_motif_dir="/home/jg2447/slayman/data/motif_manolis/meme_human_stack/")