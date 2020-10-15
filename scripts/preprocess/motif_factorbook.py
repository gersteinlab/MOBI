"""
Get motifs from factorbook. Then apply same preprocess as in motif_cisbp.py

Required:
    matrix2meme from MEME
Input:
Output:
	Cleaned stacked motifs
"""

import os
import re
import shutil
import pandas as pd
from bs4 import BeautifulSoup
from urllib.request import urlopen
import sys
sys.path.append("/home/jg2447/slayman/motif_inference/MOBI/scripts")

from src import motif_files


html = urlopen("http://www.factorbook.org/human/chipseq/tf/")
soup = BeautifulSoup(html, features='lxml')
t = soup.find_all('table', {'id': 'factorHeaderTable'})
tf_list = pd.read_html(str(t))[0].columns.values[0]

for tf in tf_list:
    html = urlopen("http://www.factorbook.org/human/chipseq/tf/%s" % tf)
    soup = BeautifulSoup(html, features='lxml')
    a_list = soup.find_all('a', {"href": re.compile(".*meme.html")})

    mid = 0
    for a in a_list:
        meme_url = a['href']
        html = urlopen("http://www.factorbook.org/"+meme_url)
        soup = BeautifulSoup(html, features='lxml')
        html_str = soup.text
        pwm_strs = re.findall(r'\"pwm\": \[((.|\n)*?)\"sites\"', html_str)
        if len(pwm_strs) > 0:
            for pwm in pwm_strs:
                pwm_str = pwm[0]
                pwm_str = pwm_str.replace("\n", "").replace(" ", "").replace("],[", "],\n[").replace("[", "").replace("],", "").replace("]", "").replace(",", " ")
                with open("./tmp_motif.txt", "w") as f:
                    f.write(pwm_str)
                os.system("matrix2meme -dna < ./tmp_motif.txt > /home/jg2447/slayman/data/motif/factorbook/meme_human/%s--%d.meme" % (tf, mid))
                mid = mid + 1

motif_files.stack_motif(
    sep_motif_dir="/home/jg2447/slayman/data/motif/factorbook/meme_human/",
    stack_motif_dir="/home/jg2447/slayman/data/motif/factorbook/meme_human_stack/")
