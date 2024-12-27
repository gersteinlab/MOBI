# MOBI
MOtif inference with advanced BInding site selection.  

### Methods overview
![fig1](https://github.com/gersteinlab/MOBI/blob/main/img/fig1.png "Fig1")

As shown in this figure, we postulate there are several TF binding modes. In the result of a typical ChIP-seq experiment targeting a certain TF (here shown in red), the identified binding sites could either contain the target motif (fig. A) or not (fig. B). The ideal binding sites for motif inference should be those from scenario one.  

Notice that among all these cases:
- Binding site in scenario one is least "crowded", i.e. there are fewest possibly interacting proteins in the binding site window. 
- The target binding motif should locate in or close to the peak summit (in this figure, it is the center).  

Therefore, we select binding sites with low "crowdness" score and trim the binding sites to a shorter length. Using such binding sites result in more accuate motif inference.

For more details, see our paper:  
[Xu, J., Gao, J., Ni, P., & Gerstein, M. (2024). Less-is-more: selecting transcription factor binding regions informative for motif inference. Nucleic Acids Research, gkad1240.](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad1240/7517490)



### Download
We applied MOBI to the 4 samples from ENCODE respectively with their best parameters and infer the motifs. You could download all our predicted motifs here.
|Sample| TFs | DREME | MEME | STREME | HOMER | DESSO |
|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| *Drosophila melanogaster* | 454 | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Fly_DREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Fly_MEME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Fly_STREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Fly_HOMER.md) | [52 TFs](https://github.com/gersteinlab/MOBI/blob/main/download/Fly_DESSO.md) |
| *Caenorhabditis elegans* | 283 | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Worm_DREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Worm_MEME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Worm_STREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/Worm_HOMER.md) | [33 TFs](https://github.com/gersteinlab/MOBI/blob/main/download/Worm_DESSO.md) |
| GM12878 | 136 | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/GM12878_DREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/GM12878_MEME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/GM12878_STREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/GM12878_HOMER.md) | [56 TFs](https://github.com/gersteinlab/MOBI/blob/main/download/GM12878_DESSO.md) |
| K562 | 336 | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/K562_DREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/K562_MEME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/K562_STREME.md) | [Download](https://github.com/gersteinlab/MOBI/blob/main/download/K562_HOMER.md) | [89 TFs](https://github.com/gersteinlab/MOBI/blob/main/download/K562_DESSO.md) |

(Note that files will be in [MEME format](http://meme-suite.org/doc/meme-format.html). File name is the TF name but any "(", ")" and ":" characters in the TF name will be omitted, e.g. TF name A(B)C will have motif file named ABC.meme)

### Examples output motifs
See [here](http://archive2.gersteinlab.org/proj/MOBI/fly.html) for a few examples of the inferred motifs.

### Download C-score files
We provide the C-score we used in the paper in bigWig format for users to download. (If the download didn't start, right click to copy the link and paste into the browser manually).  
[Fly](http://archive2.gersteinlab.org/proj/MOBI/bw/fly_C-score.bw)  
[Worm](http://archive2.gersteinlab.org/proj/MOBI/bw/worm_C-score.bw)  
[GM12878](http://archive2.gersteinlab.org/proj/MOBI/bw/GM12878_C-score.bw)  
[K562](http://archive2.gersteinlab.org/proj/MOBI/bw/K562_C-score.bw)


### Requirement of scripts
- python3
- pybedtools
- numpy
- pandas

To run the inference step, **you need to have DREME/MEME/STREME/HOMER properly installed**. See [here](https://meme-suite.org/meme/doc/download.html) for meme-suite installation. You can test by `meme -version` or other equivalent commands.  

### Run the example
1. Download the github and unzip into `MOBI/`, go to the folder `cd MOBI/`.
2. Download the example ChIP-seq data: `bash DownloadExampleData.sh`.
3. Download the genome https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz into `example_data/`, then uncompress and index the fasta (require samtools): 
```
cd example_data/  
gunzip dm6.fa.gz  
samtools faidx dm6.fa  
cd ..
```
4. Run the scripts `python MOBI.py`. This will generate a file called `joblist_inference.sh`. All intermediate files are in `example_result/`.
5. Make sure you have the inference tool (e.g. DREME) installed. Run the scripts `bash joblist_inference.sh`. All predicted motifs are now under `example_result/inference/`.

### Infer motifs for your own data
Modify line 7-14 in file MOBI.py accordingly (will be updated to argument). Then run step 4 and 5 in the previous section.
- `data_chip`(*str*): Path to the folder containing all the ENCODE ChIP-seq files. All these files are used to calculated the C-score (see below)
- `data_meta`(*str*): Path to a tab-seperated file. The first column is the basename of the file and the second column is the TF name. Notice the first column should be a subset of the basenames of files in `data_chip`
- `genome_fasta`(*str*): Path to the genome fasta file. The index file should be generated beforehand and located in the same folder
- `result_main`(*str*): Path to the main result folder.
- `width_list`(*list* of *int*): A list of binding regions length to try. 100 indicates a binding regions of +-100bp around ChIP-seq peak summit (a total of 200bp)
- `rank_list`(*list* of *str*): Choice of `RankSPP`, `RankCrowdness` or `RankLinear0.1` where 0.1 could be changed to other float. For detail of the ranking method, see the paper
- `tool_list`(*list* of *str*): Choice of `DREME`, `MEME`, `STREME` and `HOMER`.


### Finding optimal parameters for your own data
In order to find the best ranking method (weights) and binding sites width, we simply do a brute force search for the parameters that give the best result:
- Run the above section with `width_list` and `rank_list` cover all the parameters you want to try.
- Having a list of "known" motifs in the MEME format. You could download this from [Cis-BP](http://cisbp.ccbr.utoronto.ca/) or other database.
- Modify line 5-12 in `MOBI_tomtom.py`, run the scripts with `python MOBI_tomtom.py`. This will generate a file called `joblist_tomtom.sh`. Make sure you have meme-suite install by verifing `tomtom --help`. Run the script with `bash joblist_tomtom.sh`. This is to compare the inferred motifs to the known motifs.
- Modify line 5-12 in `MOBI_stats.py`, run the scripts with `python MOBI_stats.py`. The best parameters will be shown in `example_result/stats/DREME_idx.txt` if you are using DREME. You can find the result for this optimal paremeters by the file names in `example_result/inference/DREME/`.



### Contacts
For any questions, please contact Jiahao Gao(jiahao.gao@yale.edu)  
[Gerstein Lab](http://www.gersteinlab.org/) 2024  
MIT License
