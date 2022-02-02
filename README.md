# MOBI
MOtif inference with advanced BInding site selection.  

### Requirement
To run the main script:
⋅⋅* python3
⋅⋅* pybedtools
⋅⋅* numpy
⋅⋅* pandas
To run the inference step, you need to have DREME/MEME/STREME/HOMER properly installed.  
You can test this by `meme -version` or other equivalent commands.  

### Run the example
1. Download the github and unzip into `MOBI/`, go to the folder `cd MOBI/`
2. Download the example ChIP-seq data: `bash DownloadExampleData.sh`
3. Download the genome https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz into `example_data/`, then uncompress and index the fasta (require samtools): ```cd example_data/
gunzip dm6.fa.gz
samtools faidx dm6.fa
cd ..```
4. Run the scripts `python MOBI.py`. This will generate a file called `joblist_inference.sh`.
5. Run the scripts `bash joblist_inference.sh`. All predicted motifs are now under `example_result/inference/`.

### Run the scripts with your own data
1. Modify line 7-14 in file MOBI.py accordingly (will be updated to argument).
- `data_chip`: str. Path to the folder containing all the ENCODE ChIP-seq files. All these files are used to calculated the C-score (see below)
- `data_meta`: str. Path to the tab-seperated file. The first column is the basename of the file and the second column is the TF name. Notice the first column should be a subset of the basenames of files in `data_chip`
- `genome_fasta`: str. Path to the genome fasta file. The index file should be generated beforehand and located in the same folder
- `result_main`: str. Path to the main result folder.
- `width_list`: list of int. A list of binding regions length to try. 100 indicates a binding regions of +-100bp around ChIP-seq peak summit (a total of 200bp)
- `rank_list`: list of str. Choice of `RankSPP`, `RankCrowdness` or `RankLinear0.1` where 0.1 could be changed to other float. For detail of the ranking method, see the paper
- `tool_list`: list of str. Choice of `DREME`, `MEME`, `STREME` and `HOMER`.


### Download
We applied MOBI to the 4 samples from ENCODE respectively with their best parameters and infer the motifs. You could download all the motifs here.
| Sample                   | TFs  | Download Links|
|:------------------------:|:----:|:-------------:|
| *Drosophila melanogaster*| 452  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.Fly.tar.gz)|
| *Caenorhabditis elegans* | 283  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.Worm.tar.gz)|
| GM12878                  | 136  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.GM12878.tar.gz)|
| K562                     | 336  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.K562.tar.gz)|

(Note that files will be in [MEME format](http://meme-suite.org/doc/meme-format.html). File name is the TF name but any "(", ")" and ":" characters in the TF name will be omitted, e.g. TF name A(B)C will have motif file named ABC.meme)

### Methods overview
![fig1](https://github.com/gaotc200/MOBI/blob/master/img/Fig1.png "Fig1")

As shown in this figure, we postulate there are several TF binding modes. In the result of a typical ChIP-seq experiment targeting a certain TF (here shown in red), the identified binding sites could either contain the target motif (fig. A) or not (fig. B). The ideal binding sites for motif inference should be those from scenario one.  

Notice that among all these cases:
- Binding site in scenario one is least "crowded", i.e. there are fewest possibly interacting proteins in the binding site window. 
- The target binding motif should locate in or close to the peak summit (in this figure, it is the center).  

Therefore, we select binding sites with low "crowdness" score and trim the binding sites to a shorter length. Using such binding sites result in more accuate motif inference.

For more details, see our paper:  
*Jinrui Xu, Jiahao Gao, Mark Gerstein. (2020) Discovering less-is-more effects to select transcription factor binding sites informative for motif inference. (in preparation)*

### Contacts
For any questions, please contact Jiahao Gao(jiahao.gao@yale.edu)  
[Gerstein Lab](http://www.gersteinlab.org/) 2022
