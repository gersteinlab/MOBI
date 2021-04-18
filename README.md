# MOBI
MOtif inference with advanced BInding site selection.  

### Download
We applied our method and inferred motifs for samples we've collected from ENCODE. You could download all the motifs here.
| Sample                   | TFs  | Download Links|
|:------------------------:|:----:|:-------------:|
| *Drosophila melanogaster*| 452  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.Fly.tar.gz)|
| *Caenorhabditis elegans* | 283  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.Worm.tar.gz)|
| GM12878                  | 136  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.GM12878.tar.gz)|
| K562                     | 336  | [Download](https://github.com/gersteinlab/MOBI/tree/master/download/MOBI.K562.tar.gz)|

(Note that files will be in [MEME format](http://meme-suite.org/doc/meme-format.html). File name is the TF name but any "(", ")" and ":" characters in the TF name will be omitted, e.g. TF name A(B)C will have motif file named ABC.meme)

### Methods
![fig1](https://github.com/gaotc200/MOBI/blob/master/img/Fig1.png "Fig1")

As shown in this figure, we postulate there are several TF binding modes. In the final result of a typical ChIP-seq experiment targeting a certain TF (here shown in red), the identified binding sites could either contain the target motif (fig. A) or not (fig. B). The ideal binding sites for motif inference should be those from scenario one.  

Notice that among all these cases:
- Binding site in scenario one is least "crowded", i.e. there are fewest possibly interacting proteins in the binding site window. 
- The target binding motif should locate in or close to the peak summit (in this figure, it is the center).  

Therefore, we select binding sites with low "crowdness" score and trim the binding sites to a shorter length. Using such binding sites result in more accuate motif inference.

For more details, see our paper:  
*Jinrui Xu, Jiahao Gao, Mark Gerstein. (2020) Discovering less-is-more effects to select transcription factor binding sites informative for motif inference. (in preparation)*

### Contacts
For any questions, please contact Jiahao Gao(jiahao.gao@yale.edu)  
[Gerstein Lab](http://www.gersteinlab.org/) 2020
