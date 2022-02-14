import os
import pandas as pd
import numpy as np
from mobi import prepare_input, calculate_crowdness, format_length, ranking, prepare_fasta, inference, tomtom

#### user input start
data_chip   = "example_data/ChIP/"
data_meta = "example_data/TF_names.txt"
genome_fasta = "example_data/dm6.fa"
result_main = "example_result/"

width_list = [100]
rank_list = ["RankLinear1.0"]
tool_list = ["DREME"]
#### user input end


df = pd.read_csv(data_meta, sep="\t", header=None)
TF_files = df[0].values
TF_names = df[1].values
genome_size = pd.read_csv(genome_fasta+".fai", sep="\t", header=None)
genome_size = pd.Series(genome_size[1].values, index=genome_size[0]).to_dict()

mobi_data_dir = os.path.join(result_main, "input_bed/")
all_summit_file = os.path.join(result_main, "allSummit.txt")
crowdness_dir = os.path.join(result_main, "crowdness/")
rank_file_dir = os.path.join(result_main, "rank_file/")
fasta_dir = os.path.join(result_main, "fasta/")
inference_dir_raw = os.path.join(result_main, "inference_raw/")
inference_dir = os.path.join(result_main, "inference/")
tomtom_dir_raw = os.path.join(result_main, "tomtom_raw/")
tomtom_dir = os.path.join(result_main, "tomtom/")


#### prepare bed file to calculate C-score

all_TF_files = [ii.split(".")[0] for ii in os.listdir(data_chip) if not ii.startswith(".")]
os.makedirs(mobi_data_dir, exist_ok=True)

for ii in range(len(all_TF_files)):
    prepare_input.bedfile_get_subset(
        infile=os.path.join(data_chip, all_TF_files[ii]+".bed.gz"),
        outfile=os.path.join(mobi_data_dir, all_TF_files[ii]+".bed"),
        data_proportion=3000, # number of sites to use for C-score calculation
        sort_column=7, # SPP signal column
        read_method="zcat") # cat/zcat for .bed/.bed.gz


#### calculate C-score and assign C-score to binding regions

os.makedirs(crowdness_dir, exist_ok=True)
calculate_crowdness.get_crowdness(
    input_files=[os.path.join(mobi_data_dir, ii+".bed") for ii in all_TF_files],
    output_dir=crowdness_dir,
    crowdness_window=250,
    genome_size=genome_size,
    subset=[os.path.join(mobi_data_dir, ii+".bed") for ii in TF_files],
    all_summit_file_path=all_summit_file)

for ii, jj in zip(TF_files, TF_names):
    os.rename(os.path.join(crowdness_dir, ii+".bed"), os.path.join(crowdness_dir, jj+".bed")) # rename to TF names

for width in width_list:
    os.makedirs(os.path.join(crowdness_dir, str(width)), exist_ok=True)    
    for ii in TF_names:
        df = pd.read_csv(os.path.join(crowdness_dir, ii+".bed"), sep="\t", header=None)
        df[5] = 0 # add a column of 0. This column is not used in this example, but store the motif occurrence in other analysis
        start = (df[1]+df[2])//2 - width
        end = (df[1]+df[2])//2 + width
        df[1] = start
        df[2] = end
        df.to_csv(os.path.join(crowdness_dir, str(width), ii+".bed"), sep="\t", header=False, index=False)


#### rank all the binding sites based on the given ranking method

os.makedirs(rank_file_dir, exist_ok=True)
for width in width_list:
    for rank in rank_list:
        for ii in TF_names:
            ranking.get_rank_file(
                infile=os.path.join(crowdness_dir, ii+".bed"),
                outfile=os.path.join(rank_file_dir, rank+"_"+str(width)+"_"+ii+".bed"),
                rank=rank)

for width in width_list:
    for rank in rank_list:
        all_rank_file = os.path.join(rank_file_dir, rank+"_"+str(width)+"_*")
        final_rank_file = os.path.join(rank_file_dir, rank+"_"+str(width)+".bed")        
        os.system("cat %s > %s" % (all_rank_file, final_rank_file))
        os.system("rm %s" % all_rank_file)


#### retrieve sequence fasta of the selected binding regions

os.makedirs(fasta_dir, exist_ok=True)
for width in width_list:
    for rank in rank_list:
        rank_file = os.path.join(rank_file_dir, rank+"_"+str(width)+".bed") 
        output_dir = os.path.join(fasta_dir, rank+"_"+str(width))
        os.makedirs(output_dir, exist_ok=True)
        
        prepare_fasta.get_fasta(
            rank_file=rank_file,
            genome_fasta=genome_fasta,
            out_dir=output_dir, 
            allow_n_seq=250)


#### generate command lists of inference

os.makedirs(inference_dir_raw, exist_ok=True)
os.makedirs(inference_dir, exist_ok=True)
ff = open("./joblist_inference.sh", "w")
for tool in tool_list:
    os.makedirs(os.path.join(inference_dir_raw, tool), exist_ok=True)
    os.makedirs(os.path.join(inference_dir, tool), exist_ok=True)
    for width in width_list:
        for rank in rank_list:
            for motif in TF_names:
                cmd = inference.inference_cmd(
                    tool=tool,
                    out_dir=os.path.join(inference_dir_raw, tool, "%s_%s_%s/" % (rank, width, motif)),
                    out_file=os.path.join(inference_dir, tool, "%s_%s_%s.meme" % (rank, width, motif)),
                    fasta_file=os.path.join(fasta_dir, "%s_%s" % (rank, width), "%s.fasta" % motif),
                    homer_rscript="MoVRs_Motif2meme.R")
                ff.write(cmd)
ff.close()

