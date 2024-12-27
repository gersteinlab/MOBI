import os
from mobi import tomtom

####
width_list = [100]
rank_list = ["RankLinear1.0"]
tool_list = ["DREME"]
data_meta = "example_data/TF_names.txt"
result_main = "example_result/"
known_motif_dir   = "example_data/motifs/"
topN = 5
p_threshold = 0.05
####

df = pd.read_csv(data_meta, sep="\t", header=None)
TF_files = df[0].values
TF_names = df[1].values

inference_dir = os.path.join(result_main, "inference/")
tomtom_dir_raw = os.path.join(result_main, "tomtom_raw/")
tomtom_dir = os.path.join(result_main, "tomtom/")
summary_dir = os.path.join(result_main, "summary/")
stat_dir = os.path.join(result_main, "stats/")



## get summary report of tomtom results
## this include pre define measurements, including Accuracy/Coverage, median/mean logp, median/mean motif entropy, and logp/entropy for top results regardless of significance

os.makedirs(summary_dir, exist_ok=True)

for tool in tool_list:
    for width in width_list:
        for rank in rank_list:
            
            measurements_summary = []
            for motif in TF_names:
                meme_file = os.path.join(inference_dir, tool, "%s_%s_%s.meme" % (rank, width, motif))
                tomtom_file = os.path.join(tomtom_dir, tool, "%s_%s_%s.tomtom" % (rank, width, motif))
                measurements = tomtom.evaluate_statistics(
                    meme_file=meme_file,
                    tomtom_file=tomtom_file,
                    topN=topN, 
                    p_threshold=p_threshold)
                measurements_summary.append(measurements)
                
            result = pd.DataFrame(measurements_summary)
            result.index = TF_names
            result = result.reset_index()
            result.columns = ["TF", "Accuracy", "Coverage", "median_logp", "mean_logp", "median_entropy", "mean_entropy", "median_logp_top", "mean_logp_top", "median_entropy_top", "mean_entropy_top"]
            result.to_csv(os.path.join(summary_dir, "%s_%s_%s.txt" % (tool, rank, width)), sep="\t", index=False)
            
os.makedirs(stat_dir, exist_ok=True)

for tool in tool_list:
    width_result = []
    for width in width_list:
        rank_result = []
        for rank in rank_list:
            df = pd.read_csv(os.path.join(summary_dir, "%s_%s_%s.txt" % (tool, rank, width)), sep="\t")
            rank_result.append(df['Accuracy'].mean())
        width_result.append(rank_result)
    result = pd.DataFrame(width_result)

    spp_df = pd.read_csv(os.path.join(summary_dir, "%s_RankSPP_100.txt" % (tool)), sep="\t")
    spp = spp_df['Accuracy'].mean()
    result = result.T
    result.index = rank_list
    result.columns = width_list
    result_relative = (result - spp)/spp
    result.to_csv(os.path.join(stat_dir, "%s_accuracy.txt" % tool), sep="\t", float_format="%.4f")
    max_idx = np.unravel_index(result_relative.values.argmax(), [len(rank_list), len(width_list)])

    with open(os.path.join(stat_dir, "%s_idx.txt" % tool), "w") as f:
        f.write(rank_list[max_idx[0]]+"\t"+str(width_list[max_idx[1]]))
