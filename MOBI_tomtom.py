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


os.makedirs(tomtom_dir_raw, exist_ok=True)
os.makedirs(tomtom_dir, exist_ok=True)

ff = open("./joblist_tomtom", "w")
for tool in tool_list:
    os.makedirs(tomtom_dir_raw+"/"+tool, exist_ok=True)
    os.makedirs(tomtom_dir+"/"+tool, exist_ok=True)
    for width in width_list:
        for rank in rank_list:
            for motif in TF_names:
                cmd = tomtom.tomtom_cmd(
                    out_dir=os.path.join(tomtom_dir_raw, tool, "%s_%s_%s/" % (rank, width, motif)),
                    out_file=os.path.join(tomtom_dir, tool, "%s_%s_%s.tomtom" % (rank, width, motif)),
                    infer_meme=os.path.join(inference_dir, tool, "%s_%s_%s.meme" % (rank, width, motif)),
                    known_meme=os.path.join(known_motif_dir, motif+".meme"))
                ff.write(cmd)
ff.close()
