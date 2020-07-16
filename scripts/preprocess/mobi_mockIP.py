import os
import pandas as pd
from src import utils

# modify file format as the same as directly download from ENCODE
eids = [i.split("_sites.bed")[0] for i in os.listdir("/gpfs/project/fas/gerstein/jx98/projects/chipSeq/script_07072018/publishData/peakSite/ChIP-seq_wMockIP/FruitFly/") if i.endswith("sites.bed")]
df = pd.read_csv("/home/jg2447/slayman/motif_inference/result/metadata/fly-All.txt", sep="\t")
df = df[df['Experiment accession'].isin(eids)].drop_duplicates("Experiment accession")
df.to_csv("/home/jg2447/slayman/motif_inference/result/metadata/fly-All_mockIP.txt", sep="\t", index=False)

eid2fid = utils.df.colToDict(df, "Experiment accession", "File accession")

os.makedirs(para_chip_data_dir, exist_ok=True)

for ee in eids:
    if ee in eid2fid.keys():
        shutil.copyfile("/gpfs/project/fas/gerstein/jx98/projects/chipSeq/script_07072018/publishData/peakSite/ChIP-seq_wMockIP/FruitFly/"+ee+"_sites.bed", para_chip_data_dir+eid2fid[ee]+".bed")
