import os
import tempfile
import pandas as pd
import pybedtools


def bedfile_get_window_from_mid(infile, outfile, width, genome_size, keep_col=[3,4]):
    df = pd.read_csv(infile, sep="\t", header=None)
    newdf = pd.DataFrame([])
    newdf[0] = df[0]
    newdf[1] = (df[1]+df[2])//2-width
    newdf[2] = (df[1]+df[2])//2+width
    newdf.loc[newdf[1]<0, 1] = 0 # negative start will change to 0
    newdf[3] = newdf[0].map(genome_size)
    newdf.loc[newdf[2]>newdf[3], 2] = newdf.loc[newdf[2]>newdf[3], 3] # for any end > genome size, end change to genome size
    del newdf[3]
    newdf = pd.concat([newdf, df.loc[:, keep_col]], axis=1)
    newdf[1] = newdf[1].astype(int)
    newdf[2] = newdf[2].astype(int)
    newdf.to_csv(outfile, sep="\t", header=False, index=False)


def bedfile_get_window_with_motif_count(in_file, in_motif_file, out_dir, width, genome_size, place_holder=False):
    # if place holder, then just add 0 to the last column
    # else intersect with FIMO result to get the number of motifs in each binding sites
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        bedfile_get_window_from_mid(
            in_file,
            os.path.join(tmp_dir, os.path.basename(in_file)),
            width,
            genome_size,
            keep_col=[3,4])

        if place_holder:
            df = pd.read_csv(os.path.join(tmp_dir, os.path.basename(in_file)), sep="\t", header=None)
            df[5] = 0
            df.to_csv(os.path.join(out_dir, os.path.basename(in_file)), sep="\t", header=False, index=False)
        else:
            bt = pybedtools.BedTool(os.path.join(tmp_dir, os.path.basename(in_file)))
            fimo_bt = pybedtools.BedTool(in_motif_file)
            fimo_result = bt.intersect(fimo_bt, c=True)
            fimo_result.saveas(os.path.join(out_dir, os.path.basename(in_file)))