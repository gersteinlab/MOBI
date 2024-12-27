import os
import tempfile
import pandas as pd
import pybedtools


def bedfile_get_window(infile, outfile, width, genome_size, keep_col=[3,4]):
    df = pd.read_csv(infile, sep="\t", header=None)
    newdf = pd.DataFrame([])
    newdf[0] = df[0]
    newdf[1] = (df[1]+df[2])//2-width
    newdf[2] = (df[1]+df[2])//2+width
    newdf.loc[newdf[1]<0, 1] = 0
    newdf[3] = newdf[0].map(genome_size)
    newdf.loc[newdf[2]>newdf[3], 2] = newdf.loc[newdf[2]>newdf[3], 3]
    del newdf[3]
    newdf = pd.concat([newdf, df.loc[:, keep_col]], axis=1)
    newdf[1] = newdf[1].astype(int)
    newdf[2] = newdf[2].astype(int)
    newdf.to_csv(outfile, sep="\t", header=False, index=False)


def bedfile_get_window_with_motif_count(in_file, in_motif_file, out_dir, width, genome_size):
    with tempfile.TemporaryDirectory() as tmp_dir:
        bedfile_get_window(
            in_file,
            os.path.join(tmp_dir, os.path.basename(in_file)),
            width,
            genome_size,
            keep_col=[3,4])

        bt = pybedtools.BedTool(os.path.join(tmp_dir, os.path.basename(in_file)))
        fimo_bt = pybedtools.BedTool(in_motif_file)

        fimo_result = bt.intersect(fimo_bt, c=True)
        fimo_result.saveas(os.path.join(out_dir, os.path.basename(in_file)))