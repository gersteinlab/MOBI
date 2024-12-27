import os
import gzip
import pandas as pd

def bedfile_get_subset(infile, outfile, data_proportion, sort_column=7):
    """Get a subset(could be all) of sites from the given file.
    Input will be sorted (by SPP score, which is col 7, if the original file is from ENCODE ChIP-seq)
    Only the top sites will be kept
    Output files in (unzip) bed format, sorted

    Parameters
    ----------------
    infile: str
        standard bed file in .bed or .bed.gz, col 7 should be score
    outfile: str
    data_proportion: float or int
        number of bs to keep. If >1 then indicate the number. If 0-1, then is the proportion
    sort_column: int, it sort in the ascending=False way
    read_method: str
        os function to read the input file, e.g zcat, cat.

    Returns
    ----------------
    int from output of os.system
    """

    df = pd.read_csv(infile, sep="\t", header=None)
    size = df.shape[0]

    # number of regions to keep
    if data_proportion > 1:
        nsize = data_proportion
    else:
        nsize = np.floor(size * data_proportion).astype(int)

    if not sort_column:
        df = df.sort_values([0,1,2]).iloc[:nsize,] # default, sort by increasing coordinate order
    else:
        df = df.sort_values(sort_column, ascending=False).iloc[:nsize,] # sort signal from large to small

    df.to_csv(outfile, sep="\t", header=False, index=False)
    return