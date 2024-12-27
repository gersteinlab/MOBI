import os
import gzip
import pandas as pd

def bedfile_get_subset(infile, outfile, data_proportion, sort_column=7, read_method="zcat"):
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
    sort_column: int
    read_method: str
        os function to read the input file, e.g zcat, cat.

    Returns
    ----------------
    int from output of os.system
    """

    # get total lines of input file
    if read_method == "zcat":
        size = sum(1 for line in gzip.open(infile))
    else:
        size = sum(1 for line in open(infile))

    # number of sites to keep
    if data_proportion > 1:
        nsize = data_proportion
    else:
        nsize = np.floor(size * data_proportion).astype(int)

    # subset and sort
    subset_cmd = "sed %dq" % nsize
    if not sort_column:
        sort_cmd = "sort -k1,1 -k2,2n"
    else:
        sort_cmd = "sort -k%d,%dnr" % (sort_column, sort_column)
    p = os.system("%s < %s | %s | %s > %s" % (read_method, infile, sort_cmd, subset_cmd, outfile))
    return(p)