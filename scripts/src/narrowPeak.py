"""
Utils functions for ENCODE narrow peak bed format
Could be applied to other BED format
"""

import os
import gzip
import numpy as np

__all__ = ["filter_chr_name",
           "get_subset",
           "get_fix_width_region",
           "rm_column_negative"]

def filter_chr_name(infile, outfile, chr_names):
    """
    Filter records in the given chromosomes

    Parameters
    ----------------
    input_dir: str
        contain a list of .bed file. Assume first column representing the chromosome
    output_dir: str
        same name.bed as in input_dir
    chr_names: list of str
        list of chromosome names to keep

    Returns
    ----------------
    None
    """

    with open(infile, "r") as fi:
        with open(outfile, "w") as fo:
            for line in fi:
                if line.split("\t")[0] in chr_names:
                    fo.write(line)
    return


def get_subset(infile, outfile, data_proportion, sort_column=7, read_method="zcat"):
    """
    Get a subset(could be all) of sites from the given file. 
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


def get_fix_width_region(infile, outfile, width, summit_col=10, keep_col=[4,5,6,7,8,9,10], read_method="cat"):
    """
    Truncate each binding site to a fixed length or even only the summit.
    Outfile format: chr, start, end, kept_column

    Parameters
    ----------------
    infile: str
        bed, bed.gz file
    outfile: str
    width: int
        output binding site width, +- nbp around summit.
        If 0 or negative, only 1-base (a.k.a summit) is kept
    summit_col: int or None
        the nth column in infile indicating the summit. If None, use midpoint of start and end instead.
    keep_col: int, list of int or None
        the nth column(s) to keep in the output
    read_method: str
        os function to read the input file, e.g zcat, cat

    Returns
    ----------------
    int: return of os.system
    """

    # the column(s) to keep in the output file
    if (keep_col and isinstance(keep_col, int)):
        keep_col_str = ",$%s" % keep_col
    elif (keep_col and isinstance(keep_col, list)):
        keep_col_str_list = []
        for i in keep_col:
            keep_col_str_list.append(",$%s" % i)
        keep_col_str = "".join(keep_col_str_list)
    else:
        keep_col_str = ""

    # define the output coordinates
    if width > 0:
        if summit_col:
            # output is summit +- the given width
            command_str = "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2+$%s-%s,$2+$%s+%s%s}'" % (
                summit_col, width, summit_col, width, keep_col_str)
        else:
            # output is the midpoint +- the given width
            command_str = "awk 'BEGIN{OFS=\"\\t\"}{mid=int(($2+$3)/2); print $1,mid-%s,mid+%s%s}'" % (
                width, width, keep_col_str)
    else:
        if summit_col:
            # output is only the summit
            command_str = "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2+$%s,$2+$%s+1%s}'" % (
                summit_col, summit_col, keep_col_str)
        else:
            # output is only the midpoint of start and end
            command_str = "awk 'BEGIN{OFS=\"\\t\"}{mid=int(($2+$3)/2); print $1,mid,mid+1%s}'" % (
                keep_col_str)

    # output coordinates and the required kept columns
    p = os.system("%s < %s | %s > %s" % (read_method, infile, command_str, outfile))
    return p


def rm_column_negative(infile, outfile, col, read_method="cat"):
    """
    Replace all negative value to 0 in the given column

    Parameters
    ----------------
    infile: str
    outfile: str
    col: int
        the column to keep non-negative values. (1-based)

    Returns
    ----------------
    int: return of os.system
    """

    cmd = "awk -v OFS=\"\t\" '{if($%d<0){$%d=0;print $0}else{print $0}}'" % (col, col)
    p = os.system("%s < %s | %s > %s" % (read_method, infile, cmd, outfile))
    return p