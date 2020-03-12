import os
import uuid
import numpy as np
import pandas as pd
from scipy.stats import binom_test

__all__ = ["df_colToDict",
           "df_global_max_index",
           "df_splitDataFrameList",
           "os_mkdir_tmp",
           "os_mkdir_empty",
           "os_remove_size_zero_file",
           "os_nline",
           "dsq_get_sh",
           "stat_sig_test",
           "stat_SE_bootstrap"]

def df_colToDict(df, col_key, col_val):
    """
    Map tow columns in a DataFrame to dictionary

    Parameters
    ----------------
    df: DataFrame
    col_key: str
        name of the column that will be used as dict key
    col_val: str
        name of the column that will be used as dict value
    
    Returns
    ----------------
    dict
    """

    return(pd.Series(data = df[col_val].values, index=df[col_key]).to_dict())


def df_global_max_index(df):
    """
    Get the column and row index of the global max of a df

    Parameters
    ----------------
    df: DataFrame

    Returns
    ----------------
    (col_index, row_index)
    """

    idxmaxdf1 = pd.DataFrame(df.idxmax()).reset_index()
    idxmaxdf1.columns = ["rank1", "rank_common"]

    idxmaxdf2 = pd.DataFrame(df.idxmax(axis=1)).reset_index()
    idxmaxdf2.columns = ["rank_common", "rank2"]

    idxmaxdf = pd.merge(idxmaxdf1, idxmaxdf2, left_on="rank_common", right_on="rank_common")
    idxmaxdf = idxmaxdf[idxmaxdf['rank1'] == idxmaxdf['rank2']].reset_index(drop=True)

    val = np.diag(np.array(df.loc[idxmaxdf.loc[:, 'rank_common'], idxmaxdf.loc[:, 'rank1']]))

    return((idxmaxdf.loc[val.argmax(), "rank1"], idxmaxdf.loc[val.argmax(), "rank_common"]))


def df_splitDataFrameList(df,target_column,separator):
    ''' 
    from https://gist.github.com/jlln/338b4b0b55bd6984f883
    
    Parameters
    ----------------
    df: DataFrame
    target_column: str or int
        the column containing the values to split
    separator: str
        the symbol used to perform the split
    
    Returns
    ----------------
    DataFrame
        with each entry for the target column separated, with each element moved into a new row. 
        The values in the other columns are duplicated across the newly divided rows.
    '''
    def splitListToRows(row,row_accumulator,target_column,separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows,axis=1,args = (new_rows,target_column,separator))
    new_df = pd.DataFrame(new_rows)
    return new_df


def os_mkdir_tmp(root_path="./"):
    """
    Make temporary dir with name "tmp_[randomStr]"

    Parameters
    ----------------
    root_path: str

    Returns
    ---------------
    str: full path to the tmp dir
    """
    
    dirname = "tmp_"+str(uuid.uuid4())
    os.makedirs(root_path+dirname, exist_ok=False)
    return(root_path+dirname+"/")


def os_mkdir_empty(path_to_dir, msg="Folder existed!"):
    """
    Make dir if folder not exist, or exist but with no file, or will all files size zero
    
    Parameters
    ----------------
    path_to_dir: str
        The folder to create
    msg: str
        The output message if folder already exist
    
    Returns
    ----------------
    1 if folder created. 0 if not created
    """

    if os.path.exists(path_to_dir):
        if not os.listdir(path_to_dir):
            pass  # if folder existed and empty then done
        else:
            files_size = [os.path.getsize(path_to_dir+"/"+i)
                          for i in os.listdir(path_to_dir)]
            if np.all(np.array(files_size) == 0):
                for i in os.listdir(path_to_dir):
                    # if folder existed and not empyt, but all file size 0, then remove all file
                    os.remove(path_to_dir+"/"+i)
            else:
                if msg:
                    print(msg)
                return(0)
    else:
        os.makedirs(path_to_dir)  # if folder not exist, then mkdir

    return(1)


def os_remove_size_zero_file(path):
    """
    Remove all file with size 0 in the path
    Only the current level, not recursively.

    Parameters
    ----------------
    path: str

    Returns
    ----------------
    None
    """

    for i in os.listdir(path):
        f = path+"/"+i
        if os.path.isfile(f) and os.path.getsize(f) == 0:
            os.remove(f)


def os_nline(filepath):
    count = 0
    with open(filepath) as f:
        for line in f.readlines():
            count += 1
    return count


def dsq_get_sh(jobfile, outfile=None, jobname=None, mem=10240, time="3-00:00:00", no_out_log=True):
    """
    Get dSQ sh file for joblist for dSQ v0.95
    If joblist have more than 9999 file, 
    it will automatically split and generate multiple joblist and sh files

    Parameters
    ----------------
    jobfile: str
        job list file path
    outfile: str
        dSQ sh file path
    jobname: str
    mem: int
    time: str
        in slurm supported format
    """
    if not jobname:
        jobname = os.path.basename(jobfile)
    if not outfile:
        outfile = os.path.dirname(jobfile)+"/dsq-"+os.path.basename(jobfile)+".sh"

    if no_out_log:
        out_log_str = "/dev/null"
    else:
        out_log_str = "dsq-{jobname}-%A_%4a-%N.out".format(jobname=jobname)
        
    nline = os_nline(jobfile) - 1
    if nline < 9999:
        with open(outfile, "w") as f:
            f.write("""#!/bin/bash
#SBATCH --partition=pi_gerstein
#SBATCH --mem-per-cpu={mem}
#SBATCH --time={time}
#SBATCH --mail-user=jiahao.gao@yale.edu
#SBATCH --mail-type=END
#SBATCH --output {out_log_str}
#SBATCH --array 0-{nline}%200
#SBATCH --job-name {jobname}

/ysm-gpfs/apps/software/dSQ/0.95/dSQBatch.py {in_dir}/{jobname} {in_dir}
""".format(mem=mem, time=time, jobname=jobname, in_dir=os.path.dirname(jobfile), nline=nline, out_log_str=out_log_str))

    elif nline >=9999:
        os.system("split -l 9800 {jobfile} -d -a 2 {in_dir}/{jobname}__".format(jobfile=jobfile, jobname=jobname, in_dir=os.path.dirname(jobfile)))
        jobfile_list = [i for i in os.listdir(os.path.dirname(jobfile)) if i.startswith(jobname+"__")]
        counts = [i.split("__")[-1] for i in jobfile_list]
        
        for c in counts:
            nline_c = os_nline(jobfile+"__"+c) - 1
            with open(outfile+c, "w") as f:
                f.write("""#!/bin/bash
#SBATCH --partition=pi_gerstein
#SBATCH --mem-per-cpu={mem}
#SBATCH --time={time}
#SBATCH --mail-user=jiahao.gao@yale.edu
#SBATCH --mail-type=END
#SBATCH --output {out_log_str}
#SBATCH --array 0-{nline}%200
#SBATCH --job-name {jobname}__{c}

/ysm-gpfs/apps/software/dSQ/0.95/dSQBatch.py {in_dir}/{jobname}__{c} {in_dir}
""".format(mem=mem, time=time, jobname=jobname, in_dir=os.path.dirname(jobfile), nline=nline_c, out_log_str=out_log_str, c=c))
    return


def stat_sig_test(array1, array2):
    x = np.sum((array1 - array2) > 0)
    n = np.sum((array1 - array2) != 0)
    pval = binom_test(x, n, alternative='two-sided')

    if (array1.mean() - array2.mean()) > 0:
        pval = pval/2
    else:
        pval = 1
    return pval


def stat_SE_bootstrap(stat):
    """
    Calculate the standard error of bootstrap estimates of the statistics

    Parameters
    ----------------
    stat: np.array
        list of the statistics from bootstrap, len = bootstrap times

    Returns
    ----------------
    float
    """
    return np.sqrt(np.sum((stat - stat.mean())**2) / (len(stat)-1))