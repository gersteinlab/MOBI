"""dSQ array fuction
"""

import os
from . import os_func

def get_sh(jobfile, outfile=None, jobname=None, mem=10240, time="3-00:00:00", no_out_log=True):
    """Get dSQ sh file for joblist for dSQ v0.95
    If joblist have more than 9999 file,
    it will automatically split and generate multiple joblist and sh files

    Parameters
    ----------
    jobfile : str
        job list file path
    outfile : str
        dSQ sh file path
    jobname : str
    mem : int
    time : str
        in slurm supported format

    Returns
    -------
    None
    """
    if not jobname:
        jobname = os.path.basename(jobfile)
    if not outfile:
        outfile = os.path.join(os.path.dirname(jobfile), "dsq-"+os.path.basename(jobfile)+".sh")

    if no_out_log:
        out_log_str = "/dev/null"
    else:
        out_log_str = "dsq-{jobname}-%A_%4a-%N.out".format(jobname=jobname)

    nline = os_func.os_nline(jobfile) - 1
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
            nline_c = os_func.os_nline(jobfile+"__"+c) - 1
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
