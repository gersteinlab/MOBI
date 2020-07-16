"""Helper functions for os and files
"""

import os
import uuid
import numpy as np

def remove_size_zero_file(path):
    """Remove all file with size 0 in the path
    Only the current level, not recursively.

    Parameters
    ----------
    path : str
        path

    Returns
    -------
    None
    """

    for i in os.listdir(path):
        f = os.path.join(path, i)
        if os.path.isfile(f) and os.path.getsize(f) == 0:
            os.remove(f)
    return


def mkdir_tmp(root_path="./"):
    """Make temporary dir with name "tmp_[randomStr]"

    Parameters
    ----------
    root_path : str
        root path

    Returns
    -------
    str
        full path to the tmp dir
    """

    dirname = os.path.join("tmp_", str(uuid.uuid4()))
    os.makedirs(os.path.join(root_path, dirname), exist_ok=False)
    return(os.path.join(root_path, dirname))


def mkdir_empty(path_to_dir, msg="Folder existed!"):
    """Makedir if folder not exist, or exist but with no file, or all files have size zero

    Parameters
    ----------
    path_to_dir : str
        path_to_dir
    msg : str
        error message when folder existed

    Returns
    -------
    None
    """

    if os.path.exists(path_to_dir):
        if not os.listdir(path_to_dir):
            pass  # if folder existed and empty then done
        else:
        	existed_files = [os.path.join(path_to_dir, i) for i in os.listdir(path_to_dir)]
            files_size = [os.path.getsize(i) for i in existed_files]
            if np.all(np.array(files_size) == 0):
                # if folder existed and not empyt, but all file size 0, then remove all file
                __ = [os.remove(i) for i in existed_files]
            else:
                # if folder existed and have non-zero size files, then raise error
                raise OSError(msg)
    else:
        os.makedirs(path_to_dir)  # if folder not exist, then mkdir
    return


def os_nline(filepath):
    """Count file lines

    Parameters
    ----------
    filepath : str
        file

    Returns
    -------
    int
        line count
    """
    count = 0
    with open(filepath) as f:
        for line in f.readlines():
            count += 1
    return count
