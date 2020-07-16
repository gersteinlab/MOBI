"""Function to parse ENCODE metadata table
"""

import os
import shutil
import pandas as pd

def _read_df_or_str(in_file):
    if isinstance(in_file, str):
        df = pd.read_csv(in_file, sep="\t", low_memory=False)
    elif isinstance(in_file, pd.DataFrame):
        df = in_file.copy()
    else:
        raise("Unknown input data type")
    return(df)


def add_TF(encode_file, outfile, species):
    """Add two column indicating the ChIP target TF
    "ENCODE TF" and "TF"

    Parameters
    ----------
    encode_file : str or pd.DataFrame
        the metafile download from encode
        required colunmn name in input: "Experiment target"
    outfile : str or None
        output file. If none, won't save
    species : str
        eg. "dmelanogaster", "celegans", "human". The string that was used as the prefix in the "Experiment target" column in ENCODE metadata

    Returns
    ----------------
    DataFrame
        The original DataFrame with two more columns:
        1. "ENCODE TF" the name of the TF and 2. "TF" the name of the TF w/o special char "(", ")" and ":"
    """

    # read input
    df = _read_df_or_str(encode_file)

    # add TF
    df["ENCODE TF"] = [i.split("-%s" % species)[0].replace("3xFLAG-","").replace("eGFP-","") for i in df["Experiment target"]]
    # add TF w/o special char
    df["TF"] = df["ENCODE TF"].str.replace("(", "").str.replace(")", "").str.replace(":", "")

    # save to file
    if outfile:
        df.to_csv(outfile, sep="\t", index=False)
    return df


def filter_subset(encode_file, outfile, rm_poly=True, rs=200, **kwargs):
    """Filter metadata to keep certain subset
    Output will be in random order

    Parameters
    ----------
    encode_file : str or pd.DataFrame
        required column: TF
    outfile: str or None
        If None, will not save to file
    rm_poly: boolean
        whether to remove polymerase in the TF list
        All TF names that starts with Pol, POL, pol will be removed
    rs: int
        random seed
    **kwargs:
        {col_name: value}, subset the metadata base on col_name=value,
        space in col_name should be replace with "_", case sensitive

    Returns
    -------
    DataFrame
        df after all the filters
    """

    # read input
    df = _read_df_or_str(encode_file)

    # subset based on kwargs
    if kwargs is not None:
        for key, value in kwargs.items():
            key = key.replace("_", " ")
            if key in df.columns:
                df = df[df[key] == value].copy()
            else:
                raise KeyError("kwargs key not in metadata columns")

    # released data only
    df = df[df["File Status"] == "released"]

    # remove polymerase
    if rm_poly:
        df = df[~df['TF'].str.startswith("POL")]
        df = df[~df['TF'].str.startswith("pol")]
        df = df[~df['TF'].str.startswith("Pol")]

    # output will be in random order
    df = df.sample(frac=1, random_state=rs).reset_index(drop=True).copy()

    if outfile:
        df.to_csv(outfile, sep="\t", index=False)

    return df


def add_motif(encode_file, motif_dir, outfile):
    """Add motif info to ENCODE metafile based on input "TF" column

    Parameters
    ----------------
    encode_file : str or DataFrame
        required column: "TF"
    motif_dir : str
        directory containning all motifs in .meme format or FIMO output of motif in .bed format
    outfile: str or None

    Returns
    -------
    DataFrame
        ENCODE metadata table with additional column "motifs"
    """

    # read input
    df = _read_df_or_str(encode_file)

    # motif file basename
    motifs = [i.replace(".bed", "") for i in os.listdir(motif_dir) if not i.startswith(".")]
    motifs = [i.replace(".meme", "") for i in motifs]
    # motif TF name w/o special characters
    motif_TF = [i.replace("(", "").replace(")", "").replace(":", "") for i in motifs]

    # motif info
    motif_df = pd.DataFrame([motif_TF, motifs], index=["TF", "motifs"]).transpose()
    # combine with encode meta file
    meta = pd.merge(df, motif_df, left_on="TF", right_on="TF", how="outer")
    meta = meta[~meta['File accession'].isnull()]

    if outfile:
        meta.to_csv(outfile, sep="\t", index=False)
    return meta


def download_file(encode_file, download_dir, ext=".bed.gz", existed_data_dir=None):
    """Download every files in the metadata df if not already in the download directory
    Output file name format: file_accession.extension
    For bed file, example name format: ENCFFxxxID.bed.gz
    To check integrity of the downloaded files, run: gzip -tv *.bed.gz

    Parameters
    ----------------
    encode_file : str or DataFrame
        The metadata df, required column: "File download URL", "File accession"
    download_dir : str
        the dir to save the files
    ext: str
        file extension in the naming
    existed_data_dir: str
        if given, will check if the needed file has already been in this dir
        if file existed, just move it (will move, not copy)

    Returns
    ----------------
    None
    """

    # read input
    df = _read_df_or_str(encode_file)

    dls = df['File download URL'].values
    file_accessions = df['File accession'].values

    for dl, file_accession in zip(dls, file_accessions):
        dst = os.path.join(download_dir, file_accession+ext) # final file name
        if not os.path.isfile(dst):
            # if file not already existed in the output dir
            if existed_data_dir:
                existed_dst = os.path.join(existed_data_dir, file_accession+ext)
                if os.path.isfile(existed_dst):
                    # existed_data_dir given and this TF exist, just move it to the new dir
                    shutil.move(existed_dst, dst)
                else:
                    # else download
                    os.system("curl -o %s -L %s" % (dst, dl))
            else:
                # if existed_data_dir not given, download directly
                os.system("curl -o %s -L %s" % (dst, dl))
    return


def unique_TF_parser(encode_file, motif=True, prior=None):
    """Filter to get unique and non-unique TF names and files from ENCODE metafile

    Parameters
    ----------------
    meta_file : str or DataFrame
    motif : boolean
        if true, report TF names with motif as well as not with motif.
        if false, only not with motif
    prior : str
        if given, when picking TFs that have multiple experiment for unique, will first pick from 'Biosample term name'==prior

    Returns
    ----------------
    tuple: (TF_files, TF_names, TF_files_w_motif, TF_names_w_motif, motif_files, u_TF_files, u_TF_names, u_TF_files_w_motif, u_TF_names_w_motif, u_motif_files)

    if motif=False, all "*_motif" will be None
    TF_files: ENCFFxxxx
    TF_names: gene/TF names
    motif_files: gene/TF name
    u_*: unique, no redundant TF (based on name)

    """

    # read input
    df = _read_df_or_str(encode_file)

    # if given prior set, move rows that match prior to top
    if prior:
        df_p = df[df['Biosample term name'] == prior]
        df_np = df[df['Biosample term name'] != prior].sample(frac=1, random_state=200)
        df = pd.concat([df_p, df_np]).reset_index(drop=True).copy()

    # keep unique only
    u_df = df.drop_duplicates('TF', keep='first')

    TF_files = df['File accession'].values
    TF_names = df['TF'].values

    u_TF_files = u_df['File accession'].values
    u_TF_names = u_df['TF'].values

    if motif:
        df_w_motif = df[~df["motifs"].isnull()].reset_index(drop=True).copy()
        u_df_w_motif = df_w_motif.drop_duplicates('TF', keep='first')

        TF_files_w_motif = df_w_motif['File accession'].values
        TF_names_w_motif = df_w_motif['TF'].values
        motif_files = df_w_motif['motifs'].values

        u_TF_files_w_motif = u_df_w_motif['File accession'].values
        u_TF_names_w_motif = u_df_w_motif['TF'].values
        u_motif_files = u_df_w_motif['motifs'].values

        return (TF_files, TF_names, TF_files_w_motif, TF_names_w_motif, motif_files, u_TF_files, u_TF_names, u_TF_files_w_motif, u_TF_names_w_motif, u_motif_files)
    else:
        return (TF_files, TF_names, None, None, None, u_TF_files, u_TF_names, None, None, None)
