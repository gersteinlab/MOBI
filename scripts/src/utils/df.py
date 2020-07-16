"""Some helper function for manipulating pandas.DataFrame
"""

import numpy as np
import pandas as pd

def col_to_dict(df, col_key, col_val):
    """Map tow columns in a DataFrame to dictionary

    Parameters
    ----------
    df : DataFrame
        input df
    col_key : str
        name of the column that will be used as dict key
    col_val : str
        name of the column that will be used as dict value

    Returns
    -------
    dict
        {col_key: col_val}
    """

    return(pd.Series(data = df[col_val].values, index=df[col_key]).to_dict())


def global_max_index(df):
    """Get the column and row index of the global max of a df

    Parameters
    ----------
    df : DataFrame
        input df

    Returns
    -------
    tuple
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


def split_dataframe_list(df, target_column, separator):
    ''' Split the str in target_column into separated rows.
    Found from https://gist.github.com/jlln/338b4b0b55bd6984f883

    Parameters
    ----------
    df : DataFrame
        input
    target_column : str or int
        index of the column containing the values to split
    separator : str
        the symbol used to perform the split

    Returns
    -------
    DataFrame
        with each entry for the target column separated, with each element moved into a new row.
        The values in the other columns are duplicated across the newly divided rows.
    '''
    def split_list_to_rows(row, row_accumulator, target_column, separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(split_list_to_rows, axis=1, args = (new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df
