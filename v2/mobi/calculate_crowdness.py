import os
import tempfile
import pandas as pd
import pybedtools

def bedfile_get_summit(infile, outfile, summit_col=10):
    # summit_col: the column of summit position from start
        
    df = pd.read_csv(infile, sep="\t", header=None)
    df[1] = df[1] + df[summit_col-1]
    df[2] = df[1] + 1
    df[3] = os.path.basename(infile)
    df.iloc[:,:4].to_csv(outfile, sep="\t", header=False, index=False)
    return


def bedfile_get_window(infile, outfile, width, genome_size, keep_col=[6], summit_col=10):
    # resize the bed regions to +- width bp centered on the summit
    
    df = pd.read_csv(infile, sep="\t", header=None)
    newdf = pd.DataFrame([])
    newdf[0] = df[0]
    newdf[1] = df[1]+df[summit_col-1]-width
    newdf[2] = df[1]+df[summit_col-1]+width
    newdf.loc[newdf[1]<0, 1] = 0 # negative start will change to 0
    newdf[3] = newdf[0].map(genome_size)
    newdf.loc[newdf[2]>newdf[3], 2] = newdf.loc[newdf[2]>newdf[3], 3] # for any end > genome size, end change to genome size
    del newdf[3]
    newdf = pd.concat([newdf, df.loc[:, keep_col]], axis=1)
    newdf[1] = newdf[1].astype(int)
    newdf[2] = newdf[2].astype(int)
    newdf.to_csv(outfile, sep="\t", header=False, index=False)
    
    
def get_crowdness(input_files, output_dir, crowdness_window, genome_size, subset=None, all_summit_file="./tmp_allSummit.txt", c_type="count"):
    # c_type: the available choice from bedtools groupby -c option
    
    # if the summit file does NOT exist, create one by using the binding sites regions
    if not os.path.exists(all_summit_file):
        with tempfile.TemporaryDirectory() as tmp_dir:
            for ii in input_files:
                bedfile_get_summit(ii, os.path.join(tmp_dir, os.path.basename(ii)))
            os.system("cat %s/* | sort -k1,1 -k2,2n > %s" % (tmp_dir, all_summit_file))

    # for each TF, 
    with tempfile.TemporaryDirectory() as tmp_dir:
        # if given subset, only calculate crowdness for the subset
        if subset:
            files = subset
        else:
            files = input_files
            
        for ii in files:
            bedfile_get_window(
                ii,
                os.path.join(tmp_dir, os.path.basename(ii)),
                crowdness_window,
                genome_size)

            bt = pybedtools.BedTool(os.path.join(tmp_dir, os.path.basename(ii)))
            all_summit_bt = pybedtools.BedTool(all_summit_file)
            bt_summit = bt.intersect(all_summit_bt, wao=True)
            bt_summit_groupby = bt_summit.groupby(g=[1,2,3,4], c=8, o=[c_type])
            bt_summit_groupby.saveas(os.path.join(output_dir, os.path.basename(ii)))
    
    if all_summit_file == "./tmp_allSummit.txt": # if using the default tmp allSummit file
        os.remove(all_summit_file) # delete tmp file