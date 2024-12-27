def inference_cmd(tool, out_dir, out_file, fasta_file, homer_rscript="MoVRs_Motif2meme.R"):
    if tool == "DREME":
        cmd = "dreme -oc %s -k 8 -m 5 -e 1000 -verbosity 1 -p %s; cp %s/dreme.txt %s\n" % (out_dir, fasta_file, out_dir, out_file)
    elif tool == "MEME":
        cmd = "meme -oc %s -dna -nmotifs 5 -w 8 -maxsize 250000 -nostatus %s; cp %s/meme.txt %s\n" % (out_dir, fasta_file, out_dir, out_file)
    elif tool == "HOMER":
        cmd = "findMotifs.pl %s fasta %s -nocheck -nogo -len 8 -noknown -basic -S 5; Rscript %s %s/homerMotifs.motifs8 %s\n" % (fasta_file, out_dir, homer_rscript, out_dir, out_file)
    elif tool == "STREME":
        cmd = "streme --p %s --oc %s --w 8 --nmotifs 5 --verbosity 1; cp %s/streme.txt %s\n" % (fasta_file, out_dir, out_dir, out_file)
    return(cmd)