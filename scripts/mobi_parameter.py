class Parameter(object):
    def __init__(self):

    	## User input parameters and paths

        # name for this run
        self.job_id = "GM12878"
        # number of top peaks that will be used for all analysis
        self.data_n_peaks = 3000
        # window size for computing crowdness score
        self.crowdness_n_sliding = 250
        # number of sequences in the input for motif inference
        self.inference_n_seq = 250

        # genome fasta file
        self.genome_file = "/home/jg2447/slayman/data/genome/GRCh37.p13.genome.fa"
        # dir containing raw ChIP-seq bed from ENCODE
        self.chip_data_dir = "/home/jg2447/slayman/data/ENCODE_ChIP_seq_bed/human_hg19/"
        # clean ENCODE metadata file
        self.metafile = "/home/jg2447/slayman/motif_inference/result/metadata/human-GM12878.txt"
        # FIMO result of all known motifs
        self.fimo_dir = "/home/jg2447/slayman/data/FIMO/human_stack/FIMOp_0.000100_chr/"
        # meme format of all known motifs
        self.known_motif_dir = "/home/jg2447/slayman/data/motif_cisbp/meme_human_stack/"

        # main result directory
        self.result_folder = "/home/jg2447/slayman/motif_inference/result/MOBI/humanGM12878/"

        ##
        # all ChIP-seq bed files that will be used
        self.TFdata_dir = "%s/TFdata/" % self.result_folder
        # file containing all summits
        self.all_summit_file = "%s/all_summit.txt" % self.result_folder
        # containing spp, c-score for each TF
        self.crowdness_dir = "%s/crowdness/" % self.result_folder
        # containing spp, c-score, known-motif-occurance, TF-name for each TF for each desired binding site length
        self.motif_count_dir = "%s/motif_count/" % self.result_folder
        # rank files that aggregate all TF for each ranking, each length
        self.rank_file_dir = "%s/rank_file/" % self.result_folder
        # known motif enrichment directory
        self.enrichment_dir = "%s/enrichment/" % self.result_folder
        # inference input sequences in bed format
        self.bed_dir = "%s/bed_file/" % self.result_folder
        # inference input sequences in fasta format
        self.fasta_dir = "%s/fasta_file/" % self.result_folder
        # inferred motifs
        self.motif_dir = "%s/inference_motif/" % self.result_folder
        # inferred motifs vs known motifs, tomtom result
        self.tomtom_dir = "%s/inference_tomtom/" % self.result_folder
        # tomtom summary
        self.tomtom_summary_dir = "%s/tomtom_summary/" % self.result_folder
        # performance statistics
        self.performance_dir = "%s/performance/" % self.result_folder
