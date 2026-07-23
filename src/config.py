class Config:
    def __init__(self, bam_readcount_arg, reference_arg, csv_arg, bams_directory_arg, mapq_arg, extension_arg, lower_threshold_arg, upper_threshold_arg, output_dir_arg):
        self.bam_readcount_arg = bam_readcount_arg
        self.reference_arg = reference_arg
        self.csv_arg = csv_arg
        self.bams_directory_arg = bams_directory_arg
        self.mapq_arg = mapq_arg
        self.extension_arg = extension_arg
        self.lower_threshold_arg = lower_threshold_arg
        self.upper_threshold_arg = upper_threshold_arg
        self.output_dir_arg = output_dir_arg
