function apply_cdna_to_genome(se)

    fe_va = read_setting(se)

    re_ = Fastq.fastq.find()

    Fastq.fastq.align_cdna()

    Fastq.bam.call_germline_variant()

end
