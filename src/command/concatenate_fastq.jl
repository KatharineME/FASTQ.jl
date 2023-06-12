function concatenate_fastq(dna_read_directory, read_name_scheme)

    re_ = Fastq.fastq.find(dna_read_directory)

    Fastq.fastq.concatenate(re_, read_name_scheme)

end
