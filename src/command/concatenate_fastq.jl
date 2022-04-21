function concatenate_fastq(se)

    fe_va = Fastq.command.read_setting(se)

    pa = joinpath(fe_va["output_directory"], "concatenate_fastq")

    Fastq.support.error_if_directory(pa)

    re_ = Fastq.fastq.find(fe_va["dna_read_directory"])

    Fastq.fastq.concatenate(re_, fe_va["read_name_scheme"])

end
