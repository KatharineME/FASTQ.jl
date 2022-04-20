function apply_cdna_to_transcriptome(se)

    fe_va = read_setting(se)

    pou = joinpath(fe_va["output_directory"], "apply_cdna_to_transcriptome")

    Fastq.support.error_if_directory(pou)

    pap = joinpath(pou, "psuedoalign/")

    Fastq.support.error_if_directory(pap)

    re_ = Fastq.fastq.find(fe_va["cdna_read_directory"])

    Fastq.fastq.check_read(re_, joinpath(pou, "check_read"), fe_va["number_of_jobs"])

    Fastq.fastq.align_cdna_samples(
        pap,
        fe_va["cdna_read_directory"],
        fe_va["reference_transcriptome"],
        fe_va["number_of_jobs"],
        al = "transcriptome",
        fr = fe_va["fragment_length"],
        sd = fe_va["fragment_length_standard_deviation"],
    )

    Fastq.abundance.make_gene_by_sample(
        pap,
        pou,
        fe_va["organism"],
        fe_va["mouse_transcript_to_mouse_gene"],
    )

end
