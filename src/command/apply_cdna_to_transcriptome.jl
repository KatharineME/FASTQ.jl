function apply_cdna_to_transcriptome(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_transcriptome,
    fragment_length,
    fragment_length_standard_deviation,
    organism,
    mouse_transcript_to_mouse_gene,
)

    pou = joinpath(output_directory, "apply_cdna_to_transcriptome")

    Fastq.support.error_if_directory(pou)

    pap = joinpath(pou, "psuedoalign/")

    Fastq.support.error_if_directory(pap)

    re_ = Fastq.fastq.find(cdna_read_directory)

    Fastq.fastq.check_read(re_, joinpath(pou, "check_read"), number_of_jobs)

    Fastq.fastq.align_cdna_samples(
        pap,
        cdna_read_directory,
        reference_transcriptome,
        number_of_jobs,
        al = "transcriptome",
        fr = fragment_length,
        sd = fragment_length_standard_deviation,
    )

    Fastq.abundance.make_gene_by_sample(pap, pou, organism, mouse_transcript_to_mouse_gene)

end
