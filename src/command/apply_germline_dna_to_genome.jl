function apply_germline_dna_to_genome(
    output_directory,
    read1,
    read2,
    number_of_jobs,
    memory,
    sample,
    reference_genome,
    chromosome_position,
    chromosome_name,
    snpeff,
    molecule,
    exome,
    tool_directory,
    annotate_with_rsid,
    variant_database,
)

    pa = joinpath(output_directory, "apply_germline_dna_to_genome")

    Fastq.support.error_if_directory(pa)

    Fastq.fastq.examine_read(read1, read2, pa, number_of_jobs)

    for pa in [read1, read2, reference_genome, chromosome_position, chromosome_name, snpeff]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(pa, "trim/")

    Fastq.fastq.trim(read1, read2, tr, number_of_jobs)

    r1t = joinpath(tr, Fastq.TRIMMED_R1)

    r2t = joinpath(tr, Fastq.TRIMMED_R2)

    Fastq.fastq.check_read([r1t, r2t], joinpath(pa, "check_trim"), number_of_jobs)

    al = joinpath(pa, "align_dna")

    ba = joinpath(al, "$sample.bam")

    Fastq.fastq.align_dna(al, sample, ba, r1t, r2t, reference_genome, number_of_jobs, memory)

    pav = joinpath(pa, "call_germline_variant")

    Fastq.bam.call_germline_variant(
        molecule,
        exome,
        ba,
        reference_genome,
        chromosome_position,
        chromosome_name,
        pav,
        number_of_jobs,
        memory,
        tool_directory,
        snpeff,
        annotate_with_rsid,
        variant_database,
    )

end
