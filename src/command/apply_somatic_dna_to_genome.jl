function apply_somatic_dna_to_genome(
    output_directory,
    read1,
    read2,
    somatic_read1,
    somatic_read2,
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

    pa = joinpath(output_directory, "apply_somatic_dna_to_genome")

    Fastq.support.error_if_directory(pa)

    Fastq.fastq.examine_read(read1, read2, pa, number_of_jobs, somatic_read1, somatic_read2)

    for fi in [
        read1,
        read2,
        somatic_read1,
        somatic_read2,
        reference_genome,
        chromosome_position,
        chromosome_name,
        snpeff,
    ]

        if !isfile(fi)

            error("$fi does not exist.")

        end

    end

    trge = joinpath(pa, "trim", "germline")

    gr1 = joinpath(trge, Fastq.TRIMMED_R1)

    gr2 = joinpath(trge, Fastq.TRIMMED_R2)

    trso = joinpath(pa, "trim", "somatic")

    sr1 = joinpath(trso, Fastq.TRIMMED_R1)

    sr2 = joinpath(trso, Fastq.TRIMMED_R2)

    for g in [[read1, read2, trge], [somatic_read1, somatic_read2, trso]]

        Fastq.fastq.trim(g[1], g[2], g[3], number_of_jobs)

    end

    Fastq.fastq.check_read([gr1, gr2, sr1, sr2], joinpath(pa, "check_trim"), number_of_jobs)

    alg = joinpath(pa, "align_$(molecule)_germline")

    als = joinpath(pa, "align_$(molecule)_somatic")

    bage = joinpath(alg, "$sample.bam")

    baso = joinpath(als, "$sample.bam")

    for g in [[alg, bage, gr1, gr2], [als, baso, sr1, sr2]]

        Fastq.fastq.align_dna(
            g[1],
            sample,
            g[2],
            g[3],
            g[4],
            reference_genome,
            number_of_jobs,
            memory,
        )

    end

    pav = joinpath(pa, "call_somatic_variant")

    bagem = joinpath(alg, "$sample.bam")

    basom = joinpath(als, "$sample.bam")

    Fastq.bam.call_somatic_variant(
        exome,
        bagem,
        basom,
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
