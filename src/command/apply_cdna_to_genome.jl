function apply_cdna_to_genome(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_genome,
    molecule,
    exome,
    chromosome_position,
    chromosome_name,
    memory,
    tool_directory,
    snpeff,
    annotate_with_rsid,
    variant_database,
)

    pou = joinpath(output_directory, "apply_cdna_to_genome")

    Fastq.support.error_if_directory(pou)

    pac = joinpath(pou, "align_cdna/")

    Fastq.support.error_if_directory(pac)

    re_ = Fastq.fastq.find(cdna_read_directory)

    Fastq.fastq.check_read(re_, joinpath(pou, "check_read"), number_of_jobs)

    Fastq.fastq.align_cdna_samples(
        pac,
        cdna_read_directory,
        reference_genome,
        number_of_jobs,
        al = "genome",
    )

    pav = joinpath(pou, "call_germline_variant")

    for (ro, di_, fi_) in walkdir(pac)

        for fi in fi_
            if endswith(fi, ".bam")

                ba = joinpath(ro, fi)

                sa = basename(splitdir(ba)[1])

                pas = joinpath(pav, sa)

                Fastq.bam.call_germline_variant(
                    molecule,
                    exome,
                    ba,
                    reference_genome,
                    chromosome_position,
                    chromosome_name,
                    pas,
                    number_of_jobs,
                    memory,
                    tool_directory,
                    snpeff,
                    annotate_with_rsid,
                    variant_database,
                )

            end

        end

    end

end
