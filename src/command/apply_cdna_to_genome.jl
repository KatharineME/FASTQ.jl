function apply_cdna_to_genome(se)

    fe_va = read_setting(se)

    ou = fe_va["output_directory"]

    pou = joinpath(ou, "apply_cdna_to_genome")

    Fastq.support.error_if_directory(pou)

    pac = joinpath(pou, "align_cdna/")

    Fastq.support.error_if_directory(pac)

    Fastq.fastq.align_cdna_samples(
        pac,
        fe_va["cdna_read_directory"],
        fe_va["reference_genome"],
        fe_va["number_of_jobs"],
        al = "genome",
    )

    # for each sample bam file call

    pav = joinpath(pou, "call_germline_variant")

    for (ro, di_, fi_) in walkdir(pac)

        for fi in fi_
            if endswith(fi, ".bam")

                ba = joinpath(ro, fi)

                sa = basename(splitdir(ba)[1])

                pas = joinpath(pav, sa)

                println(ba)

                println(sa)

                println(pas)

                Fastq.bam.call_germline_variant(
                    fe_va["molecule"],
                    fe_va["exome"],
                    ba,
                    fe_va["reference_genome"],
                    fe_va["chromosome_position"],
                    fe_va["chromosome_name"],
                    pas,
                    fe_va["number_of_jobs"],
                    fe_va["memory"],
                    fe_va["tool_directory"],
                    fe_va["snpeff"],
                )

            end

        end

    end

end
