function apply_cdna_to_genome(se)

    fe_va = read_setting(se)

    println("this is fe_va: $fe_va")

    ou = fe_va["output_directory"]

    println("this is ou: $ou")

    pou = joinpath(ou, "apply_cdna_to_genome")

    println("this is pou: $pou")

    Fastq.support.error_if_directory(pou)

    pac = joinpath(pou, "align_cdna/")

    Fastq.support.error_if_directory(pac)

    re_ = Fastq.fastq.find(fe_va["cdna_read_directory"])

    Fastq.fastq.check_read(re_, joinpath(pou, "check_read"), fe_va["number_of_jobs"])

    Fastq.fastq.align_cdna_samples(
        pac,
        fe_va["cdna_read_directory"],
        fe_va["reference_genome"],
        fe_va["number_of_jobs"],
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
                    fe_va["annotate_with_rsid"],
                    fe_va["variant_database"],
                )

            end

        end

    end

end
