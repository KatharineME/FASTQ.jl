function apply_germline_dna_to_genome(se)

    fe_va = read_setting(se)

    pa = joinpath(fe_va["output_directory"], "apply_germline_dna_to_genome")

    Fastq.support.error_if_directory(pa)

    r1, r2, n_jo, me, sa, ge, chs, chn, sn = fe_va["read1"],
    fe_va["read2"],
    fe_va["number_of_jobs"],
    fe_va["memory"],
    fe_va["sample"],
    fe_va["reference_genome"],
    fe_va["chromosome_position"],
    fe_va["chromosome_name"],
    fe_va["snpeff"]

    for pa in [r1, r2, ge, chs, chn, sn]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(pa, "trim/")

    Fastq.fastq.trim(r1, r2, tr, n_jo)

    r1t = joinpath(tr, Fastq.TRIMMED_R1)

    r2t = joinpath(tr, Fastq.TRIMMED_R2)

    Fastq.fastq.check_read([r1t, r2t], joinpath(pa, "check_trim"), n_jo)

    al = joinpath(pa, "align_dna")

    ba = joinpath(al, "$sa.bam")

    Fastq.fastq.align_dna(al, sa, ba, r1t, r2t, ge, n_jo, me)

    pav = joinpath(pa, "call_germline_variant")

    Fastq.bam.call_germline_variant(
        fe_va["molecule"],
        fe_va["exome"],
        ba,
        ge,
        chs,
        chn,
        pav,
        n_jo,
        me,
        fe_va["tool_directory"],
        sn,
    )

end
