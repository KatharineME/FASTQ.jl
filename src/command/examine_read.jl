function examine_read(se)

    fe_va = read_setting(se)

    r1 = fe_va["read1"]

    r2 = fe_va["read2"]

    sor1 = fe_va["somatic_read1"]

    sor2 = fe_va["somatic_read2"]

    Fastq.fastq.find(dirname(r1))

    if isempty(sor1) || sor1 === nothing

        re_ = [r1, r2]

    else

        Fastq.fastq.find(dirname(sor1))

        re_ = [r1, r2, sor1, sor2]

    end

    Fastq.fastq.check_read(
        re_,
        joinpath(fe_va["output_directory"], "examine_read"),
        fe_va["number_of_jobs"],
    )

end
