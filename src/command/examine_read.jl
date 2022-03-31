function examine_read(se)

    fe_va = read_setting(se)

    r1 = fe_va["r1"]

    r2 = fe_va["r2"]

    sor1 = fe_va["sor1"]

    sor2 = fe_va["sor2"]

    Fastq.fastq.find(dirname(r1))

    println()

    if isempty(sor1) || sor1 === nothing

        re_ = [r1, r2]

    else

        Fastq.fastq.find(dirname(sor1))

        re_ = [r1, r2, sor1, sor2]

    end

    Fastq.fastq.check_read(re_, joinpath(fe_va["ou"], "check_raw"), fe_va["n_jo"])

end
