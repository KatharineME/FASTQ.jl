function examine_read(se)

    n_jo, _, _, _, _, _, _, _, ou, r1, r2, sor1, sor2, _, _, _, _, _ = read_setting(se)

    Fastq.fastq.find(dirname(r1))

    println()

    if isempty(sor1) || sor1 === nothing

        re_ = [r1, r2]

    else

        Fastq.fastq.find(dirname(sor1))

        re_ = [r1, r2, sor1, sor2]

    end

    Fastq.fastq.check_read(re_, joinpath(ou, "check_raw"), n_jo)

    return

end
