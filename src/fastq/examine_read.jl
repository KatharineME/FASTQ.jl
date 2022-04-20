function examine_read(r1, r2, pa, n_jo, sor1 = nothing, sor2 = nothing)

    if sor1 === nothing

        re_ = [r1, r2]

    else

        re_ = [r1, r2, sor1, sor2]

    end

    Fastq.fastq.check_read(re_, joinpath(pa, "check_raw"), n_jo)

end
