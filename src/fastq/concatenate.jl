function concatenate(fq_, na = "R1")

    fo_ = []

    re_ = []

    for fq in fq_

        if occursin(na, fq)

            push!(fo_, fq)

        elseif occursin(replace(na, "1" => "2"), fq)

            push!(re_, fq)

        end

    end

    n_fo = length(fo_)

    n_re = length(re_)

    println("Number of forward read files found = $n_fo\n")

    println("Number of reverse read files found = $n_re\n")

    sa = last(splitdir(dirname(fq_[1])))

    co = joinpath(dirname(dirname(fq_[1])), string(sa, "_concat"))

    if n_fo <= 1 && n_re <= 1

        println(
            "\nNothing to concatenate. Number of forward reads and reverse reads are both <= 1.\n",
        )

    else

        Fastq.support.error_if_directory(co)

        println("\nConcatenating ...\n")

        gr_su = Dict(fo_ => "_R1.fastq.gz", re_ => "_R2.fastq.gz")

        for gr in keys(gr_su)

            run(pipeline(`cat $gr`, stdout = joinpath(co, string(sa, gr_su[gr]))))

        end

        println("\nConcatenated files saved at: $co\n")

    end

end
