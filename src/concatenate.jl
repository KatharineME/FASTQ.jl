function concatenate(fq_::Vector{Any}, id::String="R1")::Nothing

    fo_ = []

    re_ = []

    for fq in fq_

        if occursin(id, fq)

            push!(fo_, fq)

        elseif occursin(replace(id, "1" => "2"), fq)

            push!(re_, fq)

        end

    end

    n_fo = length(fo_)

    n_re = length(re_)

    println("Number of forward read files found = $n_fo")

    println("Number of reverse read files found = $n_re")

    println()

    sa = last(splitdir(dirname(fq_[1])))

    co = joinpath(dirname(dirname(fq_[1])), string(sa, "_concat"))

    if ispath(co)

        println(
            "Skipping concatenation because directory already exists: $co")

    elseif n_fo <= 1 && n_re <= 1

        println(
            "Nothing to concatenate. Number of forward reads and reverse reads are both <= 1.\n")

    else

        run(`mkdir $co`)

        println("Concatenating ...")

        println()

        gr_su = Dict(fo_ => "_R1.fastq.gz", re_ => "_R2.fastq.gz")

        for gr in keys(gr_su)

            run(
                pipeline(
                    `cat $gr`,
                    stdout = joinpath(co, string(sa, gr_su[gr])),
                ),
            )

        end

        println("Concatenated files saved at: $co")

    end

    return nothing

end
