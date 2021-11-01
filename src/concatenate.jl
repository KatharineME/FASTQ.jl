function concatenate(fq_::Vector{Any}, id::String="R1")::Nothing

    if id == "R1"

        println("Using default naming scheme: \"R1\" for forward and \"R2\" for reverse")

        println()

    end

    fo_ = []

    re_ = []

    for fi in fq_

        if occursin(id, fi)

            push!(fo_, fi)

        end

        if occursin(replace(id, "1" => "2"), fi)

            push!(re_, fi)

        end

    end

    n_fo = length(fo_)

    n_re = length(re_)

    println("Number of forward read files found = $n_fo")

    println("Number of reverse read files found = $n_re")

    println()

    sa = last(splitdir(dirname(fq_[1])))

    dica = joinpath(dirname(dirname(fq_[1])), string(sa, "_concat"))
    
    if ispath(dica)

        println(
            "Skipping concatenation because directory already exists:\n $dica\n",
        )

    elseif n_fo <= 1 && n_re <= 1

        println(
            "Nothing to concatenate. Number of forward reads and reverse reads are both <= 1.\n",
        )

    else

        run(`mkdir $dica`)

        println("Concatenating read files...\n")

        gr_su = Dict(fo_ => "_R1.fastq.gz", re_ => "_R2.fastq.gz")

        for gr in keys(gr_su)

            run(
                pipeline(
                    `cat $gr`,
                    stdout = joinpath(dica, string(sa, gr_su[gr])),
                ),
            )

        end

        println("Concatenated files saved at: $dica")

    end

    return nothing

end

export concatenate
