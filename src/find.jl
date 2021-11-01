function find(di::String)::Vector{String}::Array{Any, 1}

    re_ = []

    na_n = ([(".fq", 0), (".fastq", 0), ("fq.gz", 0), ("fastq.gz", 0)])

    for (ro, di_, fi_) in walkdir(di)

       for fi in fi_

            if !occursin(".md5", fi)

                for (na, n) in na_n

                    if endswith(fi, na)

                        na_n[na] += 1

                        if endswith(fi, ".gz")

                            push!(re_, joinpath(ro, fi))

                        end

                    end

                end

            end

        end

    end

    println(File report: $na_n)

    return re_

end

export find
