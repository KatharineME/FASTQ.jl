function find(di::String)

    re_ = []

    na_n = Dict(".fq" => 0, ".fastq" => 0, "fq.gz" => 0, "fastq.gz" => 0)

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

    println("Files found in $di:\n")
    
    for na in na_n println(na)

    end

    return re_

end
