function error_if_directory(pa)

    paf = Fastq.support.get_full_path(pa)

    na = replace(basename(paf), "_" => " ")

    if ispath(paf)

        error("\nSkipping $na because directory already exists:\n $pa\n")

    else

        println("Making $paf")

        mkpath(paf)

    end

end
