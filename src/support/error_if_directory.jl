function error_if_directory(pa)

    paf = Fastq.support.get_full_path(pa)

    pr = replace(basename(paf), "_", " ")

    if !ispath(paf)

        mkpath(paf)

    else

        error("\nSkipping $pr because directory already exists:\n $pa\n")

    end

end
