function error_if_directory(pa)

    pa = Fastq.support.get_full_path(pa)

    pr = replace(basename(pa), "_", " ")

    if !ispath(pa)

        mkpath(pa)

        return true

    else

        error("\nSkipping $pr because directory already exists:\n $pa\n")

    end

end
