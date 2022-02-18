function make_directory(pa, pr)

    pa = get_full_path(pa)

    if !ispath(pa)

        mkpath(pa)

        return true

    else

        error("\nSkipping $pr because directory already exists:\n $pa\n")

    end

end
