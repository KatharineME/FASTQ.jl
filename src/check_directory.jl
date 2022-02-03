function check_directory(pa::String, pr::String)::Bool

    pa = abspath(pa)

    if ispath(pa)

        println("\nSkipping $pr because directory already exists:\n $pa\n") 

        bo = true

    else

        mkdir(pa)

        bo = false

    end

    return bo

end
