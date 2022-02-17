function examine_read(se)

    n_jo, _, _, _, _, _, _, _, ou, r1, r2, sor1, sor2, _, _, _, _, _ = read_setting(se)

    find(dirname(r1))

    println()

    if sor1 !== nothing
        
        find(dirname(sor1))

        re_ = [r1, r2, sor1, sor2]

    else

        re_ = [r1, r2]

    end

    check_read(re_, joinpath(ou, "check_raw"), n_jo)

    return

end
