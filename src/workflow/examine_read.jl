function examine_read(se)

    n_jo, _, _, _, _, _, ou, ger1, ger2, sor1, sor2, _, _, _, _, _ = read_setting(se)
    
    find(dirname(ger1))

    if ger1 !== nothing && sor1 !== nothing

        re_ = [ger1, ger2, sor1, sor2]

        find(dirname(sor1))
        

    else

        re_ = [ger1, ger2]

    end

    check_read(re_, joinpath(ou, "check_raw"), n_jo)

    return

end
