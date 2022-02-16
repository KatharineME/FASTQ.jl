function process_soma_dna(se)

    n_jo, me, mo, ta, sa, to, ou, ger1, ger2, sor1, sor2, ge, _, chs, chn, sn = read_setting(se)

    @assert make_directory(ou, "process somatic dna")

    for pa in [ger1, ger2, sor1, sor2, ge, chs, chn, sn]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    trge = joinpath(ou, "trim", "germline")

    gr1 = joinpath(trge, TRIMMED_R1)

    gr2 = joinpath(trge, TRIMMED_R2)

    trso = joinpath(ou, "trim", "somatic")

    sr1 = joinpath(trso, TRIMMED_R1) 

    sr2 = joinpath(trso, TRIMMED_R2)
   
    for g in [[ger1, ger2, trge], [sor1, sor2, trso]]

        trim(g[1], g[2], g[3], n_jo)

    end
    
    check_read([gr1, gr2, sr1, sr2], joinpath(ou, "check_trim"), n_jo)

    al = joinpath(ou, "align_$mo")

    bage = joinpath(al, "$sa.germline.bam")

    baso = joinpath(al, "$sa.somatic.bam")

    for g in [[bage, gr2, gr2], [baso, sr1, sr2]]

        align_dna(al, sa, g[1], g[2], g[3], ge, n_jo, me)

    end

    pav = joinpath(ou, "call_somatic_variant")

    call_somatic_variant(ta, bage, baso, ge, chs, chn, pav, n_jo, me, to, sn)

    return

end
