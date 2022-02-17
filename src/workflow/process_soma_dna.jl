function process_soma_dna(se)

    n_jo, me, mo, ta, _, _, sa, to, ou, r1, r2, sor1, sor2, ge, _, chs, chn, sn = read_setting(se)

    pa = joinpath(ou, "process_somatic_dna")

    @assert make_directory(pa, "process somatic dna")

    for pa in [r1, r2, sor1, sor2, ge, chs, chn, sn]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    trge = joinpath(pa, "trim", "germline")

    gr1 = joinpath(trge, TRIMMED_R1)

    gr2 = joinpath(trge, TRIMMED_R2)

    trso = joinpath(pa, "trim", "somatic")

    sr1 = joinpath(trso, TRIMMED_R1)

    sr2 = joinpath(trso, TRIMMED_R2)

    for g in [[r1, r2, trge], [sor1, sor2, trso]]

        trim(g[1], g[2], g[3], n_jo)

    end

    check_read([gr1, gr2, sr1, sr2], joinpath(pa, "check_trim"), n_jo)

    al = joinpath(pa, "align_$mo")

    bage = joinpath(al, "$sa.germline.bam")

    baso = joinpath(al, "$sa.somatic.bam")

    for g in [[bage, gr2, gr2], [baso, sr1, sr2]]

        align_dna(al, sa, g[1], g[2], g[3], ge, n_jo, me)

    end

    pav = joinpath(pa, "call_somatic_variant")

    call_somatic_variant(ta, bage, baso, ge, chs, chn, pav, n_jo, me, to, sn)

    return

end
