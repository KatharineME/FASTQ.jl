function process_dna(se)

    n_jo, me, mo, ta, sa, to, ou, ger1, ger2, _, _, ge, _, chs, chn, sn = read_setting(se)

    @assert make_directory(ou, "process dna")

    for pa in [ger1, ger2, ge, chs, chn, sn]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(ou, "trim/")

    trim(ger1, ger2, tr, n_jo)

    r1t = joinpath(tr, TRIMMED_R1)

    r2t = joinpath(tr, TRIMMED_R2)

    check_read([r1t, r2t], joinpath(ou, "check_trim"), n_jo)

    al = joinpath(ou, "align_$mo")

    ba = joinpath(al, "$sa.bam")

    if mo == "dna"

        align_dna(al, sa, ba, r1t, r2t, ge, n_jo, me)

    elseif mo == "cdna"

        align_cdna(al, sa, r1t, r2t, ge, n_jo)

    end

    pav = joinpath(ou, "call_germline_variant")

    call_germline_variant(mo, ta, ba, ge, chs, chn, pav, n_jo, me, to, sn)

    return

end
