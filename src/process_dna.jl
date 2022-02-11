function process_dna(mo, sa, fq1, fq2, ta, pao, fa, chs, chn, pas, n_jo, me, to)

    @assert make_directory(pao, "process germline dna")

    for pa in [fq1, fq2, fa, chs, chn, pas]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(pao, "trim/")

    trim(fq1, fq2, tr, n_jo)

    fq1t = joinpath(tr, TRIMMED_R1)

    fq2t = joinpath(tr, TRIMMED_R2)

    check_read([fq1t, fq2t], joinpath(pao, "check_trim"), n_jo)

    ba = joinpath(pao, "align", "germ.bam")

    if mo == "dna"

        align_dna(sa, fq1t, fq2t, fa, ba, n_jo, me)

    elseif mo == "cdna"

        align_cdna(sa, fq1t, fq2t, fa, ba, n_jo, me)

    end

    pav = joinpath(pao, "call_germline_variant")

    call_germline_variant(mo, ba, ta, fa, chs, chn, pav, n_jo, me, to, pas)

    return

end
