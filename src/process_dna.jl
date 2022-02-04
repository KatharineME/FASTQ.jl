function process_dna(
    mo::String,
    sa::String, #G2451
    fq1::String,
    fq2::String,
    ta::Bool,
    pao::String,
    fa::String,
    chs::String,
    chn::String,
    pas::String,
    n_jo::Int,
    me::Int,
    to::String,
)::Nothing

    if check_directory(pao, "process germline dna")

       return nothing

    end

    for pa in [fq1, fq2, fa, chs, chn, pas]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(pao, "trim/")

    trim(fq1, fq2, tr, n_jo)

    fq1t = joinpath(tr, "trimmed.R1.fastq.gz")

    fq2t = joinpath(tr, "trimmed.R2.fastq.gz")

    check_read([fq1t, fq2t], joinpath(pao, "check_trim"), n_jo)

    ba = joinpath(pao, "align", "germ.bam")

    if mo == "dna"

        align_dna(sa, fq1t, fq2t, fa, ba, n_jo, me)
    
    elseif mo == "cdna"

        align_cdna(sa, fq1t, fq2t, fa, ba, n_jo, me)

    end

    pav = joinpath(pao, "call_germline_variant")

    call_germline_variant(
        mo,
        ba,
        ta,
        fa,
        chs,
        chn,
        pav,
        n_jo,
        me,
        to,
        pas,
    )

    return nothing

end
