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
    met::Int,
    mej::Int,
)::Nothing

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

        align_dna(sa, fq1t, fq2t, fa, ba, n_jo, mej)
    
    elseif mo == "cdna"

        align_cdna(sa, fq1t, fq2t, fa, ba, n_jo, mej)

    end

    sp = splitext(fa)[1]

    fag = "$sp.bgz"

    if !isfile(fag)

        run(
            pipeline(
                `gzip --decompress $fa --stdout`,
                `bgzip --threads $n_jo --stdout`,
                fag,
            ),
        )

    end

    pav = joinpath(pao, "call_variant")

    call_variant(
        mo,
        ba,
        nothing,
        ta,
        fag,
        chs,
        chn,
        pav,
        n_jo,
        met,
        pas,
    )

    return nothing

end
