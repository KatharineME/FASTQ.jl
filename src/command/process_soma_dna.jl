function process_soma_dna(se)

    fe_va = read_setting(se)

    pa = joinpath(fe_va["ou"], "process_somatic_dna")

    Fastq.support.error_if_directory(pa)

    r1, r2, sor1, sor2, n_jo, me, sa, ge, chs, chn, sn, mo, ta = fe_va["r1"],
    fe_va["r2"],
    fe_va["sor1"],
    fe_va["sor2"],
    fe_va["n_jo"],
    fe_va["me"],
    fe_va["sa"],
    fe_va["ge"],
    fe_va["chs"],
    fe_va["chn"],
    fe_va["sn"],
    fe_va["mo"],
    fe_va["ta"]

    for fi in [r1, r2, sor1, sor2, ge, chs, chn, sn]

        if !isfile(fi)

            error("$fi does not exist.")

        end

    end

    trge = joinpath(pa, "trim", "germline")

    gr1 = joinpath(trge, Fastq.TRIMMED_R1)

    gr2 = joinpath(trge, Fastq.TRIMMED_R2)

    trso = joinpath(pa, "trim", "somatic")

    sr1 = joinpath(trso, Fastq.TRIMMED_R1)

    sr2 = joinpath(trso, Fastq.TRIMMED_R2)

    for g in [[r1, r2, trge], [sor1, sor2, trso]]

        Fastq.fastq.trim(g[1], g[2], g[3], n_jo)

    end

    Fastq.fastq.check_read([gr1, gr2, sr1, sr2], joinpath(pa, "check_trim"), n_jo)

    alg = joinpath(pa, "align_$(mo)_germline")

    als = joinpath(pa, "align_$(mo)_somatic")

    bage = joinpath(alg, "$sa.bam")

    baso = joinpath(als, "$sa.bam")

    for g in [[alg, bage, gr2, gr2], [als, baso, sr1, sr2]]

        Fastq.fastq.align_dna(g[1], sa, g[2], g[3], g[4], ge, n_jo, me)

    end

    pav = joinpath(pa, "call_somatic_variant")

    bagem = joinpath(alg, "$sa.markdup.bam")

    basom = joinpath(als, "$sa.markdup.bam")

    Fastq.bam.call_somatic_variant(ta, bagem, basom, ge, chs, chn, pav, n_jo, me, fe_va["to"], sn)

end
