function process_dna(se)

    n_jo, me, mo, ta, _, _, sa, to, ou, r1, r2, _, _, ge, _, chs, chn, sn = read_setting(se)

    pa = joinpath(ou, "process_dna")

    Fastq.support.error_if_directory(pa)

    for pa in [r1, r2, ge, chs, chn, sn]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(pa, "trim/")

    Fastq.fastq.trim(r1, r2, tr, n_jo)

    r1t = joinpath(tr, Fastq.TRIMMED_R1)

    r2t = joinpath(tr, Fastq.TRIMMED_R2)

    Fastq.fastq.check_read([r1t, r2t], joinpath(pa, "check_trim"), n_jo)

    al = joinpath(pa, "align_$mo")

    ba = joinpath(al, "$sa.bam")

    if mo == "dna"

        Fastq.fastq.align_dna(al, sa, ba, r1t, r2t, ge, n_jo, me)

    elseif mo == "cdna"

        Fastq.fastq.align_cdna(al, sa, r1t, r2t, ge, n_jo)

    end

    bam = joinpath(al, "$sa.markdup.bam")

    pav = joinpath(pa, "call_germline_variant")

    Fastq.bam.call_germline_variant(mo, ta, bam, ge, chs, chn, pav, n_jo, me, to, sn)

    return

end
