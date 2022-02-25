function trim(r1, r2, pa, n_jo)

    Fastq.support.error_if_directory(pa)

    ht = joinpath(pa, "fastp.html")

    js = joinpath(pa, "fastp.json")

    ou1 = joinpath(pa, Fastq.TRIMMED_R1)

    ou2 = joinpath(pa, Fastq.TRIMMED_R2)

    run(
        `fastp --detect_adapter_for_pe --thread $n_jo --json $js --html $ht --in1 $r1 --in2 $r2 --out1 $ou1 --out2 $ou2`,
    )

    println("\nTrimming finished\n")

    return

end
