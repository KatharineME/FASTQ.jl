function trim(fq1, fq2, pa, n_jo)

    @assert make_directory(pa, "trim")

    ht = joinpath(pa, "fastp.html")

    js = joinpath(pa, "fastp.json")

    ou1 = joinpath(pa, "trimmed.R1.fastq.gz")

    ou2 = joinpath(pa, "trimmed.R2.fastq.gz")

    run(
        `fastp --detect_adapter_for_pe --thread $n_jo --json $js --html $ht --in1 $fq1 --in2 $fq2 --out1 $ou1 --out2 $ou2`,
    )

    println("\nTrimming finished\n")

    return

end
