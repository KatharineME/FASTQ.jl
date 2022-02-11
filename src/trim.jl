function trim(fq1::String, fq2::String, pa::String, n_jo::Int)::Nothing

    if check_directory(pa, "trim")

        return nothing

    end

    ht = joinpath(pa, "fastp.html")

    js = joinpath(pa, "fastp.json")

    ou1 = joinpath(pa, "trimmed.R1.fastq.gz")

    ou2 = joinpath(pa, "trimmed.R2.fastq.gz")

    run(
        `fastp --detect_adapter_for_pe --json $js --html $ht --in1 $fq1 --in2 $fq2 --out1 $ou1 --out2 $ou2`,
    )

    println("\nTrimming finished\n")

    return nothing

end
