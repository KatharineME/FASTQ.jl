function trim(
    fq1::String,
    fq2::String,
    pa::String,
    n_jo::Int,
)::Nothing

    if ispath(pa)

        println("Skipping trimming because directory already exists: $pa")

        return nothing

    else
        
        mkpath(pa)

        ht = joinpath(pa, "fastp.html")

        js = joinpath(pa, "fastp.json")
        
        ou1 = joinpath(pa, basename(fq1))

        ou2 = joinpath(pa, basename(fq2))

        run(`fastp --detect_adapter_for_pe --json $js --html $ht --in1 $fq1 --in2 $fq2 --out1 $ou1 --out2 $ou2`)

        println("\nTrimming finished\n")

    end

    return nothing

end
