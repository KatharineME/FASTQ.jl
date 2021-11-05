function count_transcript(
    fa::String,
    pa::String,
    n_jo::Int64,
    fq1::String,
    fq2 = nothing,
    fr::Int64 = 51,
    sd::Float64 = 0.05,
)::Nothing

    println("Counting transcript\n")

    id = "$fa.kallisto_index"

    if !ispath(id)

        println("Creating kallisto index...\n")

        run(`kallisto index --index $id $fa`)

        println("Done running kallisto index.\n")

    end

    if ispath(id)

        println("kallisto index exists at: $id\n")

    end

    mkpath(pa)

    if fq2 !== nothing

        run(
            `kallisto quant --threads $n_jo --index $id --output-dir $pa $fq1 $fq2`,
        )

    elseif fq2 === nothing

        println("Running single end psuedoalignment")

        run(
            `kallisto quant --single --fragment-length $fr --sd $sd --threads $n_jo --index $id --output-dir $pa $fq1`,
        )

    end

    return nothing

end

export count_transcript
