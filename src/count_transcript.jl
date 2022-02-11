function count_transcript(fa, pa, n_jo, fq1, fq2 = nothing, fr = 51, sd = 0.05)

    @assert make_directory(pa, "count transcript")

    id = "$fa.kallisto_index"

    if !ispath(id)

        println("\nCreating kallisto index...")

        run(`kallisto index --index $id $fa`)

        println("\nMade kallisto index at $id\n")

    end

    if fq2 !== nothing

        println("Running paired end psuedoalignment")

        run(`kallisto quant --threads $n_jo --index $id --output-dir $pa $fq1 $fq2`)

    elseif fq2 === nothing

        println("Running single end psuedoalignment")

        run(
            `kallisto quant --single --fragment-length $fr --sd $sd --threads $n_jo --index $id --output-dir $pa $fq1`,
        )

    end

    return

end
