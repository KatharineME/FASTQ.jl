function count_transcript(se)

    n_jo, _, _, _, fr, sd, _, _, ou, ger1, ger2, _, _, _, tr, _, _, _ = read_setting(se)

    @assert make_directory(ou, "count transcript")

    id = "$tr.kallisto_index"

    if !ispath(id)

        println("\nCreating kallisto index...")

        run(`kallisto index --index $id $tr`)

        println("\nMade kallisto index at $id\n")

    end

    if ger2 !== nothing

        println("Running paired end psuedoalignment")

        run(`kallisto quant --threads $n_jo --index $id --output-dir $ou $ger1 $ger2`)

    elseif ger2 === nothing

        println("Running single end psuedoalignment")

        run(
            `kallisto quant --single --fragment-length $fr --sd $sd --threads $n_jo --index $id --output-dir $ou $ger1`,
        )

    end

    return

end
