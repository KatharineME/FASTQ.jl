function psuedoalign(tr, n_jo, ou, r1, r2, fr, sd)

    pa = joinpath(ou, "psuedoalign")

    Fastq.support.error_if_directory(pa)

    id = "$tr.kallisto_index"

    if !ispath(id)

        println("\nCreating kallisto index...")

        run(`kallisto index --index $id $tr`)

        println("\nMade kallisto index at $id\n")

    end

    fu = ["kallisto", "quant"]

    ru = ["--threads", "$n_jo", "--index", "$id", "--output-dir", "$pa"]

    if r2 !== nothing

        println("Running paired end psuedoalignment")

        run(`$fu $ru $r1 $r2`)

    else

        println("Running single end psuedoalignment")

        run(`$fu --single --fragment-length $fr --sd $sd $ru $r1`)

    end

    return

end
