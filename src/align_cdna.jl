function align_cdna(pr, fq1, fq2, fa, n_jo)

    pa = dirname(pr)

    @assert make_directory(pa, "align cdna")

    ge = joinpath(dirname(fa), "star_indexes")

    if !ispath(ge)

        mkdir(ge)

        println("\nMaking STAR indices, this may take a while...\n")

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $ge --genomeFastaFiles $fa`,
        )

    end

    println("\nRunning STAR...\n")

    run(
        `star --runThreadN $n_jo --genomeDir $ge --readFilesIn $fq1 $fq2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = string(pr, "Aligned.sortedByCoord.out.bam")

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    println("\ncDNA Alignment finished\n")

    return

end
