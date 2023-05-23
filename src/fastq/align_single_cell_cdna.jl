function align_cdna(al, sa, r1, r2, ge, n_jo)

    Fastq.support.log()

    Fastq.support.error_if_directory(al)

    # Change to star index made with gtf gene annotation file
    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        println("\nMaking STAR indices, this may take a while\n")

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(al, "$(sa).")

    # Update star run
    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end
