function align_cdna(al, sa, r1, r2, ge, n_jo)

    Fastq.support.error_if_directory(al)

    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        println("\nMaking STAR indices, this may take a while...\n")

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ge`,
        )

    end

    #log("Running STAR")

    pr = joinpath(al, "$(sa).")

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    println("\ncDNA Alignment finished\n")

end

#function log(io, st)
#
#    me = "($time) $st" 
#
#    write(io, st)
#
#    println("="^99)
#
#    println(st)
#
#    println("="^99)
#
#end
