function align_cdna(
    pr::String,
    fq1::String, #read1
    fq2::String, #read2
    fa::String, #reference
    n_jo::Int64,
    me::Int64, #memory
)::Nothing

    pa = dirname(pr)

    if check_directory(pa, "align cdna")

        return nothing

    end

    ge = joinpath(dirname(fa), "star_indexes")

    if !ispath(ge)

        mkdir(ge)

        println("\nMaking STAR indices (this may take a while)...\n")

        run(`star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $ge --genomeFastaFiles $fa`)

    end

    println("\nRunning STAR...\n")

    run(`star --runThreadN $n_jo --genomeDir $ge --readFilesIn $fq1 $fq2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`)
   
    ba = string(pr, "Aligned.sortedByCoord.out.bam")

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    println("\ncDNA Alignment finished\n")

    return nothing

end
