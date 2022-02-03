function align_cdna(
    sa::String,
    fq1::String, #read1
    fq2::String, #read2
    fa::String, #reference
    ba::String, #bam directory
    n_jo::Int64,
    me::Int64, #memory
)::Nothing

    ge = joinpath(dirname(fa), "star_indexes")

    if !ispath(ge)

        mkdir(ge)

        println("\nMaking STAR indices (this may take a while)...\n")

        run(`star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $ge --genomeFastaFiles $fa`)

    end

    println("\nRunning STAR...\n")

    run(`star --runThreadN $n_jo --genomeDir $ge --readFilesIn $fq1 $fq2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $sa`)

    println("\ncDNA Alignment finished\n")

    return nothing

end
