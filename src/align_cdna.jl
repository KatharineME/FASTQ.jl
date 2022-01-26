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

        run(`star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $ge --genomeFastaFiles $fa`)

    end

    run(`star --runThreadN $n_jo --genomeDir $ge --readFilesIn $fq1 $fq2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $sa`)

    return nothing

end
