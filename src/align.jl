function align(
    mo::String, #cdna or dna
    sa::String, #
    fq1::String, #read1
    fq2::String, #read2
    fa::String, #reference
    cr::String, #bam directory
    n_jo::Int64,
    me::Int64, #memory
)::Nothing


    fai = string(fa, ".mmi")

    if !ispath(fai)

        run(`minimap2 -t $n_jo -d $fai $fa`)

    end

    di = splitdir(cr)[1]

    mkpath(di)

    if mo == "dna"

        tm = joinpath(cr, "samtools_sort_temp")

        run(
            pipeline(
                `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" -a $fai $fq1 $fq2`,

                `samtools fixmate --threads $n_jo -u -m - -`,

                `samtools sort --threads $n_jo -u -T $tm`,

                `samtools markdup --threads $n_jo --reference $fa --output-fmt CRAM $cr`
               )
           )

        run(`samtools index -@ $n_jo $cr`)

        run(`samtools stats --threads $n_jo $cr`, "$cr.stat")

    elseif mo == "cdna"

        # Directory where STAR will generate its genome indexes 
        ge = joinpath(dirname(fa), "star_indexes")

        if !ispath(ge)

            mkdir(ge)

            # Generate genome indexes for STAR
            run(`star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $ge --genomeFastaFiles $fa`)

        end

        # Run STAR
        run(`star --runThreadN $n_jo --genomeDir $ge --readFilesIn $fq1 $fq2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $sa`)

    end

    return nothing

end
