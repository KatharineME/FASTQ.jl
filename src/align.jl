function align(
    mo::String, #cdna or dna
    sa::String, #
    fq1::String, #read1
    fq2::String, #read2
    fa::String, #reference
    ba::String, #bam directory
    n_jo::Int64,
    me::Int64, #memory
)::Nothing


    fai = string(fa, ".mmi")

    if !ispath(fai)

        run(`minimap2 -t $n_jo -d $fai $fa`)

    end

    di = splitdir(ba)[1]

    mkpath(di)

    if mo == "dna"

        run(
            pipeline(
                `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" -a $fai $fq1 $fq2`,

                `samtools fixmate --threads $n_jo -m - -`,

                `samtools sort --threads $n_jo`,

                "$ba.tmp",
            ),
        )

    ### minimap2 -ax sr -t 10 -K 4G -R "@RG\tID:$test\tSM:$test" /home/jovyan/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.mmi test_dna_4k.R1.fastq.gz test_dna_4k.R2.fastq.gz | samtools fixmate -u -m - - | samtools sort -u -@2 -T ./tmp/test_run_dna_tmp - | samtools markdup -@8 --reference /home/jovyan/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz -O BAM - test_dna_final.bam


        run(`samtools markdup --threads $n_jo -s $ba.tmp $ba`)

        rm("$ba.tmp")

        run(`samtools index -@ $n_jo $ba`)

        run(
            pipeline(`samtools flagstat --threads $n_jo $ba`, "$ba.flagstat"),
        )

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
