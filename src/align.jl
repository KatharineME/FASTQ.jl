function align(
    mo::String,
    sa::String,
    fq1::String,
    fq2::String,
    fa::String,
    ba::String,
    n_jo::Int64,
    me::Int64,
)::Nothing


    fai = string(fa, ".mmi")

    if !ispath(fai)

        run(`minimap2 -t $n_jo -d $fai $fa`)

    end

    di = splitdir(ba)[1]

    mkpath(di)

    if mo == "dna"

        md = "sr"


    elseif mo == "cdna"

        md = "splice -uf"

    end

    run(
        pipeline(
            `minimap2 -ax $md -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" -a $fai $fq1 $fq2`,

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


    return nothing

end

export align
