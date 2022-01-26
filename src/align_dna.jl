function align_dna(
    sa::String,
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

    tm = joinpath(ba, "samtools_sort_temp")

    run(
        pipeline(
            `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" -a $fai $fq1 $fq2`,

            `samtools fixmate --threads $n_jo -u -m - -`,

            `samtools sort --threads $n_jo -u -T $tm`,

            `samtools markdup --threads $n_jo --reference $fa --output-fmt BAM $ba`
           )
       )

    run(`samtools index -@ $n_jo $ba`)

    run(`samtools stats --threads $n_jo $ba`, "$ba.stat")

    return nothing

end
