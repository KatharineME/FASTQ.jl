function align_dna(
    sa,
    fq1, #read1
    fq2, #read2
    fa, #reference
    ba, #final bam file
    n_jo,
    me, #memory
)

    di = dirname(ba)

    @assert make_directory(di, "align dna")

    fai = string(fa, ".mmi")

    if !ispath(fai)

        run(`minimap2 -t $n_jo -d $fai $fa`)

    end

    tm = joinpath(di, "samtools_sort.bam")

    run(
        pipeline(
            `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" $fai $fq1 $fq2`,
            `samtools fixmate --threads $n_jo -u -m - -`,
            `samtools sort --threads $n_jo -T $tm -u -`,
            "$ba",
        ),
    )

    run(`samtools markdup --threads $n_jo --reference $fa --output-fmt BAM $ba final.bam`)

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    println("\nDNA Alignment finished\n")

    return

end
