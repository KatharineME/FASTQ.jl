function align_dna(al, sa, ba, fq1, fq2, ge, n_jo, me)

    @assert make_directory(al, "align dna")

    gei = string(ge, ".mmi")

    if !ispath(gei)

        run(`minimap2 -t $n_jo -d $gei $ge`)

    end

    tm = joinpath(al, "samtools_sort.bam")

    run(
        pipeline(
            `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" $gei $fq1 $fq2`,
            `samtools fixmate --threads $n_jo -u -m - -`,
            `samtools sort --threads $n_jo -T $tm -u -`,
            "$ba",
        ),
    )

    ma = "$ba.markdup.bam"

    run(`samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $ba $ma`)

    run(`samtools index -@ $n_jo $ma`)

    run(pipeline(`samtools stats --threads $n_jo $ma`, "$ma.stat"))

    println("\nDNA Alignment finished\n")

    return

end
