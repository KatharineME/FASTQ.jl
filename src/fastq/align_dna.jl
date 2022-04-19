function align_dna(al, sa, ba, r1, r2, ge, n_jo, me)

    Fastq.support.log()

    Fastq.support.error_if_directory(al)

    gei = "$ge.mmi"

    if !ispath(gei)

        run(`minimap2 -t $n_jo -d $gei $ge`)

    end

    tm = joinpath(al, "samtools_sort")

    du = joinpath(al, "$sa.unmarked_duplicates.bam")

    run(
        pipeline(
            `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" $gei $r1 $r2`,
            `samtools fixmate --threads $n_jo -u -m - -`,
            `samtools sort --threads $n_jo -T $tm -o $du`,
        ),
    )

    run(`samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $du $ba`)

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    rm(du)

end
