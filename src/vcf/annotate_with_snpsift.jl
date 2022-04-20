function annotate_with_snpsift(pa, sn, va, papa, n_jo)

    Fastq.support.log()

    pass = joinpath(pa, "snpsift")

    Fastq.support.error_if_directory(pass)

    vc = joinpath(pass, "snpsift.vcf.gz")

    ss = joinpath(dirname(sn), "SnpSift.jar")

    run(
        pipeline(
            `java -jar $ss annotate -tabix -id -v $va $papa`,
            `bgzip --threads $n_jo --stdout`,
            vc,
        ),
    )

    run(`tabix $vc`)

end
