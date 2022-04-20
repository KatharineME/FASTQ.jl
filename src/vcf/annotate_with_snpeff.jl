function annotate_with_snpeff(pao, me, sn, paco, n_jo)

    Fastq.support.log()

    pasn = joinpath(pao, "snpeff")

    Fastq.support.error_if_directory(pasn)

    vc = joinpath(pasn, "snpeff.vcf.gz")

    run(
        pipeline(
            `java -Xmx$(me)g -jar $sn GRCh38.99 -noLog -verbose -csvStats $(joinpath(pasn, "stats.csv")) -htmlStats $(joinpath(pasn, "stats.html")) $paco`,
            `bgzip --threads $n_jo --stdout`,
            vc,
        ),
    )

    run(`tabix $vc`)

    papa = joinpath(pao, "pass.vcf.gz")

    run(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $vc`,
            `bgzip --threads $n_jo --stdout`,
            papa,
        ),
    )

    run(`tabix $papa`)

    return papa

end
