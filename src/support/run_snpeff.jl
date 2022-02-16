function run_snpeff(pao, me, sn, paco, n_jo)

    pasn = joinpath(pao, "snpeff")

    vc = joinpath(pasn, "snpeff.vcf.gz")

    mkpath(pasn)

    run(
        pipeline(
            `java -Xmx$(me)g -jar $sn GRCh38.99 -noLog -verbose -csvStats $(joinpath(pasn, "stats.csv")) -htmlStats $(joinpath(pasn, "stats.html")) $paco`,
            `bgzip --threads $n_jo --stdout`,
            vc,
        ),
    )

    run(`tabix $vc`)

    ps = joinpath(pao, "pass.vcf.gz")

    run(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $vc`,
            `bgzip --threads $n_jo --stdout`,
            ps,
        ),
    )

    run(`tabix $ps`)

    return

end
