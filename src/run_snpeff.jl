function run_snpeff(pao, me, pas, paco, n_jo)

    sn = joinpath(pao, "snpeff")

    snvc = joinpath(sn, "snpeff.vcf.gz")

    mkpath(sn)

    run(
        pipeline(
            `java -Xmx$(me)g -jar $pas GRCh38.99 -noLog -verbose -csvStats $(joinpath(sn, "stats.csv")) -htmlStats $(joinpath(sn, "stats.html")) $paco`,
            `bgzip --threads $n_jo --stdout`,
            snvc,
        ),
    )

    run(`tabix $snvc`)

    ps = joinpath(pao, "pass.vcf.gz")

    run(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $snvc`,
            `bgzip --threads $n_jo --stdout`,
            ps,
        ),
    )

    run(`tabix $ps`)

    return nothing

end
