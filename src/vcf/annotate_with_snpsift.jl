function annotate_with_snpsift(pa, sn)

    pass = joinpath(pa, "snpsift")

    vc = joinpath(pass, "snpsift.vcf.gz")

    mkpath(pass)

    run(
        pipeline(
            `java -jar ~/craft/guardiome/tool/snpEff/SnpSift.jar annotate -tabix -id -v ~/craft/guardiome/tool/ensembl/homo_sapiens-chr1_y.vcf.gz pass.vcf.gz > pass.rsids.vcf`,
            vc,
        ),
    )

    run(`tabix $vc`)

end
