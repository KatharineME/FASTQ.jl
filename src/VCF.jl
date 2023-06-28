module VCF

using FASTQ

function annotate_with_snpeff(pao, me, sn, paco, n_jo)

    FASTQ.Support.log()

    pasn = joinpath(pao, "snpeff")

    FASTQ.Support.error_if_directory(pasn)

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

    papa

end

function annotate_with_snpsift(pa, sn, va, papa, n_jo)

    FASTQ.Support.log()

    pass = joinpath(pa, "snpsift")

    FASTQ.Support.error_if_directory(pass)

    vc = joinpath(pass, "snpsift.vcf.gz")

    run(
        pipeline(
            `java -jar $(joinpath(dirname(sn), "SnpSift.jar")) annotate -tabix -id -v $va $papa`,
            `bgzip --threads $n_jo --stdout`,
            vc,
        ),
    )

    run(`tabix $vc`)

end

function combine_vcf(n_jo, vc_, chn, paco)

    run(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chn`,
            `bgzip --threads $n_jo --stdout`,
            paco,
        ),
    )

end

function reheader_vcf(pa, n_jo, sa)

    par = joinpath(dirname(pa), "$(split(basename(pa), "vcf.gz")[1]) reheader.vcf.gz")

    run(pipeline(`bcftools reheader --threads $n_jo --samples $sa $pa`, "$par"))

    run(`tabix --force $par`)

    par

end

end
