module VCF

using Fastq

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

    run(
        pipeline(
            `java -Xmx$(me)g -jar $sn GRCh37.75 -noLog -verbose -csvStats $(joinpath(pasn, "stats.csv")) -htmlStats $(joinpath(pasn, "stats.html")) $paco`,
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

function combine_vcf(vc_, chn, paco, n_jo)

    run(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chn`,
            `bgzip --threads $n_jo --stdout`,
            paco,
        ),
    )

end

function reheader_vcf(sa, pa, n_jo)

    na = "$(split(basename(pa), "vcf.gz")[1]) reheader.vcf.gz"

    par = joinpath(dirname(pa), na)

    run(pipeline(`bcftools reheader --threads $n_jo --samples $sa $pa`, "$par"))

    run(`tabix --force $par`)

    return par

end

end
