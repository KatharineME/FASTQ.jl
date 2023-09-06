module VCF

using FASTQ

function annotate_with_snpeff(pa, me, sn, pac, n_jo)

    FASTQ.Support.log_sub_level_function()

    pas = joinpath(pa, "snpeff")

    FASTQ.Support.trash_remake_directory(pas)

    stc = joinpath(pas, "stats.csv")

    sth = joinpath(pas, "stats.html")

    vc = joinpath(pas, "snpeff.vcf.gz")

    run(
        pipeline(
            `java -Xmx$(me)g -jar $sn GRCh38.99 -noLog -verbose -csvStats $stc -htmlStats $sth $pac`,
            `bgzip --threads $n_jo --stdout`,
            vc,
        ),
    )

    run(`tabix $vc`)

    pap = joinpath(pa, "pass.vcf.gz")

    run(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $vc`,
            `bgzip --threads $n_jo --stdout`,
            pap,
        ),
    )

    run(`tabix $pap`)

    pap

end

function annotate_with_snpsift(pa, sn, va, pap, n_jo)

    FASTQ.Support.log_sub_level_function()

    pas = joinpath(pa, "snpsift")

    FASTQ.Support.trash_remake_directory(pas)

    vc = joinpath(pas, "snpsift.vcf.gz")

    ss = joinpath(dirname(sn), "SnpSift.jar")

    run(
        pipeline(
            `java -jar $ss annotate -tabix -v $va $pap`,
            `bgzip --threads $n_jo --stdout`,
            vc,
        ),
    )

    run(`tabix $vc`)

end

function combine_vcf(pa, n_jo, vc_, chn)

    FASTQ.Support.log_sub_level_function()

    run(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chn`,
            `bgzip --threads $n_jo --stdout`,
            pa,
        ),
    )

end

function reheader_vcf(pa, n_jo, sa)

    FASTQ.Support.log_sub_level_function()

    na = split(basename(pa), "vcf.gz")[1]

    par = joinpath(dirname(pa), "$na reheader.vcf.gz")

    run(pipeline(`bcftools reheader --threads $n_jo --samples $sa $pa`, "$par"))

    run(`tabix --force $par`)

    par

end

end
