module VCF

using ..FASTQ

function combine_vcf(pa, vc_, chn, n_jo)

    FASTQ.Support.log_sub_level_function()

    paco = joinpath(pa, "concat.vcf.gz")

    run(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chn`,
            `bgzip --threads $n_jo --stdout`,
            paco,
        ),
    )

    run(`tabix $paco`)

    paco

end

function reheader_vcf(pa, sa, n_jo)

    FASTQ.Support.log_sub_level_function()

    na = split(basename(pa), "vcf.gz")[1]

    par = joinpath(dirname(pa), "$na.reheader.vcf.gz")

    run(pipeline(`bcftools reheader --threads $n_jo --samples $sa $pa`, "$par"))

    run(`tabix --force $par`)

    par

end

function _check_genome_version(re)

    ren = basename(re)

    ve = Int64

    if occursin("38", ren)

        ve = 38

    elseif occursin("37", ren)

        ve = 37

    else

        @error "Unkown genome version: $re"

    end

    ve

end

function annotate_with_snpeff(pa, vc, re, se, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    ve = _check_genome_version(re)

    sg = "GRCh38.99"

    if ve == 37

        sg = "GRCh37.75"

    end

    stc = joinpath(pa, "stats.csv")

    sth = joinpath(pa, "stats.html")

    vcs = joinpath(pa, "snpeff.vcf.gz")

    run(
        pipeline(
            `java -Xmx$(me)g -jar $se $sg -noLog -verbose -csvStats $stc -htmlStats $sth $vc`,
            `bgzip --threads $n_jo --stdout`,
            vcs,
        ),
    )

    run(`tabix $vcs`)

    vcs

end

function annotate_with_snpsift(pa, vc, va, cl, se, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    s1 = joinpath(pa, "snpsift.vcf.gz")

    ss = joinpath(dirname(se), "SnpSift.jar")

    s2 = joinpath(pa, "clinvar.vcf.gz")

    run(
        pipeline(
            `java -Xmx$(me)g -jar $ss annotate -tabix -v $va $vc`,
            `bgzip --threads $n_jo --stdout`,
            s1,
        ),
    )

    run(
        pipeline(
            `java -Xmx$(me)g -jar $ss annotate -tabix -v $cl $s1`,
            `bgzip --threads $n_jo --stdout`,
            s2,
        ),
    )

    s2

end

function filter_vcf(pa, vc, n_jo)

    FASTQ.Support.log_sub_level_function()

    pap = joinpath(pa, "pass.vcf.gz")

    run(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $vc`,
            `bgzip --threads $n_jo --stdout`,
            pap,
        ),
    )

    run(`tabix -f $pap`)

    pap

end

end
