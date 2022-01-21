function call_variant(
    mo::String,
    ge::Union{String, Nothing},
    so::Union{String, Nothing},
    ta::Bool,
    fa::String,
    chs::String,
    chn::String,
    pao::String,
    n_jo::Int,
    me::Int,
    pas::String,
)

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    # Set config parameters

    co = "--referenceFasta $fa --callRegions $chs"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end

    if ge != nothing && so != nothing

        co = "$co --normalBam $ge --tumorBam $so"

    elseif ge != nothing

        co = "$co --bam $ge"

    else

        error(
            "specify germ and soma .bam (for processing soma) or germ .bam (for processing germ).",
        )

    end

    # Set run parameters

    ruString = "--mode local --jobs $n_jo --memGb $me --quiet"

    pav = joinpath("results", "variants")

    if mo == "dna"

        pam = joinpath(pao, "manta")

        pamr = joinpath(pam, "runWorkflow.py")

        run(
            `bash -c "source activate py2 && configManta.py $co --outputContig --runDir $pam && $pamr $ru"`,
        )

        # ./configManta.py --bam /home/craft/data/738/germ.bam --referenceFasta /home/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz --callRegions /home/craft/data/grch/chromosome/chromosome.bed.gz --outputContig --runDir /home/craft/tool/Fastq.jl/tool/manta_config/

    end

    pas = joinpath(pao, "strelka")


    # Configure strelka

    if ge != nothing && so != nothing

        st = "configureStrelkaSomaticWorkflow.py $co --indelCandidates $(joinpath(pam, pav, "candidateSmallIndels.vcf.gz")) --runDir $pas"

    else

        st = "configureStrelkaGermlineWorkflow.py $co --runDir $pas"

        # ./configureStrelkaGermlineWorkflow.py --referenceFasta /home/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz --callRegions /home/craft/data/grch/chromosome/chromosome.bed.gz --bam /home/craft/data/738/germ.bam --runDir /home/craft/tool/Fastq.jl/tool/strelka_config/

    end


    # Run strelka

    pasr = joinpath(pas, "runWorkflow.py")

    run(`bash -c "source activate py2 && $st && $pasr $ru"`)


    if ge != nothing && so != nothing

        sa = joinpath(pao, "sample.txt")

        # TODO: get sample names (maybe from .bam) and use them instead of "Germ" and "Soma"

        open(io -> write(io, "Germ\nSoma"), sa; write = true)

        pain = joinpath(pas, pav, "somatic.indels.vcf.gz")

        run(
            pipeline(
                `bcftools reheader --threads $n_jo --samples $sa $pain`,
                "$pain.tmp",
            ),
        )

        mv("$pain.tmp", pain; force = true)

        run(`tabix --force $pain`)

        pasv = joinpath(pas, pav, "somatic.snvs.vcf.gz")

        run(
            pipeline(
                `bcftools reheader --threads $n_jo --samples $sa $pasv`,
                "pasv.tmp",
            ),
        )

        mv("$pasv.tmp", pasv; force = true)

        run(`tabix --force $pasv`)

        vc_ = [joinpath(pam, pav, "somaticSV.vcf.gz"), pain, pasv]

    elseif mo == "cdna"

        vc_ = [joinpath(pas, pav, "variants.vcf.gz")]

        println("this is vc_: $vc_")

    else

        vc_ = [
            joinpath(pam, pav, "diploidSV.vcf.gz"),
            joinpath(pas, pav, "variants.vcf.gz"),
        ]

    end

    paco = joinpath(pao, "concat.vcf.gz")

    run(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chn`,
            `bgzip --threads $n_jo --stdout`,
            paco,
        ),
    )

    run(`tabix $paco`)

    sn = joinpath(pao, "snpeff")

    mkpath(sn)

    snvc = joinpath(sn, "snpeff.vcf.gz")

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

    return run(`tabix $ps`)

end
