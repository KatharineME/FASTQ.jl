function call_germline_variant(
    mo::String,
    ba::Union{String, Nothing},
    ta::Bool,
    fa::String,
    chs::String,
    chn::String,
    pao::String,
    n_jo::Int,
    me::Int,
    to::String, #path to tools: strelka and manta
    pas::String, #path to snpeff
)

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end


    # Set config parameters

    co = "--referenceFasta $fa --callRegions $chs --bam $ba"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end


    # Set run parameters

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"

    pav = joinpath("results", "variants")


    # Configure manta

    pam = joinpath(pao, "manta")

    pamr = joinpath(pam, "runWorkflow.py")

    vo = last(split(to, "/"))
    
    id = run_docker_container(to, vo)




    
    re =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`))

    println("$(join(re, " "))\n")

    run(
        `bash -c "source activate py2 && configManta.py $co --outputContig --runDir $pam && $pamr $ru"`,
        )


    # Configure strelka

    past = joinpath(pao, "strelka")
    
    st = "configureStrelkaGermlineWorkflow.py $co --runDir $past"


    # Run strelka

    pasr = joinpath(past, "runWorkflow.py")

    run(`bash -c "source activate py2 && $st && $pasr $ru"`)

    if mo == "cdna"

        vc_ = [joinpath(pas, pav, "variants.vcf.gz")]

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


    # SNPEFF

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
