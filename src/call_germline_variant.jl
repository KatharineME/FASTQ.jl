function call_germline_variant(
    mo::String,
    ge::String,
    ta::Bool,
    fa::String,
    chs::String,
    chn::String,
    pao::String,
    n_jo::Int,
    me::Int,
    to::String, #path to tools: strelka and manta
    pas::String, #path to snpeff
)::Nothing

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    if check_directory(pao, "call germline variant")

        return nothing

    end

    # Run docker container

    id, voo, vof, voc, voge, vot = run_docker_container(to, fa, chs, ge, pao)

    

    # Set config parameters
    
    co = "--referenceFasta /home/$vof --callRegions home/$voc --bam home/$voge"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end



    # Set run parameters

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"



    # Configure and run manta

    pam = joinpath(voo, "manta")

    pamr = joinpath(pam, "runWorkflow.py")
    
    sc = "manta-1.6.0.centos6_x86_64/bin/configManta.py" 

    re =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$pam && ./home/$pamr $ru"`))

    println("$(join(re, " "))\n")

   

    # Configure and run strelka

    past = joinpath(voo, "strelka")
    
    pasr = joinpath(past, "runWorkflow.py")

    sc = "strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py"

    re =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --runDir /home/$past && ./home/$pasr $ru"`))

    println("$(join(re, " "))\n")
    
    

    # Remove docker container

    remove_docker_container(id)



    ## Bcftools

    pav = joinpath("results", "variants")
    
    if mo == "cdna"

        vc_ = [joinpath(past, pav, "variants.vcf.gz")]

    else

        vc_ = [
            joinpath(pam, pav, "diploidSV.vcf.gz"),
            joinpath(past, pav, "variants.vcf.gz"),
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


    # Snpeff

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
