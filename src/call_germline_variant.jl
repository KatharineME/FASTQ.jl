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
)::Nothing

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    if !ispath(pao)

        mkdir(pao)

    end


    # Set config parameters

    vot = basename(to)

    pab = dirname(abspath(ba))

    vob = basename(pab)

    pag = dirname(abspath(fa))

    vog = basename(pag)
    
    vof = joinpath(vog, basename(fa))

    voc = joinpath(vog, "chromosome", basename(chs))

    voba = joinpath(vob, basename(ba))

    pao = abspath(pao)

    voo = basename(pao)
    
    pam = joinpath(voo, "manta")

    pamr = joinpath(pam, "runWorkflow.py")
    
    
    co = "--referenceFasta /home/$vof --callRegions home/$voc --bam home/$voba"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end


    # Set run parameters

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"

    pav = joinpath("results", "variants")



    # Run the docker container

    id = readlines(pipeline(
                 `docker run --interactive --detach --tty --user root --memory=30g --volume $to:/home/$vot --volume $pab:/home/$vob --volume $pag:/home/$vog --volume $pao:/home/$voo centos:centos6 bash`,
               )
       )


    # Configure and run manta

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

    return

    

    ## Bcftools

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
