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

    index_genome_files(fa, chs)

    if check_directory(pao, "call germline variant")

        return nothing

    end



    # Run docker container

    id, voo, vof, voc, vogefi, vot = run_docker_container(to, fa, chs, ge, pao)

    

    # Set config parameters
    
    co = "--referenceFasta /home/$vof --callRegions home/$voc --bam home/$vogefi"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end



    # Set run parameters

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"



    # Configure and run manta
    
    pam = configure_and_run_manta(pao, id, vot, co, ru)


    # Configure and run strelka

    past = joinpath(voo, "strelka")
    
    pasr = joinpath(past, "runWorkflow.py")

    sc = "strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py"

    re =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --runDir /home/$past && ./home/$pasr $ru"`))

    println("$(join(re, " "))\n")
    
    

    # Remove docker container

    remove_docker_container(id)



    ## bcftools

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


    # snpeff

    run_snpeff(pao, me, pas, paco, n_jo)

    return nothing 

end
