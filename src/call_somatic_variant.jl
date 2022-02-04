function call_somatic_variant(
    mo::String,
    ge::String,
    so::String,
    fa::String,
    chs::String,
    chn::String,
    pao::String,
    n_jo::Int,
    me::Int,
    to::String,
    pas::String,
)::Nothing

    index_genome_files(fa, chs)

    if check_directory(pao, "call somatic variant")

        return nothing

    end


    # Run docker container
   
    id, voo, vof, voc, vogefi, vosofi, vot = run_docker_container(to, fa, chs, ge, pao, so)


    # Set config parameters

    co = "--referenceFasta /home/$vof --callRegions /home/$voc --normalBam /home/$vogefi --tumorBam /home/$vosofi"


    # Set run parameters

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"



    # Configure and run manta

    pam = configure_and_run_manta(voo, id, vot, co, ru)

    println(pam)



    # Configure and run strelka

    pav = joinpath("results", "variants")

    past = joinpath(voo, "strelka")
    
    pasr = joinpath(past, "runWorkflow.py")

    sc = "strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"

    re =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --indelCandidates $(joinpath("home", pam, pav, "candidateSmallIndels.vcf.gz")) --runDir /home/$past && ./home/$pasr $ru"`))

    println("$(join(re, " "))\n")

   

    # Remove docker container

    remove_docker_container(id)


    # bcftools

    sa = joinpath(pao, "sample.txt")
    
    open(io -> write(io, "Germline\nSomatic"), sa; write = true)

    ie = joinpath(past, pav, "somatic.indels.vcf.gz")

    ier = reheader_vcf(sa, ie, n_jo)

    sv = joinpath(past, pav, "somatic.snvs.vcf.gz")

    svr = reheader_vcf(sa, sv, n_jo)

    svm = joinpath(pam, pav, "somaticSV.vcf.gz")

    svmr = reheader_vcf(sa, svm, n_jo)

    vc_ = [ier, svr, svmr]

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
