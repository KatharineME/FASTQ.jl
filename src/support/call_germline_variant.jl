function call_germline_variant(mo, ta, ge, fa, chs, chn, pao, n_jo, me, to, sn)

    index_genome_files(fa, chs)

    @assert make_directory(pao, "call germline variant")


    # Run docker container

    id, voo, vof, voc, vogefi, vot = run_docker_container(to, fa, chs, ge, pao, nothing)


    # Set config parameters

    co = "--referenceFasta /home/$vof --callRegions home/$voc --bam home/$vogefi"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end


    # Configure and run manta

    configure_and_run_manta(voo, id, vot, co, RU)


    # Configure and run strelka

    vost = joinpath(voo, "strelka")

    vosr = joinpath(vost, "runWorkflow.py")

    sc = "$STRELKA/bin/configureStrelkaGermlineWorkflow.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --runDir /home/$vost && ./home/$vosr $RU"`,
        ),
    )

    println("$(join(re, " "))\n")


    # Remove docker container

    remove_docker_container(id)


    ## bcftools

    pav = joinpath("results", "variants")

    if mo == "cdna"

        vc_ = [joinpath(pao, "strelka", pav, "variants.vcf.gz")]

    else

        vc_ = [
            joinpath(pao, "manta", pav, "diploidSV.vcf.gz"),
            joinpath(pao, "strelka", pav, "variants.vcf.gz"),
        ]

    end

    paco = joinpath(pao, "concat.vcf.gz")

    combine_vcf(vc_, chn, paco, n_jo)

    run(`tabix $paco`)


    # snpeff

    run_snpeff(pao, me, sn, paco, n_jo)

    return

end
