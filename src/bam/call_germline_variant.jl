function call_germline_variant(mo, ta, ge, fa, chs, chn, pao, n_jo, me, to, sn)

    Fastq.support.log()

    Fastq.support.index_genome_files(fa, chs)

    Fastq.support.error_if_directory(pao)


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


    # Set runtime parameters

    ru = set_strelka_manta_run(n_jo, me)


    # Configure and run manta

    configure_and_run_manta(voo, id, vot, co, ru)


    # Configure and run strelka

    vost = joinpath(voo, "strelka")

    vosr = joinpath(vost, "runWorkflow.py")

    sc = "$(Fastq.STRELKA)/bin/configureStrelkaGermlineWorkflow.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --runDir /home/$vost && ./home/$vosr $ru"`,
        ),
    )

    println("$(join(re, " "))\n")


    # Remove docker container

    remove_docker_container(id)


    ## Combine vcfs

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

    Fastq.vcf.combine_vcf(vc_, chn, paco, n_jo)

    run(`tabix $paco`)


    # snpeff

    Fastq.vcf.run_snpeff(pao, me, sn, paco, n_jo)

end
