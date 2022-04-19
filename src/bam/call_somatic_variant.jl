function call_somatic_variant(ta, ge, so, fa, chs, chn, pao, n_jo, me, to, sn)

    Fastq.support.log()

    Fastq.support.index_genome_files(fa, chs)

    Fastq.support.error_if_directory(pao)


    # Run docker container

    id, voo, vof, voc, vogefi, vosofi, vot = run_docker_container(to, fa, chs, ge, pao, so)


    # Set config parameters

    co = "--referenceFasta /home/$vof --callRegions /home/$voc --normalBam /home/$vogefi --tumorBam /home/$vosofi"

    if ta

        co = "$co --exome"

    end

    # Set runtime parameters

    ru = set_strelka_manta_run(n_jo, me)


    # Configure and run manta

    vom = configure_and_run_manta(voo, id, vot, co, ru)


    # Configure and run strelka

    pav = joinpath("results", "variants")

    vost = joinpath(voo, "strelka")

    vosr = joinpath(vost, "runWorkflow.py")

    sc = "$(Fastq.STRELKA)/bin/configureStrelkaSomaticWorkflow.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --indelCandidates $(joinpath("home", vom, pav, "candidateSmallIndels.vcf.gz")) --runDir /home/$vost && ./home/$vosr $ru"`,
        ),
    )

    println("$(join(re, " "))\n")


    # Remove docker container

    remove_docker_container(id)


    # Combine vcfs

    past = joinpath(pao, "strelka")

    sa = joinpath(pao, "sample.txt")

    open(io -> write(io, "Germline\nSomatic"), sa; write = true)

    ie = joinpath(past, pav, "somatic.indels.vcf.gz")

    ier = Fastq.vcf.reheader_vcf(sa, ie, n_jo)

    sv = joinpath(past, pav, "somatic.snvs.vcf.gz")

    svr = Fastq.vcf.reheader_vcf(sa, sv, n_jo)

    svm = joinpath(pao, "manta", pav, "somaticSV.vcf.gz")

    svmr = Fastq.vcf.reheader_vcf(sa, svm, n_jo)

    vc_ = [ier, svr, svmr]

    paco = joinpath(pao, "concat.vcf.gz")

    Fastq.vcf.combine_vcf(vc_, chn, paco, n_jo)

    run(`tabix $paco`)

    # snpeff

    run_snpeff(pao, me, sn, paco, n_jo)

end
