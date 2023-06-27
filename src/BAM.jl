module BAM

using FASTQ

function _configure_and_run_manta(voo, id, vot, co, ru)

    vom = joinpath(voo, "manta")

    sc = "$(FASTQ.MANTA)/bin/configManta.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$vom && ./home/$(joinpath(vom, "runWorkflow.py")) $ru"`,
        ),
    )

    println("$(join(re, " "))\n")

    vom

end

function _set_strelka_manta_run(n_jo, me)

    return "--mode local --jobs $n_jo --memGb $me --quiet"

end

function _run_strelka_manta_docker_container(to, fa, chs, ge, pao; so = nothing)

    vot = basename(to)

    page = dirname(BioLab.Path.make_absolute(ge))

    sp = split(page, "/")

    voge = joinpath(sp[length(sp) - 1], last(sp))

    vogefi = joinpath(voge, basename(ge))

    pag = dirname(BioLab.Path.make_absolute(fa))

    vog = basename(pag)

    vof = joinpath(vog, basename(fa))

    voc = joinpath(vog, "chromosome", basename(chs))

    pao = BioLab.Path.make_absolute(pao)

    voo = basename(pao)

    com = ["docker", "run", "--interactive", "--detach", "--tty", "--user", "root", "--memory=30g"]

    vo = [
        "--volume",
        "$to:/home/$vot",
        "--volume",
        "$page:/home/$voge",
        "--volume",
        "$pag:/home/$vog",
        "--volume",
        "$pao:/home/$voo",
    ]

    con = ["centos:centos6", "bash"]

    if so !== nothing

        paso = dirname(BioLab.Path.make_absolute(so))

        voso = basename(paso)

        vosofi = joinpath(voso, basename(so))

        id = readlines(pipeline(`$com $vo --volume $paso:/home/$voso $con`))

        return id, voo, vof, voc, vogefi, vosofi, vot

    else

        id = readlines(pipeline(`$com $vo $con`))

        return id, voo, vof, voc, vogefi, vot

    end

end

function call_germline_variant(fa, chs, pao, to, ge, ta, mo, n_jo, me, chn, sn, rs, va)

    FASTQ.Support.log()

    FASTQ.Support.index_genome_files(fa, chs)

    FASTQ.Support.error_if_directory(pao)

    id, voo, vof, voc, vogefi, vot = _run_strelka_manta_docker_container(to, fa, chs, ge, pao)

    co = "--referenceFasta /home/$vof --callRegions home/$voc --bam home/$vogefi"

    if ta

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

    end

    ru = _set_strelka_manta_run(n_jo, me)

    _configure_and_run_manta(voo, id, vot, co, ru)

    vost = joinpath(voo, "strelka")

    sc = "$(FASTQ.STRELKA)/bin/configureStrelkaGermlineWorkflow.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --runDir /home/$vost && ./home/$(joinpath(vost, "runWorkflow.py")) $ru"`,
        ),
    )

    println("$(join(re, " "))\n")

    FASTQ.Support.remove_docker_container(id)

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

    FASTQ.VCF.combine_vcf(vc_, chn, paco, n_jo)

    run(`tabix $paco`)

    papa = FASTQ.VCF.annotate_with_snpeff(pao, me, sn, paco, n_jo)

    if rs

        FASTQ.VCF.annotate_with_snpsift(pao, sn, va, papa, n_jo)

    end

end

function call_somatic_variant(ta, ge, so, fa, chs, chn, pao, n_jo, me, to, sn, rs, va)

    FASTQ.Support.log()

    FASTQ.Support.index_genome_files(fa, chs)

    FASTQ.Support.error_if_directory(pao)


    # Run docker container

    id, voo, vof, voc, vogefi, vosofi, vot =
        _run_strelka_manta_docker_container(to, fa, chs, ge, pao, so = so)


    # Set config parameters

    co = "--referenceFasta /home/$vof --callRegions /home/$voc --normalBam /home/$vogefi --tumorBam /home/$vosofi"

    if ta

        co = "$co --exome"

    end

    # Set runtime parameters

    ru = _set_strelka_manta_run(n_jo, me)


    # Configure and run manta

    vom = _configure_and_run_manta(voo, id, vot, co, ru)


    # Configure and run strelka

    pav = joinpath("results", "variants")

    vost = joinpath(voo, "strelka")

    sc = "$(FASTQ.STRELKA)/bin/configureStrelkaSomaticWorkflow.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --indelCandidates $(joinpath("home", vom, pav, "candidateSmallIndels.vcf.gz")) --runDir /home/$vost && ./home/$(joinpath(vost, "runWorkflow.py")) $ru"`,
        ),
    )

    println("$(join(re, " "))\n")


    # Remove docker container

    FASTQ.Support.remove_docker_container(id)


    # Combine vcfs

    past = joinpath(pao, "strelka")

    sa = joinpath(pao, "sample.txt")

    open(io -> write(io, "Germline\nSomatic"), sa; write = true)

    ie = joinpath(past, pav, "somatic.indels.vcf.gz")

    ier = FASTQ.VCF.reheader_vcf(sa, ie, n_jo)

    sv = joinpath(past, pav, "somatic.snvs.vcf.gz")

    svr = FASTQ.VCF.reheader_vcf(sa, sv, n_jo)

    svm = joinpath(pao, "manta", pav, "somaticSV.vcf.gz")

    svmr = FASTQ.VCF.reheader_vcf(sa, svm, n_jo)

    vc_ = [ier, svr, svmr]

    paco = joinpath(pao, "concat.vcf.gz")

    FASTQ.VCF.combine_vcf(vc_, chn, paco, n_jo)

    run(`tabix $paco`)

    # snpeff

    papa = FASTQ.VCF.annotate_with_snpeff(pao, me, sn, paco, n_jo)


    # snpsift

    if rs

        FASTQ.VCF.annotate_with_snpsift(pao, sn, va, papa, n_jo)

    end

end

end
