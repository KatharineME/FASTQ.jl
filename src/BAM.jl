module BAM

using Fastq

function configure_and_run_manta(voo, id, vot, co, ru)

    vom = joinpath(voo, "manta")

    vomr = joinpath(vom, "runWorkflow.py")

    sc = "$(Fastq.MANTA)/bin/configManta.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$vom && ./home/$vomr $ru"`,
        ),
    )

    println("$(join(re, " "))\n")

    return vom

end


####################################################################

function set_strelka_manta_run(n_jo, me)

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"

    return ru

end


####################################################################


function run_strelka_manta_docker_container(to, fa, chs, ge, pao, so)

    vot = basename(to)

    page = dirname(Fastq.support.get_full_path(ge))

    sp = split(page, "/")

    le = length(sp)

    voge = joinpath(sp[le - 1], last(sp))

    vogefi = joinpath(voge, basename(ge))

    pag = dirname(Fastq.support.get_full_path(fa))

    vog = basename(pag)

    vof = joinpath(vog, basename(fa))

    voc = joinpath(vog, "chromosome", basename(chs))

    pao = Fastq.support.get_full_path(pao)

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

        paso = dirname(Fastq.support.get_full_path(so))

        voso = basename(paso)

        vosofi = joinpath(voso, basename(so))

        id = readlines(pipeline(`$com $vo --volume $paso:/home/$voso $con`))

        return id, voo, vof, voc, vogefi, vosofi, vot

    else

        id = readlines(pipeline(`$com $vo $con`))

        return id, voo, vof, voc, vogefi, vot

    end

end

####################################################################

function call_germline_variant(mo, ta, ge, fa, chs, chn, pao, n_jo, me, to, sn, rs, va)

    Fastq.support.log()

    Fastq.support.index_genome_files(fa, chs)

    Fastq.support.error_if_directory(pao)


    # Run docker container

    id, voo, vof, voc, vogefi, vot =
        run_strelka_manta_docker_container(to, fa, chs, ge, pao, nothing)


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

    Fastq.support.remove_docker_container(id)


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

    papa = Fastq.vcf.annotate_with_snpeff(pao, me, sn, paco, n_jo)

    # snpsift 

    if rs

        Fastq.vcf.annotate_with_snpsift(pao, sn, va, papa, n_jo)

    end

end

####################################################################

function call_somatic_variant(ta, ge, so, fa, chs, chn, pao, n_jo, me, to, sn, rs, va)
    Fastq.support.log()

    Fastq.support.index_genome_files(fa, chs)

    Fastq.support.error_if_directory(pao)


    # Run docker container

    id, voo, vof, voc, vogefi, vosofi, vot =
        run_strelka_manta_docker_container(to, fa, chs, ge, pao, so)


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

    Fastq.support.remove_docker_container(id)


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

    papa = Fastq.vcf.annotate_with_snpeff(pao, me, sn, paco, n_jo)


    # snpsift

    if rs

        Fastq.vcf.annotate_with_snpsift(pao, sn, va, papa, n_jo)

    end

end

end