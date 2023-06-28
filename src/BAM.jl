module BAM

using FASTQ

function _prepare_for_variant_calling(fa, chs, pa)

    FASTQ.Support.index_genome_files(fa, chs)

    FASTQ.Support.error_if_directory(pa)

end


function _configure_and_run_manta(voo, id, vot, co, ru)

    vom = joinpath(voo, "manta")

    vomr = joinpath(vom, "runWorkflow.py")

    sc = "$(FASTQ.MANTA)/bin/configManta.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$vom && ./home/$vomr $ru"`,
        ),
    )

    @info join(re, " ")

    vom

end

function _set_strelka_manta_run(n_jo, me)

    return "--mode local --jobs $n_jo --memGb $me --quiet"

end

function _run_strelka_manta_docker_container(to, ge, fa, chs, pao; so = nothing)

    com = ["docker", "run", "--interactive", "--detach", "--tty", "--user", "root", "--memory=30g"]

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

    vo_ = [
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

        id = readlines(pipeline(`$com $vo_ --volume $paso:/home/$voso $con`))

        return id, voo, vof, voc, vogefi, vosofi, vot

    else

        id = readlines(pipeline(`$com $vo_ $con`))

        return id, voo, vof, voc, vogefi, vot

    end

end

function call_germline_variant(pao, to, ge, fa, chs, ta, mo, n_jo, me, chn, sn, rs, va)

    FASTQ.Support.log()

    _prepare_for_variant_calling(pao)

    id, voo, vof, voc, vogefi, vot = _run_strelka_manta_docker_container(to, ge, fa, chs, pao)

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

    @info join(re, " ")

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

    FASTQ.VCF.combine_vcf(n_jo, vc_, chn, paco)

    run(`tabix $paco`)

    papa = FASTQ.VCF.annotate_with_snpeff(pao, me, sn, paco, n_jo)

    if rs

        FASTQ.VCF.annotate_with_snpsift(pao, sn, va, papa, n_jo)

    end

end

function call_somatic_variant(pao, to, ge, fa, chs, so, ta, n_jo, me, chn, sn, rs, va)

    FASTQ.Support.log()

    _prepare_for_variant_calling(pao)

    id, voo, vof, voc, vogefi, vosofi, vot =
        _run_strelka_manta_docker_container(to, ge, fa, chs, pao, so = so)

    co = "--referenceFasta /home/$vof --callRegions /home/$voc --normalBam /home/$vogefi --tumorBam /home/$vosofi"

    if ta

        co = "$co --exome"

    end

    ru = _set_strelka_manta_run(n_jo, me)

    vom = _configure_and_run_manta(voo, id, vot, co, ru)

    pav = joinpath("results", "variants")

    vost = joinpath(voo, "strelka")

    sc = "$(FASTQ.STRELKA)/bin/configureStrelkaSomaticWorkflow.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --indelCandidates $(joinpath("home", vom, pav, "candidateSmallIndels.vcf.gz")) --runDir /home/$vost && ./home/$(joinpath(vost, "runWorkflow.py")) $ru"`,
        ),
    )

    @info join(re, " ")

    FASTQ.Support.remove_docker_container(id)

    past = joinpath(pao, "strelka")

    sa = joinpath(pao, "sample.txt")

    open(io -> write(io, "Germline\nSomatic"), sa; write = true)

    vc_ = [
        FASTQ.VCF.reheader_vcf(joinpath(past, pav, "somatic.indels.vcf.gz"), n_jo, sa),
        FASTQ.VCF.reheader_vcf(joinpath(past, pav, "somatic.snvs.vcf.gz"), n_jo, sa),
        FASTQ.VCF.reheader_vcf(joinpath(pao, "manta", pav, "somaticSV.vcf.gz"), n_jo, sa),
    ]

    paco = joinpath(pao, "concat.vcf.gz")

    FASTQ.VCF.combine_vcf(n_jo, vc_, chn, paco)

    run(`tabix $paco`)

    papa = FASTQ.VCF.annotate_with_snpeff(pao, me, sn, paco, n_jo)

    if rs

        FASTQ.VCF.annotate_with_snpsift(pao, sn, va, papa, n_jo)

    end

end

end
