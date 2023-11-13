module BAM

using ..FASTQ

const HO = "home"

function _run_strelka_manta_docker_container(pa, bage, re, chs, to; baso = nothing)

    bage, re, pa = [FASTQ.Support.make_path_absolute(pa) for pa in (bage, re, pa)]

    vot = basename(to)

    page = dirname(bage)

    voge = joinpath(basename(dirname(page)), basename(page))

    par = dirname(re)

    vor = basename(par)

    vorfi = joinpath(vor, basename(re))

    voo = basename(pa)

    vo = "volume"

    vo_ = (
        "--$vo",
        "$to:/$HO/$vot",
        "--$vo",
        "$page:/$HO/$voge",
        "--$vo",
        "$par:/$HO/$vor",
        "--$vo",
        "$pa:/$HO/$voo",
    )

    com = ("docker", "run", "--interactive", "--detach", "--tty", "--user", "root", "--memory=30g")

    con = ("centos:centos6", "bash")

    vost = joinpath(voo, "strelka")

    vostr = joinpath(vost, "runWorkflow.py")

    voc = joinpath(vor, basename(chs))

    vogefi = joinpath(voge, basename(bage))

    if baso != nothing

        paso = dirname(FASTQ.Support.make_path_absolute(baso))

        voso = basename(paso)

        vosofi = joinpath(voso, basename(baso))

        id = readlines(pipeline(`$com $vo_ --$vo $paso:/$HO/$voso $con`))

        return id, voo, vost, vostr, vorfi, voc, vogefi, vosofi, vot

    else

        id = readlines(pipeline(`$com $vo_ $con`))

        return id, voo, vost, vostr, vorfi, voc, vogefi, vot

    end

end

function _set_output_path(pa)

    past, pama = [joinpath(pa, st) for st in ("strelka", "manta")]

    pav = joinpath("results", "variants")

    past, pama, pav

end

function _set_strelka_manta_run(n_jo, me)

    "--mode local --jobs $n_jo --memGb $me --quiet"

end

function _configure_and_run_manta(id, co, ru, voo, vot)

    vom = joinpath(voo, "manta")

    vomr = joinpath(vom, "runWorkflow.py")

    sc = "$(FASTQ._MA)/bin/configManta.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./$HO/$vot/$(sc) $co --outputContig --runDir /$HO/$vom && ./$HO/$vomr $ru"`,
        ),
    )

    @info join(re, " ")

    vom

end

function call_germline_variant(pa, ba, mo, ex, re, chs, chn, va, to, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    id, voo, vost, vostr, vorfi, voc, voba, vot =
        _run_strelka_manta_docker_container(pa, ba, re, chs, to)

    _, _, pav = _set_output_path(pa)

    co = "--referenceFasta /$HO/$vorfi --callRegions $HO/$voc --bam $HO/$voba"

    vc_ = (
        joinpath(pa, "manta", pav, "diploidSV.vcf.gz"),
        joinpath(pa, "strelka", pav, "variants.vcf.gz"),
    )

    if ex

        co = "$co --exome"

    end

    if mo == "cdna"

        co = "$co --rna"

        vc_ = [joinpath(pa, "strelka", pav, "variants.vcf.gz")]

    end

    ru = _set_strelka_manta_run(n_jo, me)

    _configure_and_run_manta(id, co, ru, voo, vot)

    sc = joinpath(FASTQ._ST, "bin", "configureStrelkaGermlineWorkflow.py")

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./$HO/$vot/$(sc) $co --runDir /$HO/$vost && ./$HO/$vostr $ru"`,
        ),
    )

    FASTQ.Support.remove_docker_container(id)

    vc_

end

function call_somatic_variant(pa, bage, baso, ex, re, chs, chn, va, to, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    id, voo, vost, vostr, vorfi, voc, vogefi, vosofi, vot =
        _run_strelka_manta_docker_container(pa, bage, re, chs, to; baso = baso)

    past, pama, pav = _set_output_path(pa)

    co = "--referenceFasta /$HO/$vorfi --callRegions /$HO/$voc --normalBam /$HO/$vogefi --tumorBam /$HO/$vosofi"

    if ex

        co = "$co --exome"

    end

    ru = _set_strelka_manta_run(n_jo, me)

    vom = _configure_and_run_manta(id, co, ru, voo, vot)

    sc = joinpath(FASTQ._ST, "bin", "configureStrelkaSomaticWorkflow.py")

    ca = joinpath("$HO", vom, pav, "candidateSmallIndels.vcf.gz")

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./$HO/$vot/$(sc) $co --indelCandidates $ca --runDir /$HO/$vost && ./$HO/$vostr $ru"`,
        ),
    )

    FASTQ.Support.remove_docker_container(id)

    sa = joinpath(pa, "sample.txt")

    open(io -> write(io, "Germline\nSomatic"), sa; write = true)

    vc_ = (
        FASTQ.VCF.reheader_vcf(joinpath(past, pav, "somatic.indels.vcf.gz"), sa, n_jo),
        FASTQ.VCF.reheader_vcf(joinpath(past, pav, "somatic.snvs.vcf.gz"), sa, n_jo),
        FASTQ.VCF.reheader_vcf(joinpath(pama, pav, "somaticSV.vcf.gz"), sa, n_jo),
    )

    vc_

end

end
