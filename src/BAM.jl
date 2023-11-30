module BAM

using ..FASTQ

const _HO = "home"

function _run_strelka_manta_docker_container(pa, bage, re, chs, to; baso = nothing)

    bage, re, pa = [FASTQ.Support.make_path_absolute(pan) for pan in (bage, re, pa)]

    vot, voo = [basename(pan) for pan in (to, pa)]

    page = dirname(bage)

    voge = joinpath(basename(dirname(page)), basename(page))

    par = dirname(re)

    vor = basename(par)

    vorfi = joinpath(vor, basename(re))

    vo = "volume"

    vo_ = (
        "--$vo",
        "$to:/$_HO/$vot",
        "--$vo",
        "$page:/$_HO/$voge",
        "--$vo",
        "$par:/$_HO/$vor",
        "--$vo",
        "$pa:/$_HO/$voo",
    )

    com = ("docker", "run", "--interactive", "--detach", "--tty", "--user", "root", "--memory=30g")

    con = ("centos:centos6", "bash")

    vost = joinpath(voo, "strelka")

    vostr = joinpath(vost, "runWorkflow.py")

    voc = joinpath(vor, basename(chs))

    vogefi = joinpath(voge, basename(bage))

    if baso === nothing

        id = readchomp(`$com $vo_ $con`)

        return id, voo, vost, vostr, vorfi, voc, vogefi, vot

    else

        paso = dirname(FASTQ.Support.make_path_absolute(baso))

        voso = basename(paso)

        vosofi = joinpath(voso, basename(baso))

        id = readchomp(`$com $vo_ --$vo $paso:/$_HO/$voso $con`)

        return id, voo, vost, vostr, vorfi, voc, vogefi, vosofi, vot

    end

end

function _set_output_path(pa)

    [joinpath(pa, st) for st in ("strelka", "manta")]..., joinpath("results", "variants")

end

function _set_strelka_manta_run(n_jo, me)

    "--mode local --jobs $n_jo --memGb $me --quiet"

end

function _configure_and_run_manta(id, co, ru, voo, vot)

    vom = joinpath(voo, "manta")

    run(
        `docker exec --interactive $id bash -c "./$_HO/$vot/$(FASTQ._MA)/bin/configManta.py $co --outputContig --runDir /$_HO/$vom && ./$_HO/$(joinpath(vom, "runWorkflow.py")) $ru"`,
    )

    vom

end

function call_germline_variant(pa, ba, mo, ex, re, chs, to, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    id, voo, vost, vostr, vorfi, voc, voba, vot =
        _run_strelka_manta_docker_container(pa, ba, re, chs, to)

    past, pama, pav = _set_output_path(pa)

    co = "--referenceFasta /$_HO/$vorfi --callRegions $_HO/$voc --bam $_HO/$voba"

    vc_ = (joinpath(pama, pav, "diploidSV.vcf.gz"), joinpath(past, pav, "variants.vcf.gz"))

    ex ? co = "$co --exome" : nothing

    if mo == "cdna"

        co = "$co --rna"

        vc_ = [joinpath(pa, "strelka", pav, "variants.vcf.gz")]

    end

    ru = _set_strelka_manta_run(n_jo, me)

    _configure_and_run_manta(id, co, ru, voo, vot)

    run(
        `docker exec --interactive $id bash -c "./$_HO/$vot/$(joinpath(FASTQ._ST, "bin", "configureStrelkaGermlineWorkflow.py")) $co --runDir /$_HO/$vost && ./$_HO/$vostr $ru"`,
    )

    FASTQ.Support.remove_docker_container(id)

    vc_

end

function call_somatic_variant(pa, bage, baso, ex, re, chs, to, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    id, voo, vost, vostr, vorfi, voc, vogefi, vosofi, vot =
        _run_strelka_manta_docker_container(pa, bage, re, chs, to; baso = baso)

    past, pama, pav = _set_output_path(pa)

    co = "--referenceFasta /$_HO/$vorfi --callRegions /$_HO/$voc --normalBam /$_HO/$vogefi --tumorBam /$_HO/$vosofi"

    ex ? co = "$co --exome" : nothing

    ru = _set_strelka_manta_run(n_jo, me)

    vom = _configure_and_run_manta(id, co, ru, voo, vot)

    run(
        `docker exec --interactive $id bash -c "./$_HO/$vot/$(joinpath(FASTQ._ST, "bin", "configureStrelkaSomaticWorkflow.py")) $co --indelCandidates $(joinpath(_HO, vom, pav, "candidateSmallIndels.vcf.gz")) --runDir /$_HO/$vost && ./$_HO/$vostr $ru"`,
    )

    FASTQ.Support.remove_docker_container(id)

    sa = joinpath(pa, "sample.txt")

    open(io -> write(io, "Germline\nSomatic"), sa; write = true)

    [
        FASTQ.VCF.reheader_vcf(joinpath(el...), sa, n_jo) for el in (
            (past, pav, "somatic.indels.vcf.gz"),
            (past, pav, "somatic.snvs.vcf.gz"),
            (pama, pav, "somaticSV.vcf.gz"),
        )
    ]

end

end
