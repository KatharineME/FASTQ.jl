function run_docker_container(to, fa, chs, ge, pao, so)

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
