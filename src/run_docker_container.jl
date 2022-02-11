function run_docker_container(to, fa, chs, ge, pao, so)

    vot = basename(to)

    page = dirname(abspath(ge))

    voge = basename(page)

    vogefi = joinpath(voge, basename(ge))

    pag = dirname(abspath(fa))

    vog = basename(pag)

    vof = joinpath(vog, basename(fa))

    voc = joinpath(vog, "chromosome", basename(chs))

    pao = abspath(pao)

    voo = basename(pao)

    if so !== nothing

        paso = dirname(abspath(so))

        voso = basename(paso)

        vosofi = joinpath(voso, basename(so))

        id = readlines(
            pipeline(
                `docker run --interactive --detach --tty --user root --memory=30g --volume $to:/home/$vot --volume $page:/home/$voge --volume $pag:/home/$vog --volume $paso:/home/$voso --volume $pao:/home/$voo centos:centos6 bash`,
            ),
        )

        return id, voo, vof, voc, vogefi, vosofi, vot

    else

        id = readlines(
            pipeline(
                `docker run --interactive --detach --tty --user root --memory=30g --volume $to:/home/$vot --volume $page:/home/$voge --volume $pag:/home/$vog --volume $pao:/home/$voo centos:centos6 bash`,
            ),
        )

        return id, voo, vof, voc, vogefi, vot

    end

end
