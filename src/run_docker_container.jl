function run_docker_container(to::String, fa::String, chs::String, ba::String, pao::String)

    vot = basename(to)

    pab = dirname(abspath(ge))

    vob = basename(pab)

    pag = dirname(abspath(fa))

    vog = basename(pag)
    
    vof = joinpath(vog, basename(fa))

    voc = joinpath(vog, "chromosome", basename(chs))

    voge = joinpath(vob, basename(ge))

    pao = abspath(pao)

    voo = basename(pao)
    
    id = readlines(pipeline(
                 `docker run --interactive --detach --tty --user root --memory=30g --volume $to:/home/$vot --volume $pab:/home/$vob --volume $pag:/home/$vog --volume $pao:/home/$voo centos:centos6 bash`,
               )
       )

    return id, voo, vof, voc, voge, vot

end

