function test_strelka_and_manta(pa)

    for pr in [MANTA, STRELKA]

        if !(pr in readdir(pa))

            error("You dont have the correct version ($pr)")

        end

    end

    vo = last(split(pa, "/"))

    id = readlines(
        pipeline(
            `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
        ),
    )

    for sc in [
        joinpath(MANTA, "bin", "runMantaWorkflowDemo.py"),
        joinpath(STRELKA, "bin", "runStrelkaGermlineWorkflowDemo.bash"),
        joinpath(STRELKA, "bin", "runStrelkaSomaticWorkflowDemo.bash"),
    ]

        re = readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`))

        println("$(join(re, " "))\n")

    end

    remove_docker_container(id)

    return

end
