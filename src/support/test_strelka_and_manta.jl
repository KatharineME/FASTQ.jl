function test_strelka_and_manta(pa)

    for pr in (Fastq.MANTA, Fastq.STRELKA)

        if !(pr in readdir(pa))

            error("you dont have the correct version ($pr)")

        end

    end

    vo = last(split(pa, "/"))

    id = readlines(
        pipeline(
            `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
        ),
    )

    for sc in [
        joinpath(Fastq.MANTA, "bin", "runMantaWorkflowDemo.py"),
        joinpath(Fastq.STRELKA, "bin", "runStrelkaGermlineWorkflowDemo.bash"),
        joinpath(Fastq.STRELKA, "bin", "runStrelkaSomaticWorkflowDemo.bash"),
    ]

        re = readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`))

        println("$(join(re, " "))\n")

    end

    Fastq.bam.remove_docker_container(id)

    return

end
