function test_strelka_and_manta(pa::String)::Nothing

    ma = "manta-1.6.0.centos6_x86_64"

    st = "strelka-2.9.10.centos6_x86_64"

    for pr in [ma, st]

        if !(pr in readdir(pa))

            println("You dont have the correct version ($pr)")

            return

        end

    end
    
    vo = last(split(pa, "/"))

    id = readlines(pipeline(
                 `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
               )
       )

    for sc in [
        joinpath(ma, "bin", "runMantaWorkflowDemo.py"),
        joinpath(st, "bin", "runStrelkaGermlineWorkflowDemo.bash"),
        joinpath(st, "bin", "runStrelkaSomaticWorkflowDemo.bash"),
    ]

        re =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`))

        println("$(join(re, " "))\n")

    end

    run(`docker kill $id`)
    
    run(`docker rm $id`)

    return nothing

end
