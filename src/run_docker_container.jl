function run_docker_container(pa::String, vo::SubString{String})

    vo = last(split(pa, "/"))

    id = readlines(pipeline(
                 `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
               )
       )

    return id

end
