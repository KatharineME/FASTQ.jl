function remove_docker_container(id)

    run(`docker kill $id`)

    run(`docker rm $id`)

    println("\nDocker container was removed\n")

    return

end
