function remove_docker_container(id::Vector{String})::Nothing

    run(`docker kill $id`)

    run(`docker rm $id`)

    println("\nDocker container was removed\n")

    return nothing

end
