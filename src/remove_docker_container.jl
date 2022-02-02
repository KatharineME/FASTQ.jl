function remove_docker_container(id::Vector{String})::Nothing

    run(`docker kill $id`)
    
    run(`docker rm $id`)

    println("Docker container was removed")

    return nothing

end
