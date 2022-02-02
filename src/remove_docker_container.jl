function remove_docker_container(id::Vector{String})::Nothing

    run(`docker kill $id`)
    
    run(`docker rm $id`)

    return nothing

end
