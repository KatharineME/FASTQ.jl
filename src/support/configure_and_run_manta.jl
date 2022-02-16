function configure_and_run_manta(voo, id, vot, co, ru)

    vom = joinpath(voo, "manta")

    vomr = joinpath(vom, "runWorkflow.py")

    sc = "$MANTA/bin/configManta.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$vom && ./home/$vomr $ru"`,
        ),
    )

    println("$(join(re, " "))\n")

    return vom

end
