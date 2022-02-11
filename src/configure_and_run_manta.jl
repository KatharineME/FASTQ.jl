function configure_and_run_manta(voo, id, vot, co, ru)

    pam = joinpath(voo, "manta")

    pamr = joinpath(pam, "runWorkflow.py")

    sc = "$MANTA/bin/configManta.py"

    re = readlines(
        pipeline(
            `docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$pam && ./home/$pamr $ru"`,
        ),
    )

    println("$(join(re, " "))\n")

    return pam

end
