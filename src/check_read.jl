function check_read(fq_::Array, di::String, n_jo::Int)::Nothing

    if ispath(di)

        println(
            "Skipping because directory already exists:\n $di\n",
        )

    else

        mkpath(di)

        println("Running FastQC ...")

        run(
            `fastqc --threads $(minimum((length(fq_), n_jo))) --outdir $di $fq_`,
           )

        println("Running MultiQC ...")

        run(`multiqc --outdir $di $di`)

    end

    return nothing 

end
