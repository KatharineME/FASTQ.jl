function check_read(fq_::Array, di::String, n_jo::Int)::Nothing

    if check_directory(di, "check read")

        return nothing

    end

    println("Running FastQC ...")

    run(`fastqc --threads $(minimum((length(fq_), n_jo))) --outdir $di $fq_`)

    println("Running MultiQC ...")

    run(`multiqc --outdir $di $di`)

    return nothing

end
