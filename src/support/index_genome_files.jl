function index_genome_files(ge, chs)

    if !(isfile("$ge.fai") && ispath("$ge.gzi"))

        run(`samtools faidx $ge`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    return

end
