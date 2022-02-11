function index_genome_files(fa, chs)

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    return

end
