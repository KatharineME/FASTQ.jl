function index_genome_files(fa::String, chs::String)::Nothing

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end
    
    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    return nothing

end

