module Reference

using ..FASTQ

function get_chromosome_file_path(re)

    gr = dirname(re)

    joinpath(gr, "chromosome.bed.gz"), joinpath(gr, "chrn_n.tsv")

end

function index_genome_file(ge, chs)

    geu = ge

    if endswith(ge, ".gz")

        un = splitext(ge)[1]

        isfile(un) ? geu = un : nothing

    end

    for fi in (ge, geu)

        !(isfile("$fi.fai") && ispath("$fi.gzi")) ? run(`samtools faidx $fi`) : nothing

    end

    !ispath("$chs.tbi") ? run(`tabix --force $chs`) : nothing

    nothing

end

function generate_star_genome_file(ge, n_jo; ga = nothing)

    id = joinpath(dirname(ge), "StarIndex")

    if ga === nothing

        ga = joinpath(
            FASTQ._DA,
            "GRCh38",
            "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set",
            "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf",
        )

    end

    if !ispath(joinpath(id, "geneInfo.tab"))

        ispath(id) ? rm(id, recursive = true) : nothing

        mkdir(id)

        ged = splitext(ge)[1]

        !isfile(ged) ? run(pipeline(`bgzip --decompress --stdout $ge`, ged)) : nothing

        run(
            `star --runMode genomeGenerate --runThreadN $n_jo --genomeDir $id --genomeFastaFiles $ged --sjdbGTFfile $ga`,
        )

    end

    id

end

end
