module Support

using Dates

using BioLab

using FASTQ

function error_if_directory(pa)

    paf = BioLab.Path.make_absolute(pa)

    if ispath(paf)

        error(
            "\nSkipping $(replace(basename(paf), "_" => " ")) because directory already exists:\n $pa\n",
        )

    else

        mkpath(paf)

    end

end

function index_genome_files(ge, chs)

    if !(isfile("$ge.fai") && ispath("$ge.gzi"))

        run(`samtools faidx $ge`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

end

function log()

    println("\n", "="^99)

    println("$(StackTraces.stacktrace()[2].func)\n")

    println(BioLab.Time.stamp())

    println("="^99, "\n")

end

function remove_docker_container(id)

    run(`docker kill $id`)

    run(`docker rm $id`)

    println("\nDocker container was removed\n")

end

function test_local_environment()

    println("Checking for programs\n")

    for pr in [
        "fastp",
        "fastqc",
        "multiqc",
        "bgzip",
        "tabix",
        "minimap2",
        "samtools",
        "bcftools",
        "kallisto",
        "star",
        "docker",
    ]

        run(`which $pr`)

    end

end

function test_strelka_and_manta(pa)

    for pr in (FASTQ.MANTA, FASTQ.STRELKA)

        if !(pr in readdir(pa))

            error("you dont have the correct version ($pr)")

        end

    end

    vo = last(split(pa, "/"))

    id = readlines(
        pipeline(
            `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
        ),
    )

    for sc in [
        joinpath(FASTQ.MANTA, "bin", "runMantaWorkflowDemo.py"),
        joinpath(FASTQ.STRELKA, "bin", "runStrelkaGermlineWorkflowDemo.bash"),
        joinpath(FASTQ.STRELKA, "bin", "runStrelkaSomaticWorkflowDemo.bash"),
    ]

        re = readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`))

        println("$(join(re, " "))\n")

    end

    FASTQ.BAM.remove_docker_container(id)

end

end
