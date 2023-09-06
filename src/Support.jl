module Support

using Dates

using BioLab

using FASTQ


function make_path_absolute(pa)

    rstrip(abspath(expanduser(pa)), '/')

end



function error_if_file_missing(fi)

    if !isfile(fi)

        error("File does not exist: $fi")

    end

end

function trash_remake_directory(di)

    dia = rstrip(abspath(expanduser(di)), '/')

    if ispath(dia)

        @warn "Trashing $dia"

        mv(dia, joinpath(joinpath(homedir(), ".Trash"), basename(dia)), force = true)

    elseif ispath(di)

        error("$di is not a directory.")

    end

    @info "Making $dia"

    mkpath(dia)

end

function index_genome_files(ge, chs)

    if !(isfile("$ge.fai") && ispath("$ge.gzi"))

        run(`samtools faidx $ge`)

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

end

function log_top_level_function()

    fu = replace(string(StackTraces.stacktrace()[2].func), "_" => " ")

    ti = BioLab.Time.stamp()

    @info "==============================
    $fu 
    $ti
    ================================"
end

function log_sub_level_function()

    fu = replace(string(StackTraces.stacktrace()[2].func), "_" => " ")

    ti = BioLab.Time.stamp()

    @info "------------------------------
    $fu 
    $ti
    --------------------------------"

end

function remove_docker_container(id)

    run(`docker kill $id`)

    run(`docker rm $id`)

    @info "Docker container was removed."

end

function test_local_environment()

    @info "Checking for programs"

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

    for pr in (FASTQ._MA, FASTQ._ST)

        if !(pr in readdir(pa))

            error("You dont have the correct version ($pr).")

        end

    end

    vo = last(split(pa, "/"))

    id = readlines(
        pipeline(
            `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
        ),
    )

    for sc in [
        joinpath(FASTQ._MA, "bin", "runMantaWorkflowDemo.py"),
        joinpath(FASTQ._ST, "bin", "runStrelkaGermlineWorkflowDemo.bash"),
        joinpath(FASTQ._ST, "bin", "runStrelkaSomaticWorkflowDemo.bash"),
    ]

        re = readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`))

        @info join(re, " ")

    end

    remove_docker_container(id)

end

end
