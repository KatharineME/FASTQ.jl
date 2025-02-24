module Support

using Dates: now

using ..FASTQ

function check_os()

    readchomp(`uname`)

end

function start_docker()

    try

        run(`docker ps -a`)

    catch

        run(`open -a docker`)

        @info "Started Docker."

    else

        @info "Docker is already running."

    end

end

function test_local_environment()

    @info "Checking for programs"

    for pr in (
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
    )

        run(`which $pr`)

    end

end

function log_top_level_function()

    fu = _clean_function_name(string(StackTraces.stacktrace()[2].func))

    @info "============================================================
    $fu 
    $(now())
    =============================================================="

end

function log_sub_level_function()

    fu = _clean_function_name(string(StackTraces.stacktrace()[2].func))

    @info "------------------------------
    $fu 
    $(now())
    --------------------------------"

end

function calculate_size(fi)

    parse(Float64, split(Base.format_bytes(stat(fi).size), " ")[1])

end

function _clean_function_name(na)

    sp = split(na, "#")

    id = findmax([lastindex(st) for st in sp])[2]

    replace(sp[id], "_" => " ")

end

function error_if_file_missing(fi_)

    for fi in fi_

        !isfile(fi) ? error("File does not exist: $fi") : nothing

    end

    nothing

end

function make_path_absolute(pa)

    rstrip(abspath(expanduser(pa)), '/')

end

function make_sample_to_fastq_dictionary(di, na)

    dia = FASTQ.Support.make_path_absolute(di)

    sa_fq_ = Dict{String, Tuple{String, String}}()

    for sa in readdir(dia; join = true)

        if isdir(sa)

            fq_ = FASTQ.Raw.find(sa)

            le = lastindex(fq_)

            if le > 2

                error("Sample $sa has more than two fastq files. Concatenate first.")

            elseif 2 > le

                error("Sample $sa has less than two fastq files.")

            else

                fq1 = fq_[1]

                if occursin(na, fq1)

                    fq2 = replace(fq1, na => replace(na, "1" => "2"))

                else

                    fq2, fq1 = fq1, replace(fq1, na => replace(na, "2" => "1"))

                end

                sa_fq_[sa] = (fq1, fq2)

            end

        end

    end

    sa_fq_

end

function trash_remake_directory(di)

    dia = make_path_absolute(di)

    if ispath(dia)

        @warn "Trashing $dia"

        mv(dia, joinpath(joinpath(homedir(), ".Trash"), basename(dia)); force = true)

    end

    @info "Making $dia"

    mkdir(dia)

end

function make_analysis_directory(di, to, su_; sa_fq_ = nothing)

    toa = trash_remake_directory(joinpath(di, to))

    sua_ = Vector{String}()

    for su in su_

        sua = mkpath(joinpath(toa, su))

        push!(sua_, sua)

        if sa_fq_ !== nothing

            for sa in keys(sa_fq_)

                mkdir(joinpath(sua, basename(sa)))

            end

        end

    end

    sua_

end

function remove_docker_container(id)

    run(`docker kill $id`)

    run(`docker rm $id`)

    @info "Docker container was removed."

    nothing

end

function test_strelka_and_manta(pa)

    for pr in (FASTQ._MA, FASTQ._ST)

        !(pr in readdir(pa)) ? error("You dont have the correct version ($pr).") : nothing

    end

    vo = basename(pa)

    id = readchomp(
        `docker run --interactive --detach --tty --user root --volume $pa:/home/$vo centos:centos6 bash`,
    )

    for sc in (
        joinpath(FASTQ._MA, "bin", "runMantaWorkflowDemo.py"),
        joinpath(FASTQ._ST, "bin", "runStrelkaGermlineWorkflowDemo.bash"),
        joinpath(FASTQ._ST, "bin", "runStrelkaSomaticWorkflowDemo.bash"),
    )

        run(`docker exec --interactive $id bash -c "./home/$vo/$(sc)"`)

    end

    remove_docker_container(id)

    nothing

end

end
