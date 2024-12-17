using Trixi
using MPI
using TimerOutputs
using CUDA

function main(elixir_path)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    isroot = rank == 0

    # pin rank to device?
    #if machine == "jedi"
    #    CUDA.device!(rank % 4)
    #end
    print("Rank $rank has device: $(CUDA.device())\n")

    # setup
    maxiters = 400

    if isroot
        println("Warming up...")
    end

    # start simulation with tiny final time to trigger precompilation
    duration_precompile = @elapsed trixi_include(elixir_path,
        tspan=(0.0, 1e-14))

    if isroot
        println("Finished warm-up in $duration_precompile seconds\n")
        println("Starting simulation...")
    end

    # start the real simulation
    duration_elixir = @elapsed trixi_include(elixir_path, maxiters=maxiters)

    # store metrics (on every rank!)
    metrics = Dict{String, Float64}("elapsed time" => duration_elixir)

    # read TimerOutputs timings
    timer = Trixi.timer()
    metrics["total time"] = 1.0e-9 * TimerOutputs.tottime(timer)
    metrics["rhs! time"] = 1.0e-9 * TimerOutputs.time(timer["rhs!"])

    # compute performance index
    nrhscalls = Trixi.ncalls(semi.performance_counter)
    walltime = 1.0e-9 * take!(semi.performance_counter)
    metrics["PID"] = walltime * Trixi.mpi_nranks() / (Trixi.ndofsglobal(semi) * nrhscalls)

    # gather metrics from all ranks
    gathered_metrics = MPI.gather(metrics, comm)

    if isroot
        # reduce metrics per rank
        open("metrics.out", "w") do io
            for (key, _) in gathered_metrics[1]
                println(io, key, ": ", mapreduce(x->x[key], min, gathered_metrics))
            end
        end
    end
end

# hardcoded elixir
elixir_path = joinpath(@__DIR__(), "elixirs/elixir_euler_taylor_green_vortex.jl")

main(elixir_path)
