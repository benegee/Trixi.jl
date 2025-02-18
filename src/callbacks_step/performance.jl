# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    PerformanceDataCallback(interval)

Gather some performance metrics every `interval` steps.
"""
mutable struct PerformanceDataCallback
    interval::Int
end

function PerformanceDataCallback(; interval = 0)
    performance_callback = PerformanceDataCallback(interval)

    DiscreteCallback(performance_callback, # the first one is the condition,
                     performance_callback, # the second the affect!
                     save_positions = (false, false),
                     initialize = initialize!)
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:PerformanceDataCallback})
    @nospecialize cb # reduce precompilation time

    performance_callback = cb.affect!
    @unpack interval = performance_callback
    print(io, "PerformanceDataCallback(interval=", interval, ")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:PerformanceDataCallback})
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        performance_callback = cb.affect!

        setup = [
            "Interval" => performance_callback.interval
        ]
        summary_box(io, "PerformanceDataCallback", setup)
    end
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u, t,
                     integrator) where {Condition, Affect! <: PerformanceDataCallback}
    cb.affect!(integrator)
end

# This method is called to determine whether the callback should be activated
function (performance_callback::PerformanceDataCallback)(u, t, integrator)
    @unpack interval = performance_callback

    # With error-based step size control, some steps can be rejected. Thus,
    #   `integrator.iter >= integrator.stats.naccept`
    #    (total #steps)       (#accepted steps)
    # We need to check the number of accepted steps since callbacks are not
    # activated after a rejected step.
    return interval > 0 && (integrator.stats.naccept % interval == 0 ||
                            isfinished(integrator))
end

# This method is called when the callback is activated
function (performance_callback::PerformanceDataCallback)(integrator)
    @trixi_timeit timer() "performance data" begin
        # General information
        mpi_println()
        mpi_println("─"^100)
        mpi_println(" Performance data")
        for rank in 0:mpi_nranks()-1
            if mpi_rank() == rank
                gpu = CUDA.device()
                println("Rank $rank, device $(gpu), ID $(CUDA.uuid(gpu))")
                CUDA.pool_status()
            end
            MPI.Barrier(mpi_comm())
        end
        mpi_println("─"^100)
    end
    return nothing
end
end # @muladd
