# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

mutable struct VisualizationCallback{SolutionVariables, VariableNames, PlotDataCreator,
                                     PlotCreator}
    interval::Int
    solution_variables::SolutionVariables
    variable_names::VariableNames
    show_mesh::Bool
    plot_data_creator::PlotDataCreator
    plot_creator::PlotCreator
    plot_arguments::Dict{Symbol, Any}
end

function Base.show(io::IO,
                   cb::DiscreteCallback{Condition, Affect!}) where {Condition,
                                                                    Affect! <:
                                                                    VisualizationCallback
                                                                    }
    visualization_callback = cb.affect!
    @unpack interval, plot_arguments, solution_variables, variable_names, show_mesh, plot_creator, plot_data_creator = visualization_callback
    print(io, "VisualizationCallback(",
          "interval=", interval, ", ",
          "solution_variables=", solution_variables, ", ",
          "variable_names=", variable_names, ", ",
          "show_mesh=", show_mesh, ", ",
          "plot_data_creator=", plot_data_creator, ", ",
          "plot_creator=", plot_creator, ", ",
          "plot_arguments=", plot_arguments, ")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{Condition, Affect!}) where {Condition,
                                                                    Affect! <:
                                                                    VisualizationCallback
                                                                    }
    if get(io, :compact, false)
        show(io, cb)
    else
        visualization_callback = cb.affect!

        setup = [
            "interval" => visualization_callback.interval,
            "plot arguments" => visualization_callback.plot_arguments,
            "solution variables" => visualization_callback.solution_variables,
            "variable names" => visualization_callback.variable_names,
            "show mesh" => visualization_callback.show_mesh,
            "plot creator" => visualization_callback.plot_creator,
            "plot data creator" => visualization_callback.plot_data_creator
        ]
        summary_box(io, "VisualizationCallback", setup)
    end
end

"""
    VisualizationCallback(; interval=0,
                            solution_variables=cons2prim,
                            variable_names=[],
                            show_mesh=false,
                            plot_data_creator=PlotData2D,
                            plot_creator=show_plot,
                            plot_arguments...)

Create a callback that visualizes results during a simulation, also known as *in-situ
visualization*.

!!! warning "Experimental implementation"
    This is an experimental feature and may change in any future releases.

The `interval` specifies the number of time step iterations after which a new plot is generated. The
available variables to plot are configured with the `solution_variables` parameter, which acts the
same way as for the [`SaveSolutionCallback`](@ref). The variables to be actually plotted can be
selected by providing a single string or a list of strings to `variable_names`, and if `show_mesh`
is `true`, an additional plot with the mesh will be generated.

To customize the generated figure, `plot_data_creator` allows to use different plot data types. With
`plot_creator` you can further specify an own function to visualize results, which must support the
same interface as the default implementation [`show_plot`](@ref). All remaining
keyword arguments are collected and passed as additional arguments to the plotting command.
"""
function VisualizationCallback(; interval = 0,
                               solution_variables = cons2prim,
                               variable_names = [],
                               show_mesh = false,
                               plot_data_creator = PlotData2D,
                               plot_creator = show_plot,
                               plot_arguments...)
    mpi_isparallel() && error("this callback does not work in parallel yet")

    if variable_names isa String
        variable_names = String[variable_names]
    end

    visualization_callback = VisualizationCallback(interval,
                                                   solution_variables, variable_names,
                                                   show_mesh,
                                                   plot_data_creator, plot_creator,
                                                   Dict{Symbol, Any}(plot_arguments))

    # Warn users if they create a visualization callback without having loaded the Plots package
    #
    # Note: This warning is added for convenience, as Plots is the only "officially" supported
    #       visualization package right now. However, in general nothing prevents anyone from using
    #       other packages such as Makie, Gadfly etc., given that appropriate `plot_creator`s are
    #       passed. This is also the reason why the visualization callback is not included via
    #       Requires.jl only when Plots is present.
    #       In the future, we should update/remove this warning if other plotting packages are
    #       starting to be used.
    if !(:Plots in names(@__MODULE__, all = true))
        @warn "Package `Plots` not loaded but required by `VisualizationCallback` to visualize results"
    end

    DiscreteCallback(visualization_callback, visualization_callback, # the first one is the condition, the second the affect!
                     save_positions = (false, false),
                     initialize = initialize!)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u, t,
                     integrator) where {Condition, Affect! <: VisualizationCallback}
    visualization_callback = cb.affect!

    visualization_callback(integrator)

    return nothing
end

# this method is called to determine whether the callback should be activated
function (visualization_callback::VisualizationCallback)(u, t, integrator)
    @unpack interval = visualization_callback

    # With error-based step size control, some steps can be rejected. Thus,
    #   `integrator.iter >= integrator.stats.naccept`
    #    (total #steps)       (#accepted steps)
    # We need to check the number of accepted steps since callbacks are not
    # activated after a rejected step.
    return interval > 0 && (integrator.stats.naccept % interval == 0 ||
            isfinished(integrator))
end

# this method is called when the callback is activated
function (visualization_callback::VisualizationCallback)(integrator)
    u_ode = integrator.u
    semi = integrator.p
    mesh, equations, solver, cache = mesh_equations_solver_cache(integrator.p)
    @unpack plot_arguments, solution_variables, variable_names, show_mesh, plot_data_creator, plot_creator = visualization_callback
    if ndims(mesh) == 3 
        plot_data_creator = plot_data_creator == PlotData2D ? PlotData3D : plot_data_creator
        plot_creator = plot_creator == show_plot ? show_plot3D : plot_creator
    end

    # Extract plot data
    plot_data = plot_data_creator(u_ode, semi, solution_variables = solution_variables)

    # If variable names were not specified, plot everything
    if isempty(variable_names)
        variable_names = String[keys(plot_data)...]
    end

    # Create plot
    plot_creator(plot_data, variable_names;
                 show_mesh = show_mesh, plot_arguments = plot_arguments,
                 time = integrator.t, timestep = integrator.stats.naccept)

    # avoid re-evaluating possible FSAL stages
    u_modified!(integrator, false)
    return nothing
end

"""
    show_plot(plot_data, variable_names;
              show_mesh=true, plot_arguments=Dict{Symbol,Any}(),
              time=nothing, timestep=nothing)

Visualize the plot data object provided in `plot_data` and display result, plotting only the
variables in `variable_names` and, optionally, the mesh (if `show_mesh` is `true`).  Additionally,
`plot_arguments` will be unpacked and passed as keyword arguments to the `Plots.plot` command.

This function is the default `plot_creator` argument for the [`VisualizationCallback`](@ref).
`time` and `timestep` are currently unused by this function.

!!! warning "Experimental implementation"
    This is an experimental feature and may change in future releases.

See also: [`VisualizationCallback`](@ref), [`save_plot`](@ref)
"""
function show_plot(plot_data, variable_names;
                   show_mesh = true, plot_arguments = Dict{Symbol, Any}(),
                   time = nothing, timestep = nothing)
    # Gather subplots
    plots = []
    for v in variable_names
        push!(plots, Plots.plot(plot_data[v]; plot_arguments...))
    end
    if show_mesh
        push!(plots, Plots.plot(getmesh(plot_data); plot_arguments...))
    end

    # Note, for the visualization callback to work for general equation systems
    # this layout construction would need to use the if-logic below.
    # Currently, there is no use case for this so it is left here as a note.
    #
    # Determine layout
    if length(plots) <= 3
      cols = length(plots)
      rows = 1
    else
      cols = ceil(Int, sqrt(length(plots)))
      rows = div(length(plots), cols, RoundUp)
    end
    layout = (rows, cols)

    # Determine layout
    cols = ceil(Int, sqrt(length(plots)))
    rows = div(length(plots), cols, RoundUp)
    layout = (rows, cols)

    # Show plot
    display(Plots.plot(plots..., layout = layout))
end

#converts a single int into a tuple of ints, to get a square arrangement for example f(1) = (1,1) f(2) = (2,1) f(3) = (2,2) f(4) = (1,2)
function intTo2DInt(n)
    if n == 1
        return (1,1)
    end
    t = intTo2DInt(n-1)
    if t[1] == 1
        return (t[2] + 1, 1)
    elseif t[1] > t[2]
        return (t[1], t[2] + 1)
    elseif t[2] >= t[1]
        return (t[1] - 1, t[2])
    end
end

function show_plot3D(plot_data, variable_names;
                    show_mesh = false, plot_arguments = Dict{Symbol, Any}(),
                    time = nothing, timestep = nothing)
    # Gather subplots
    # println(size(plot_data))
    # println(plot_data)

    xs = plot_data.x
    ys = plot_data.y
    zs = plot_data.z
    gsx = size(xs)[1]
    gsy = size(ys)[1]
    gsz = size(zs)[1]
    x = [xs[i] for k in 1:gsz for j in 1:gsy for i in 1:gsx]
    y = [ys[j] for k in 1:gsz for j in 1:gsy for i in 1:gsx]
    z = [zs[k] for k in 1:gsz for j in 1:gsy for i in 1:gsx]

    # connectivity = [0, 1, gsx + 1, gsx, gsy * gsx, gsy * gsx + 1, gsy * gsx + gsx + 1, gsy * gsx + gsx]
    #     for c_x in 0:(gsx - 2), c_y in 0:(gsy - 2), c_z in 0:(gsz - 2)
    #         if !(c_x == 0 && c_y == 0 && c_z == 0)
    #             push!(connectivity, ((c_z * gsy * gsx) + (c_y * gsx) + c_x
    #             , (c_z * gsy * gsx) + (c_y * gsx) + c_x + 1
    #             , (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
    #             , (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x
    #             , ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x
    #             , ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x + 1
    #             , ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
    #             , ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x))
    #         end
    #     end
    fig = Figure()
    for v in 1:size(variable_names)[1]
        volume(fig[intTo2DInt(v)...], x, y, z, vec(plot_data.data[v]))
    end

    fig
end

"""
    save_plot(plot_data, variable_names;
              show_mesh=true, plot_arguments=Dict{Symbol,Any}(),
              time=nothing, timestep=nothing)

Visualize the plot data object provided in `plot_data` and save result as a PNG file in the `out`
directory, plotting only the variables in `variable_names` and, optionally, the mesh (if `show_mesh`
is `true`).  Additionally, `plot_arguments` will be unpacked and passed as keyword arguments to the
`Plots.plot` command.

The `timestep` is used in the filename. `time` is currently unused by this function.

!!! warning "Experimental implementation"
    This is an experimental feature and may change in future releases.

See also: [`VisualizationCallback`](@ref), [`show_plot`](@ref)
"""
function save_plot(plot_data, variable_names;
                   show_mesh = true, plot_arguments = Dict{Symbol, Any}(),
                   time = nothing, timestep = nothing)
    # Gather subplots
    plots = []
    for v in variable_names
        push!(plots, Plots.plot(plot_data[v]; plot_arguments...))
    end
    if show_mesh
        push!(plots, Plots.plot(getmesh(plot_data); plot_arguments...))
    end

    # Determine layout
    cols = ceil(Int, sqrt(length(plots)))
    rows = div(length(plots), cols, RoundUp)
    layout = (rows, cols)

    # Create plot
    Plots.plot(plots..., layout = layout)

    # Determine filename and save plot
    filename = joinpath("out", @sprintf("solution_%09d.png", timestep))
    Plots.savefig(filename)
end
end # @muladd
