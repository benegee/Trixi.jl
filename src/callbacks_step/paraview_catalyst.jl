# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

mutable struct ParaviewCatalystCallback
    interval::Int
end

function Base.show(io::IO,
                   cb::DiscreteCallback{Condition, Affect!}) where {Condition,
                                                                    Affect! <:
                                                                    ParaviewCatalystCallback
                                                                    }
    visualization_callback = cb.affect!
    @unpack interval = visualization_callback
    print(io, "ParaviewCatalystCallback(",
          "interval=", interval,")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{Condition, Affect!}) where {Condition,
                                                                    Affect! <:
                                                                    ParaviewCatalystCallback
                                                                    }
    if get(io, :compact, false)
        show(io, cb)
    else
        visualization_callback = cb.affect!

        setup = [
            "interval" => visualization_callback.interval,
            
        ]
        summary_box(io, "ParaviewCatalystCallback", setup)
    end
end

"""
    ParaviewCatalystCallback(; interval=0,
                            )

Create a callback that visualizes results during a simulation, also known as *in-situ
visualization*. Make sure to set the PARAVIEW_CATALYST_PATH Environment Variable to the path of the Catalyst lib.

!!! warning "Experimental implementation"
    This is an experimental feature and may change in any future releases.
"""
function ParaviewCatalystCallback(; interval = 0, catalyst_pipeline=nothing
                               )
    mpi_isparallel() && error("this callback does not work in parallel yet")

    # ParaviewCatalyst.catalyst_initialize(libpath="/home/nico/Paraview/ParaView-5.13.0-MPI-Linux-Python3.10-x86_64/lib/catalyst")
    if catalyst_pipeline === nothing
        ParaviewCatalyst.catalyst_initialize()
    else
        ParaviewCatalyst.catalyst_initialize(catalyst_pipeline=catalyst_pipeline)
    end

    visualization_callback = ParaviewCatalystCallback(interval)

    # Warn users if they create a ParaviewCatalystCallback without having loaded the ParaviewCatalyst package
    if !(:ParaviewCatalyst in nameof.(Base.loaded_modules |> values))
        @warn "Package `ParaviewCatalyst` not loaded but required by `ParaviewCatalystCallback` to visualize results"
    end

    DiscreteCallback(visualization_callback, visualization_callback, # the first one is the condition, the second the affect!
                     save_positions = (false, false),
                     initialize = initialize!)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u, t,
                     integrator) where {Condition, Affect! <: ParaviewCatalystCallback}
    visualization_callback = cb.affect!

    visualization_callback(integrator)

    return nothing
end

# this method is called to determine whether the callback should be activated
function (visualization_callback::ParaviewCatalystCallback)(u, t, integrator)
    @unpack interval = visualization_callback

    # With error-based step size control, some steps can be rejected. Thus,
    #   `integrator.iter >= integrator.stats.naccept`
    #    (total #steps)       (#accepted steps)
    # We need to check the number of accepted steps since callbacks are not
    # activated after a rejected step.
    return interval > 0 && (integrator.stats.naccept % interval == 0 ||
            isfinished(integrator))
end

function create_conduit_node(integrator, mesh::TreeMesh)
    node = ParaviewCatalyst.ConduitNode()
    timestep = integrator.stats.naccept
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "uniform"

    #deklare variables
    pd = nothing
    c_i = 0
    c_j = 0
    c_k = 0
    x0 = 0
    y0 = 0
    z0 = 0
    dx = 0
    dy = 0
    dz = 0

    #get the data to be plotted via a PlotData function corresponding to the dimension. Then determine point count (c_i, c_j, c_k), start values (x0, y0, z0) and step size (dx, dy, dz) for each dimension
    if ndims(mesh) == 1
        pd = PlotData1D(integrator.u, integrator.p)
        c_i = length(pd.x)
        x0 = min(pd.x...)
        uniq_ind = unique(i -> pd.x[i], eachindex(pd.x))  # indices of unique elements in pd.x
        dx = min([pd.x[uniq_ind[i + 1]] - pd.x[uniq_ind[i]] for i in 1:(length(uniq_ind) - 1)]...)
    elseif ndims(mesh) == 2
        pd = PlotData2D(integrator.u, integrator.p)
        
        c_i = length(pd.x)
        x0 = min(pd.x...)
        dx = min([pd.x[i + 1] - pd.x[i] for i in 1:(c_i - 1)]...)

        c_j = length(pd.y)
        y0 = min(pd.y...)
        dy = min([pd.y[i + 1] - pd.y[i] for i in 1:(c_j - 1)]...)
    elseif ndims(mesh) == 3
        pd = PlotData3D(integrator.u, integrator.p; grid_lines=false)
        c_i = length(pd.x)
        x0 = min(pd.x...)
        dx = min([pd.x[i + 1] - pd.x[i] for i in 1:(c_i - 1)]...)

        c_j = length(pd.y)
        y0 = min(pd.y...)
        dy = min([pd.y[i + 1] - pd.y[i] for i in 1:(c_j - 1)]...)

        c_k = length(pd.z)
        z0 = min(pd.z...)
        dz = min([pd.z[i + 1] - pd.z[i] for i in 1:(c_k - 1)]...)
    end

    #telling catalyst the measurements of the uniform grid
    node["catalyst/channels/input/data/coordsets/coords/dims/i"] = c_i
    node["catalyst/channels/input/data/coordsets/coords/origin/x"] = x0
    node["catalyst/channels/input/data/coordsets/coords/spacing/dx"] = dx
    if ndims(mesh) > 1
        node["catalyst/channels/input/data/coordsets/coords/dims/j"] = c_j
        node["catalyst/channels/input/data/coordsets/coords/origin/y"] = y0
        node["catalyst/channels/input/data/coordsets/coords/spacing/dy"] = dy
        if ndims(mesh) > 2
            node["catalyst/channels/input/data/coordsets/coords/dims/k"] = c_k
            node["catalyst/channels/input/data/coordsets/coords/origin/z"] = z0
            node["catalyst/channels/input/data/coordsets/coords/spacing/dz"] = dz
        end
    end

    #creating a topology
    node["catalyst/channels/input/data/topologies/mesh/type"] = "uniform"
    node["catalyst/channels/input/data/topologies/mesh/coordset"] = "coords"

    #creating a field for the data from the simulation and passing the data to catalyst
    node["catalyst/channels/input/data/fields/solution/association"] = "vertex"
    node["catalyst/channels/input/data/fields/solution/topology"] = "mesh"
    node["catalyst/channels/input/data/fields/solution/volume_dependent"] = "false"
    if ndims(mesh) == 1
        node["catalyst/channels/input/data/fields/solution/values"] = pd.data
    elseif ndims(mesh) == 2
        node["catalyst/channels/input/data/fields/solution/values"] = vec(pd.data[1])
    elseif ndims(mesh) == 3
        node["catalyst/channels/input/data/fields/solution/values"] = vec(pd.data[1])
    end
    return node
end

function create_conduit_node(integrator, mesh::P4estMesh)
    node = ParaviewCatalyst.ConduitNode()
    timestep = integrator.stats.naccept
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "uniform"

    vtk_points, vtk_cells = calc_vtk_points_cells(mesh.tree_node_coordinates)
    println()
    println(vtk_points)
    println()
    println(vtk_cells)
    println()

    #creating a topology
    node["catalyst/channels/input/data/topologies/mesh/type"] = "uniform"
    node["catalyst/channels/input/data/topologies/mesh/coordset"] = "coords"

    #creating a field for the data from the simulation and passing the data to catalyst
    node["catalyst/channels/input/data/fields/solution/association"] = "vertex"
    node["catalyst/channels/input/data/fields/solution/topology"] = "mesh"
    node["catalyst/channels/input/data/fields/solution/volume_dependent"] = "false"
    if ndims(mesh) == 1
        node["catalyst/channels/input/data/fields/solution/values"] = pd.data
    elseif ndims(mesh) == 2
        println()
        println(pd)
        println()
        println(pd.data[1])
        println()
        println(vec(pd.data[1]))
        println()
        node["catalyst/channels/input/data/fields/solution/values"] = vec(pd.data[1])
    elseif ndims(mesh) == 3
        node["catalyst/channels/input/data/fields/solution/values"] = vec(pd.data[1])
    end
    return node
end

# this method is called when the callback is activated
function (visualization_callback::ParaviewCatalystCallback)(integrator)
    u_ode = integrator.u
    mesh, equations, solver, cache = mesh_equations_solver_cache(integrator.p)
    time = integrator.t
    timestep = integrator.stats.naccept

    println()
    println("*** Catalyst Callback activated")
    println("*** Time ", time)
    println("*** Step ", timestep)
    println("*** u[1] ", u_ode[1])
    println()

    # avoid re-evaluating possible FSAL stages
    u_modified!(integrator, false)
    # Conduit.node_info(node) do info_node
    #    Conduit.node_print(info_node, detailed = true)
    # end
    ParaviewCatalyst.catalyst_execute(create_conduit_node(integrator, mesh))

    return nothing
end

# Copy from Trixi2Vtk
# Convert coordinates and level information to a list of points and VTK cells for `StructuredMesh` (2D version)
function calc_vtk_points_cells(node_coordinates::AbstractArray{<:Any,4})
    n_elements = size(node_coordinates, 4)
    size_ = size(node_coordinates)
    n_points = prod(size_[2:end])
    # Linear indices to access points by node indices and element id
    linear_indices = LinearIndices(size_[2:end])
  
    # Use lagrange nodes as VTK points
    vtk_points = reshape(node_coordinates, (2, n_points))
    vtk_cells = Vector{Vector{Float64}}(undef, n_elements)
  
    # Create cell for each element
    for element in 1:n_elements
      vertices = [linear_indices[1, 1, element],
                  linear_indices[end, 1, element],
                  linear_indices[end, end, element],
                  linear_indices[1, end, element]]
  
      edges = vcat(linear_indices[2:end-1, 1, element],
                   linear_indices[end, 2:end-1, element],
                   linear_indices[2:end-1, end, element],
                   linear_indices[1, 2:end-1, element])
  
      faces = vec(linear_indices[2:end-1, 2:end-1, element])
  
      point_ids = vcat(vertices, edges, faces)
      vtk_cells[element] = node_coordinates[point_ids]
    end
  
    return vtk_points, vtk_cells
  end

# Copy from Trixi2Vtk
# Convert coordinates and level information to a list of points and VTK cells for `StructuredMesh` (3D version)
function calc_vtk_points_cells(node_coordinates::AbstractArray{<:Any,5})
    n_elements = size(node_coordinates, 5)
    size_ = size(node_coordinates)
    n_points = prod(size_[2:end])
    # Linear indices to access points by node indices and element id
    linear_indices = LinearIndices(size_[2:end])
  
    # Use lagrange nodes as VTK points
    vtk_points = reshape(node_coordinates, (3, n_points))
    vtk_cells = Vector{Vector{Float64}}(undef, n_elements)
  
    # Create cell for each element
    for element in 1:n_elements
      vertices = [linear_indices[1, 1, 1, element],
                  linear_indices[end, 1, 1, element],
                  linear_indices[end, end, 1, element],
                  linear_indices[1, end, 1, element],
                  linear_indices[1, 1, end, element],
                  linear_indices[end, 1, end, element],
                  linear_indices[end, end, end, element],
                  linear_indices[1, end, end, element]]
  
      # This order doesn't make any sense. This is completely different
      # from what is shown in
      # https://blog.kitware.com/wp-content/uploads/2018/09/Source_Issue_43.pdf
      # but this is the way it works.
      edges = vcat(linear_indices[2:end-1, 1, 1, element],
                   linear_indices[end, 2:end-1, 1, element],
                   linear_indices[2:end-1, end, 1, element],
                   linear_indices[1, 2:end-1, 1, element],
                   linear_indices[2:end-1, 1, end, element],
                   linear_indices[end, 2:end-1, end, element],
                   linear_indices[2:end-1, end, end, element],
                   linear_indices[1, 2:end-1, end, element],
                   linear_indices[1, 1, 2:end-1, element],
                   linear_indices[end, 1, 2:end-1, element],
                   linear_indices[1, end, 2:end-1, element],
                   linear_indices[end, end, 2:end-1, element])
  
      # See above
      faces = vcat(vec(linear_indices[1, 2:end-1, 2:end-1, element]),
                   vec(linear_indices[end, 2:end-1, 2:end-1, element]),
                   vec(linear_indices[2:end-1, 1, 2:end-1, element]),
                   vec(linear_indices[2:end-1, end, 2:end-1, element]),
                   vec(linear_indices[2:end-1, 2:end-1, 1, element]),
                   vec(linear_indices[2:end-1, 2:end-1, end, element]))
  
      volume = vec(linear_indices[2:end-1, 2:end-1, 2:end-1, element])
  
      point_ids = vcat(vertices, edges, faces, volume)
      vtk_cells[element] = node_coordinates[point_ids]
    end
  
    return vtk_points, vtk_cells
end
  

end # @muladd
