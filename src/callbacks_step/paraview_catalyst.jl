# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

mutable struct ParaviewCatalystCallback
    interval::Int
    nvisnodes
    catalyst_pipeline
    interpolation
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
            "nvisnodes" => visualization_callback.nvisnodes,
            "catalyst_pipeline" => visualization_callback.catalyst_pipeline,
            "interpolation" => visualization_callback.interpolation,
            
        ]
        summary_box(io, "ParaviewCatalystCallback", setup)
    end
end

"""
    ParaviewCatalystCallback(; interval=0, nvisnodes = nothing, catalyst_pipeline = nothing
                            )

Create a callback that visualizes results during a simulation using Paraview, also known as *in-situ
visualization*. Make sure that your downloaded Paraview Installation includes the Catalyst library and
make sure to set the PARAVIEW_CATALYST_PATH Environment Variable to the path of the Catalyst lib.
You can also specify a path for a custom catalyst pipeline to be used instead of the default by ParaviewCatalyst.jl

!!! warning "Experimental implementation"
    This is an experimental feature and may change in any future releases.
"""
function ParaviewCatalystCallback(; interval = 0, nvisnodes = nothing, catalyst_pipeline = nothing, interpolation = true
                               )
    mpi_isparallel() && error("this callback does not work in parallel yet")

    # ParaviewCatalyst.catalyst_initialize(libpath="/home/nico/Paraview/ParaView-5.13.0-MPI-Linux-Python3.10-x86_64/lib/catalyst")
    if catalyst_pipeline === nothing
        ParaviewCatalyst.catalyst_initialize()
    else
        ParaviewCatalyst.catalyst_initialize(catalyst_pipeline=catalyst_pipeline)
    end

    visualization_callback = ParaviewCatalystCallback(interval, nvisnodes, catalyst_pipeline, interpolation)

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

function create_conduit_node(integrator, mesh::TreeMesh, equations, solver, caches, nvisnodes)
    #creating the conduit node, that will later be passed to paraview
    node = ParaviewCatalyst.ConduitNode() 

    #information about the problem being solved
    nvars = nvariables(equations)
    solution_variables_ = digest_solution_variables(equations, nothing)
    varnames = Trixi.varnames(solution_variables_, equations)
    timestep = integrator.stats.naccept

    #declare variables
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
        pd = PlotData1D(integrator.u, integrator.p, nvisnodes=nvisnodes)
        @trixi_timeit timer() "uniform grid generation" begin
            uniq_ind = unique(i -> pd.x[i], eachindex(pd.x))  # indices of unique elements in pd.x
            c_i = length(uniq_ind)
            x0 = pd.x[1]
            dx = pd.x[uniq_ind[2]] - pd.x[uniq_ind[1]]
        end
    elseif ndims(mesh) == 2
        pd = PlotData2D(integrator.u, integrator.p, nvisnodes=nvisnodes)
        @trixi_timeit timer() "uniform grid generation" begin
        
            c_i = length(pd.x)
            x0 = pd.x[1]
            dx = pd.x[2] - pd.x[1]

            c_j = length(pd.y)
            y0 = pd.y[1]
            dy = pd.y[2] - pd.y[1]
        end
    elseif ndims(mesh) == 3
        pd = PlotData3D(integrator.u, integrator.p; grid_lines=false, nvisnodes=nvisnodes)
        @trixi_timeit timer() "uniform grid generation" begin
            c_i = length(pd.x)
            x0 = pd.x[1]
            dx = pd.x[2] - pd.x[1]

            c_j = length(pd.y)
            y0 = pd.y[1]
            dy = min([pd.y[i + 1] - pd.y[i] for i in 1:(c_j - 1)]...)

            c_k = length(pd.z)
            z0 = pd.z[1]
            dz = min([pd.z[i + 1] - pd.z[i] for i in 1:(c_k - 1)]...)
        end
    end

    #passing general information to the conduit node
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "uniform"

    #telling the conduit node the measurements of the uniform grid
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

    #creating a field for the data from the simulation and passing the data to the conduit node
    @trixi_timeit timer() "passing simulation data to conduit node" begin
        for i in 1:nvars
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/association"] = "vertex"
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/topology"] = "mesh"
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/volume_dependent"] = "false"
            if ndims(mesh) == 1
                node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = pd.data[i]
            elseif ndims(mesh) == 2
                node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(pd.data[i])
            elseif ndims(mesh) == 3
                node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(pd.data[i])
            end
        end
    end
    
    return node
end

function create_conduit_node(integrator, mesh::P4estMesh, equations, solver, cache, nvisnodes)
    #creating the conduit node, that will later be passed to paraview
    node = ParaviewCatalyst.ConduitNode()

    #information about the problem being solved
    n_visnodes = (nvisnodes === nothing) ? 2 * nnodes(solver) : nvisnodes
    ndims_ = ndims(mesh)
    nvars = nvariables(equations)
    solution_variables_ = digest_solution_variables(equations, nothing)
    varnames = Trixi.varnames(solution_variables_, equations)
    timestep = integrator.stats.naccept

    
    #node data, location and interpolation
    @trixi_timeit timer() "node/data interpolation" begin
        nodes = collect(range(-1, 1, length=n_visnodes))
        node_coordinates = Array{Float64, ndims_+2}(undef, ndims_, ntuple(_ -> n_visnodes, ndims_)..., Trixi.ncells(mesh))
        interpolation_node_coordinates = calc_node_coordinates!(node_coordinates, mesh, nodes)
        unstructured_data = get_unstructured_data(Trixi.wrap_array(integrator.u, integrator.p), solution_variables_, mesh, equations, solver, cache)
        interpolated_data = raw2interpolated(unstructured_data, nodes)
    end

    #information about the topology underneath the data
    grid_size = size(interpolation_node_coordinates)
    gsx = grid_size[2]
    gsy = grid_size[3]
    gsz =(ndims_ == 3) ? grid_size[4] : nothing
    gstr =(ndims_ == 3) ? grid_size[5] : grid_size[4]

    #passing time and position data over to the conduit node
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "explicit"
    x = (ndims_ == 2) ? vec(interpolation_node_coordinates[1, :, :, :]) : vec(interpolation_node_coordinates[1, :, :, :, :])
    node["catalyst/channels/input/data/coordsets/coords/values/x"] = x
    y = (ndims_ == 2) ? vec(interpolation_node_coordinates[2, :, :, :]) : vec(interpolation_node_coordinates[2, :, :, :, :])
    node["catalyst/channels/input/data/coordsets/coords/values/y"] = y
    z = nothing
    if ndims_ == 3
        z = vec(interpolation_node_coordinates[3, :, :, :, :])
        node["catalyst/channels/input/data/coordsets/coords/values/z"] = z
    end

    #creating a topology
    node["catalyst/channels/input/data/topologies/mesh/type"] = "unstructured"
    node["catalyst/channels/input/data/topologies/mesh/coordset"] = "coords"
    node["catalyst/channels/input/data/topologies/mesh/elements/shape"] = (ndims_ == 2) ? "quad" : "hex"
    if ndims(mesh) == 2
        #The array in the form of [[tree1_Cell1_lower_left_corner, tree1_Cell1_upper_left_corner, tree1_Cell1_upper_right_corner, tree1_Cell1_lower_right_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
        @trixi_timeit timer() "generating connectivity array" begin
            node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
                [(c_tr * gsy * gsx) + (c_y * gsx) + c_x
                (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x
                (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
                (c_tr * gsy * gsx) + (c_y * gsx) + c_x + 1]
                for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_tr in 0:(gstr - 1)]...), :)
        end
    else
        #The array in the form of [[tree1_Cell1_lower_left_front_corner, tree1_Cell1_lower_right_front_corner, tree1_Cell1_upper_right_front_corner, tree1_Cell1_upper_left_front_corner, tree1_Cell1_lower_left_back_corner, tree1_Cell1_lower_right_back_corner, tree1_Cell1_upper_right_back_corner, tree1_Cell1_upper_left_back_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
        @trixi_timeit timer() "generating connectivity array" begin
        connectivity = [0, 1, gsx + 1, gsx, gsy * gsx, gsy * gsx + 1, gsy * gsx + gsx + 1, gsy * gsx + gsx]
        for c_x in 0:(gsx - 2), c_y in 0:(gsy - 2), c_z in 0:(gsz - 2), c_tr in 0:(gstr - 1)
            if !(c_x == 0 && c_y == 0 && c_z == 0 && c_tr == 0)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x + 1)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x + 1)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1)
                push!(connectivity, (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x)
            end
        end
        node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = connectivity
            # node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
            #     [(c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x
            #     (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x + 1
            #     (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
            #     (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x
            #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x
            #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x + 1
            #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
            #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x]
            #     for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_z in 0:(gsz - 2) for c_tr in 0:(gstr - 1)]...), :)
        end
    end

    #passing the simulation solution data for each variable to the conduit node
    for i in 1:nvars
        node["catalyst/channels/input/data/fields/" * varnames[i] * "/association"] = "vertex"
        node["catalyst/channels/input/data/fields/" * varnames[i] * "/topology"] = "mesh"
        node["catalyst/channels/input/data/fields/" * varnames[i] * "/volume_dependent"] = "false"
        if ndims(mesh) == 2
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(interpolated_data[:, i])
        else
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(interpolated_data[:, i])
        end
    end

    #returns the completed conduit node
    return node
end

function create_conduit_node_no_interpolation(integrator, mesh::TreeMesh, equations, solver, cache)
    #creating the conduit node, that will later be passed to paraview
    node = ParaviewCatalyst.ConduitNode() 

    #information about the problem being solved
    nvars = nvariables(equations)
    solution_variables_ = digest_solution_variables(equations, nothing)
    varnames = Trixi.varnames(solution_variables_, equations)
    timestep = integrator.stats.naccept

    #node data and location
    leaf_cell_ids = leaf_cells(mesh.tree)
    unstructured_data = get_unstructured_data(Trixi.wrap_array(integrator.u, integrator.p), solution_variables_, mesh, equations, solver, cache)
    uds = size(unstructured_data)
    reshaped_data = (ndims(mesh) == 3) ? reshape(vcat([
        unstructured_data[i,j,k,l,m]
        for m in 1:uds[5] for k in uds[3] for j in 1:uds[2] for l in 1:uds[4] for i in 1:uds[1]
    ]...), :) : reshape(vcat([
        unstructured_data[i,j,k,l]
        for l in 1:uds[4] for j in 1:uds[2] for k in 1:uds[3] for i in 1:uds[1]
    ]...), :)
    coordinates = mesh.tree.coordinates[:, leaf_cell_ids]
    levels = mesh.tree.levels[leaf_cell_ids]
    length_level_0 = mesh.tree.length_level_0
    center_level_0 = mesh.tree.center_level_0
    max_level = maximum(levels)
    max_nvisnodes = nnodes(solver)
    resolution = max_nvisnodes * 2^max_level
    xs = collect(range(-1, 1, length = resolution + 1)) .* length_level_0 / 2 .+
         center_level_0[1]
    ys = collect(range(-1, 1, length = resolution + 1)) .* length_level_0 / 2 .+
         center_level_0[2]
    zs = (ndims(mesh) == 3) ? collect(range(-1, 1, length = resolution + 1)) .* length_level_0 / 2 .+
         center_level_0[3] : nothing

    #information about the topology underneath the data
    gsx = size(xs)[1]
    gsy = size(ys)[1]
    gsz =(ndims(mesh) == 3) ? size(zs)[1] : nothing

    #passing general information to the conduit node
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "explicit"

    if ndims(mesh) == 3
        x = [xs[i] for k in 1:gsz for j in 1:gsy for i in 1:gsx]
        node["catalyst/channels/input/data/coordsets/coords/values/x"] = x
        y = [ys[j] for k in 1:gsz for j in 1:gsy for i in 1:gsx]
        node["catalyst/channels/input/data/coordsets/coords/values/y"] = y
        z = [zs[k] for k in 1:gsz for j in 1:gsy for i in 1:gsx]
        node["catalyst/channels/input/data/coordsets/coords/values/z"] = z
    else
        x = [xs[i] for j in 1:gsy for i in 1:gsx]
        node["catalyst/channels/input/data/coordsets/coords/values/x"] = x
        y = [ys[j] for j in 1:gsy for i in 1:gsx]
        node["catalyst/channels/input/data/coordsets/coords/values/y"] = y
    end

    #creating a topology
    node["catalyst/channels/input/data/topologies/mesh/type"] = "unstructured"
    node["catalyst/channels/input/data/topologies/mesh/coordset"] = "coords"
    node["catalyst/channels/input/data/topologies/mesh/elements/shape"] = (ndims(mesh) == 2) ? "quad" : "hex"
    if ndims(mesh) == 2
        #The array in the form of [[tree1_Cell1_lower_left_corner, tree1_Cell1_upper_left_corner, tree1_Cell1_upper_right_corner, tree1_Cell1_lower_right_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
        @trixi_timeit timer() "generating connectivity array" begin
            node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
                [(c_y * gsx) + c_x
                ((c_y + 1) * gsx) + c_x
                ((c_y + 1) * gsx) + c_x + 1
                (c_y * gsx) + c_x + 1]
                for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2)]...), :)
        end
    else
        #The array in the form of [[tree1_Cell1_lower_left_front_corner, tree1_Cell1_lower_right_front_corner, tree1_Cell1_upper_right_front_corner, tree1_Cell1_upper_left_front_corner, tree1_Cell1_lower_left_back_corner, tree1_Cell1_lower_right_back_corner, tree1_Cell1_upper_right_back_corner, tree1_Cell1_upper_left_back_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
        @trixi_timeit timer() "generating connectivity array" begin
            node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
                [(c_z * gsy * gsx) + (c_y * gsx) + c_x
                (c_z * gsy * gsx) + (c_y * gsx) + c_x + 1
                (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
                (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x
                ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x
                ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x + 1
                ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
                ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x]
                for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_z in 0:(gsz - 2)]...), :)
        end
    end

    #creating a field for the data from the simulation and passing the data to the conduit node
    @trixi_timeit timer() "passing simulation data to conduit node" begin
        for i in 1:nvars
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/association"] = "element"
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/topology"] = "mesh"
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/volume_dependent"] = "false"
            if ndims(mesh) == 2
                node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(unstructured_data[:, :, :, i])
            elseif ndims(mesh) == 3
                println("soll:" * string(size(x)))
                println(size(vec(unstructured_data[:, :, :, :, i])))
                node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(unstructured_data[:, :, :, :, i])
            end
        end
    end
    
    return node
end

function create_conduit_node_no_interpolation(integrator, mesh::P4estMesh, equations, solver, cache)
    #creating the conduit node, that will later be passed to paraview
    node = ParaviewCatalyst.ConduitNode()

    #information about the problem being solved
    n_visnodes = nnodes(solver)
    ndims_ = ndims(mesh)
    nvars = nvariables(equations)
    solution_variables_ = digest_solution_variables(equations, nothing)
    varnames = Trixi.varnames(solution_variables_, equations)
    timestep = integrator.stats.naccept

    
    #node data, location and interpolation
    @trixi_timeit timer() "node/data interpolation" begin
        node_coordinates = mesh.tree_node_coordinates
        # node_coordinates = Array{Float64, ndims_+2}(undef, ndims_, ntuple(_ -> n_visnodes, ndims_)..., Trixi.ncells(mesh))
        unstructured_data = get_unstructured_data(Trixi.wrap_array(integrator.u, integrator.p), solution_variables_, mesh, equations, solver, cache)
    end

    #information about the topology underneath the data
    grid_size = size(node_coordinates)
    gsx = grid_size[2]
    gsy = grid_size[3]
    gsz =(ndims_ == 3) ? grid_size[4] : nothing
    gstr =(ndims_ == 3) ? grid_size[5] : grid_size[4]

    #passing time and position data over to the conduit node
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "explicit"
    x = (ndims_ == 2) ? vec(node_coordinates[1, :, :, :]) : vec(node_coordinates[1, :, :, :, :])
    node["catalyst/channels/input/data/coordsets/coords/values/x"] = x
    y = (ndims_ == 2) ? vec(node_coordinates[2, :, :, :]) : vec(node_coordinates[2, :, :, :, :])
    node["catalyst/channels/input/data/coordsets/coords/values/y"] = y
    z = nothing
    if ndims_ == 3
        z = vec(node_coordinates[3, :, :, :, :])
        node["catalyst/channels/input/data/coordsets/coords/values/z"] = z
    end

    #creating a topology
    node["catalyst/channels/input/data/topologies/mesh/type"] = "unstructured"
    node["catalyst/channels/input/data/topologies/mesh/coordset"] = "coords"
    node["catalyst/channels/input/data/topologies/mesh/elements/shape"] = (ndims_ == 2) ? "quad" : "hex"
    if ndims(mesh) == 2
        #The array in the form of [[tree1_Cell1_lower_left_corner, tree1_Cell1_upper_left_corner, tree1_Cell1_upper_right_corner, tree1_Cell1_lower_right_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
        @trixi_timeit timer() "generating connectivity array" begin
            node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
                [(c_tr * gsy * gsx) + (c_y * gsx) + c_x
                (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x
                (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
                (c_tr * gsy * gsx) + (c_y * gsx) + c_x + 1]
                for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_tr in 0:(gstr - 1)]...), :)
        end
    else
        #The array in the form of [[tree1_Cell1_lower_left_front_corner, tree1_Cell1_lower_right_front_corner, tree1_Cell1_upper_right_front_corner, tree1_Cell1_upper_left_front_corner, tree1_Cell1_lower_left_back_corner, tree1_Cell1_lower_right_back_corner, tree1_Cell1_upper_right_back_corner, tree1_Cell1_upper_left_back_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
        @trixi_timeit timer() "generating connectivity array" begin
            node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
                [(c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x
                (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x + 1
                (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
                (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x
                (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x
                (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x + 1
                (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
                (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x]
                for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_z in 0:(gsz - 2) for c_tr in 0:(gstr - 1)]...), :)
        end
    end

    #passing the simulation solution data for each variable to the conduit node
    for i in 1:nvars
        node["catalyst/channels/input/data/fields/" * varnames[i] * "/association"] = "vertex"
        node["catalyst/channels/input/data/fields/" * varnames[i] * "/topology"] = "mesh"
        node["catalyst/channels/input/data/fields/" * varnames[i] * "/volume_dependent"] = "false"
        if ndims(mesh) == 2
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(unstructured_data[:, :, :, i])
        else
            node["catalyst/channels/input/data/fields/" * varnames[i] * "/values"] = vec(unstructured_data[:, :, :, :, i])
        end
    end

    #returns the completed conduit node
    return node
end

# this method is called when the callback is activated
function (visualization_callback::ParaviewCatalystCallback)(integrator)
    #extracting problem data from the integrator
    mesh, equations, solver, cache = mesh_equations_solver_cache(integrator.p)

    # avoid re-evaluating possible FSAL stages
    u_modified!(integrator, false)
    # Conduit.node_info(node) do info_node
    #    Conduit.node_print(info_node, detailed = true)
    # end
    node = (visualization_callback.interpolation ? 
    create_conduit_node(integrator, mesh, equations, solver, cache, visualization_callback.nvisnodes) : 
    create_conduit_node_no_interpolation(integrator, mesh, equations, solver, cache))
    #passing the problem data to a function fitting the right mesh, to create a conduit node, which is then passed to paraview
    @trixi_timeit timer() "catalyst execute" begin
        ParaviewCatalyst.catalyst_execute(node)
    end

    return nothing
end 

#Stolen Code from Trixi2Vtk

function raw2interpolated(data_gl::AbstractArray{Float64}, nodes_out)
    # Extract number of spatial dimensions
    ndims_ = ndims(data_gl) - 2
  
    # Extract data shape information
    n_nodes_in = size(data_gl, 1)
    n_nodes_out = length(nodes_out)
    n_elements = size(data_gl, ndims_ + 1)
    n_variables = size(data_gl, ndims_ + 2)
  
    # Get node coordinates for DG locations on reference element
    nodes_in, _ = gauss_lobatto_nodes_weights(n_nodes_in)
  
    # Calculate Vandermonde matrix
    vandermonde = polynomial_interpolation_matrix(nodes_in, nodes_out)
  
    if ndims_ == 2
      # Create output data structure
      data_vis = Array{Float64}(undef, n_nodes_out, n_nodes_out, n_elements, n_variables)
  
      # For each variable, interpolate element data and store to global data structure
      for v in 1:n_variables
        # Reshape data array for use in interpolate_nodes function
        @views reshaped_data = reshape(data_gl[:, :, :, v], 1, n_nodes_in, n_nodes_in, n_elements)
  
        # Interpolate data to visualization nodes
        for element_id in 1:n_elements
          @views data_vis[:, :, element_id, v] .= reshape(
              interpolate_nodes(reshaped_data[:, :, :, element_id], vandermonde, 1),
              n_nodes_out, n_nodes_out)
        end
      end
    elseif ndims_ == 3
      # Create output data structure
      data_vis = Array{Float64}(undef, n_nodes_out, n_nodes_out, n_nodes_out, n_elements, n_variables)
  
      # For each variable, interpolate element data and store to global data structure
      for v in 1:n_variables
        # Reshape data array for use in interpolate_nodes function
        @views reshaped_data = reshape(data_gl[:, :, :, :, v],
                                       1, n_nodes_in, n_nodes_in, n_nodes_in, n_elements)
  
        # Interpolate data to visualization nodes
        for element_id in 1:n_elements
          @views data_vis[:, :, :, element_id, v] .= reshape(
              interpolate_nodes(reshaped_data[:, :, :, :, element_id], vandermonde, 1),
              n_nodes_out, n_nodes_out, n_nodes_out)
        end
      end
    else
      error("Unsupported number of spatial dimensions: ", ndims_)
    end
  
    # Return as one 1D array for each variable
    return reshape(data_vis, n_nodes_out^ndims_ * n_elements, n_variables)
end

# Interpolate data using the given Vandermonde matrix and return interpolated values (3D version).
function interpolate_nodes(data_in::AbstractArray{T, 4},
    vandermonde, n_vars) where T
    n_nodes_out = size(vandermonde, 1)
    data_out = zeros(eltype(data_in), n_vars, n_nodes_out, n_nodes_out, n_nodes_out)
    interpolate_nodes!(data_out, data_in, vandermonde, n_vars)
end

function interpolate_nodes!(data_out::AbstractArray{T, 4}, data_in::AbstractArray{T, 4},
    vandermonde, n_vars) where T
    n_nodes_out = size(vandermonde, 1)
    n_nodes_in  = size(vandermonde, 2)

    for k in 1:n_nodes_out, j in 1:n_nodes_out, i in 1:n_nodes_out
        for v in 1:n_vars
            acc = zero(eltype(data_out))
            for kk in 1:n_nodes_in, jj in 1:n_nodes_in, ii in 1:n_nodes_in
                acc += vandermonde[i, ii] * vandermonde[j, jj] * vandermonde[k, kk] * data_in[v, ii, jj, kk]    
            end
            data_out[v, i, j, k] = acc
        end
    end

    return data_out
end

# Interpolate data using the given Vandermonde matrix and return interpolated values (2D version).
function interpolate_nodes(data_in::AbstractArray{T, 3},
    vandermonde, n_vars) where T
    n_nodes_out = size(vandermonde, 1)
    data_out = zeros(eltype(data_in), n_vars, n_nodes_out, n_nodes_out)
    interpolate_nodes!(data_out, data_in, vandermonde, n_vars)
end


function interpolate_nodes!(data_out::AbstractArray{T, 3}, data_in::AbstractArray{T, 3},
     vandermonde, n_vars) where T
    n_nodes_out = size(vandermonde, 1)
    n_nodes_in  = size(vandermonde, 2)

    for j in 1:n_nodes_out
        for i in 1:n_nodes_out
            for v in 1:n_vars
                acc = zero(eltype(data_out))
                for jj in 1:n_nodes_in
                    for ii in 1:n_nodes_in
                        acc += vandermonde[i, ii] * data_in[v, ii, jj] * vandermonde[j, jj]
                    end
                end
                data_out[v, i, j] = acc
            end
        end
    end

    return data_out
end

end # @muladd
