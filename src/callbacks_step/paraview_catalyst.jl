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

function create_conduit_node(integrator, mesh::TreeMesh, equations, solver, cache)
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

function create_conduit_node(integrator, mesh::P4estMesh, equations, solver, cache)
    node = ParaviewCatalyst.ConduitNode()
    timestep = integrator.stats.naccept
    node["catalyst/state/timestep"] = timestep
    node["catalyst/state/time"] = timestep
    node["catalyst/channels/input/type"] = "mesh"
    node["catalyst/channels/input/data/coordsets/coords/type"] = "explicit"

    solution_variables_ = digest_solution_variables(equations, nothing)

    u = Trixi.wrap_array(integrator.u, integrator.p)
    # unstructured_data = get_unstructured_data(u, solution_variables_, mesh, equations,
    # solver, cache)
    # println()
    # println(unstructured_data[:, : , : , : , 1])
    # println(size(unstructured_data))
    # println()
    ndims_ = ndims(mesh)
    n_visnodes = 2 * nnodes(solver)
    nodes, _ = gauss_lobatto_nodes_weights(n_visnodes)
    nvars = nvariables(equations)

    node_coordinates = Array{Float64, ndims_+2}(undef, ndims_,
                                              ntuple(_ -> n_visnodes, ndims_)...,
                                              Trixi.ncells(mesh))

    interpolation_node_coordinates = calc_node_coordinates!(node_coordinates, mesh, nodes)
    interpolated_data = interpolate_data(u, mesh, 4)
    data = vec(interpolated_data)
    grid_size = size(interpolation_node_coordinates)
    gsx = grid_size[2]
    gsy = grid_size[3]
    gsz =(ndims_ == 3) ? grid_size[4] : nothing
    gstr =(ndims_ == 3) ? grid_size[5] : grid_size[4]
    # println()
    # println(size(interpolation_node_coordinates))
    # println(interpolation_node_coordinates)
    # println()

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
        node["catalyst/channels/input/data/topologies/mesh/elements/connectivity"] = reshape(vcat([
            [(c_tr * gsy * gsx) + (c_y * gsx) + c_x
            (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x
            (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
            (c_tr * gsy * gsx) + (c_y * gsx) + c_x + 1]
            for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_tr in 0:(gstr - 1)]...), :)
    else
        #The array in the form of [[tree1_Cell1_lower_left_front_corner, tree1_Cell1_lower_right_front_corner, tree1_Cell1_upper_right_front_corner, tree1_Cell1_upper_left_front_corner, tree1_Cell1_lower_left_back_corner, tree1_Cell1_lower_right_back_corner, tree1_Cell1_upper_right_back_corner, tree1_Cell1_upper_left_back_corner], ... (iterating first over cells, then trees)]
        #gets reshaped to a 1D Array
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
    #creating a field for the data from the simulation and passing the data to catalyst
    # node["catalyst/channels/input/data/fields/solution/association"] = "element"
    node["catalyst/channels/input/data/fields/solution/association"] = "vertex"
    node["catalyst/channels/input/data/fields/solution/topology"] = "mesh"
    node["catalyst/channels/input/data/fields/solution/volume_dependent"] = "false"
    if ndims(mesh) == 2
        node["catalyst/channels/input/data/fields/solution/values"] = [sin(i) for i in 1:(size(x)[1])] #TODO: echte DatenS
        # println()
        # println("soll:" * string(size(x)[1]) * "(vertex)")
        # println("soll:" * string(size(reshape(vcat([
        #     [(c_tr * gsy * gsx) + (c_y * gsx) + c_x
        #     (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x
        #     (c_tr * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
        #     (c_tr * gsy * gsx) + (c_y * gsx) + c_x + 1]
        #     for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_tr in 0:(gstr - 1)]...), :))[1]/4) * "(element)")
        # println(size(data))
        # println()
        # node["catalyst/channels/input/data/fields/solution/values"] = data
    else
        node["catalyst/channels/input/data/fields/solution/values"] = [sin(i) for i in 1:(size(x)[1])]
        # println()
        # println("soll:" * string(size(x)[1]) * "(vertex)")
        # println("soll:" * string(size(reshape(vcat([
        #     [(c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x
        #     (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + (c_y * gsx) + c_x + 1
        #     (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
        #     (c_tr * gsz * gsy * gsx) + (c_z * gsy * gsx) + ((c_y + 1) * gsx) + c_x
        #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x
        #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + (c_y * gsx) + c_x + 1
        #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x + 1
        #     (c_tr * gsz * gsy * gsx) + ((c_z + 1) * gsy * gsx) + ((c_y + 1) * gsx) + c_x]
        #     for c_x in 0:(gsx - 2) for c_y in 0:(gsy - 2) for c_z in 0:(gsz - 2) for c_tr in 0:(gstr - 1)]...), :))[1]/8) * "(element)")
        # println(size(data))
        # println()
        # node["catalyst/channels/input/data/fields/solution/values"] = data
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
    ParaviewCatalyst.catalyst_execute(create_conduit_node(integrator, mesh, equations, solver, cache))

    return nothing
end 

#Stolen Code from Trixi2Vtk

# Interpolate data from input format to desired output format (StructuredMesh or UnstructuredMesh2D version)
function interpolate_data(input_data,
    mesh::Union{StructuredMesh, UnstructuredMesh2D, P4estMesh},
    n_visnodes)
    # Calculate equidistant output nodes
    nodes_out = collect(range(-1, 1, length=n_visnodes))

    return raw2interpolated(input_data, nodes_out)
end

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
