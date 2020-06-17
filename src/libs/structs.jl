using Agents

mutable struct delta_agent <: AbstractAgent
    #=
    generic cell structure, all nodes of domain contain one (and only one) of these
    suffix _p denotes holding value for next step
    id: cell id in abm
    pos: position of cell in grid as integer tuple, also used for flow direction calculations
    w_level: water level
    s_elev: surface elevation
    in_sed: sediment carried at start of step
    is_sea_boundary: used to check if height should be held at sea level
    is_inlet_boundary: used to define inlet for water parcels and sediment parcels
    =#
    id::Int64
    pos::Tuple{Int64,Int64}
    w_lev::Float64
    w_lev_p::Float64
    s_ele::Float64
    s_ele_p::Float64
    w_flo::Float64
    w_flo_p::Float64
    i_sed::Float64
    i_sed_p::Float64
    is_sea_boundary::Bool
    is_inlet::Bool
end

function delta_model_abm_f(
    height_matrix,
    depth_matrix,
    sbx,
    sby;
    timestep = 0.1,
    eps_smooth = 1e-3,
    I_0 = 1.7e-4,
    J_0 = 2.5e-4,
    c_sigma = 8.5,
    I_ero = 4e-6,
    V_ero = 4e-6,
    theta_ero = -5e-5,
    c1 = 0.1,
    c2 = 0.1,
    moore = false,
    saturate = 0.05,
)
    properties = Dict(
        :dt => timestep,
        :ϵ => eps_smooth,
        :I_0 => I_0,
        :s_0 => J_0,
        :c_sigma => c_sigma,
        :I_ero => I_ero,
        :V_ero => V_ero,
        :ero_max => theta_ero,
        :dep_max => -theta_ero,
        :c1 => c1,
        :c2 => c2,
        :moore => moore,
        :saturate => saturate,
    )

    grid_dims = size(height_matrix)
    space = GridSpace(grid_dims, moore = moore)

    model = AgentBasedModel(delta_agent, space; properties = properties)

    _idx = 1
    for x = 1:grid_dims[1]
        for y = 1:grid_dims[2]
            if (x == 1 && y == 1)
                add_agent_pos!(
                    delta_agent(
                        _idx,
                        (x, y),
                        depth_matrix[x, y],
                        depth_matrix[x, y],
                        height_matrix[x, y],
                        height_matrix[x, y],
                        -I_0,
                        -I_0,
                        J_0,
                        J_0,
                        false,
                        true,
                    ),
                    model,
                )
            elseif ((x < grid_dims[1] && y < grid_dims[2]) || (x == grid_dims[1] && y < sby) || (y == grid_dims[2] && x < sbx))
                add_agent_pos!(
                    delta_agent(
                        _idx,
                        (x, y),
                        depth_matrix[x, y],
                        depth_matrix[x, y],
                        height_matrix[x, y],
                        height_matrix[x, y],
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        false,
                        false,
                    ),
                    model,
                )
            else
                add_agent_pos!(
                    delta_agent(
                        _idx,
                        (x, y),
                        depth_matrix[x, y],
                        depth_matrix[x, y],
                        height_matrix[x, y],
                        height_matrix[x, y],
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        true,
                        false,
                    ),
                    model,
                )
            end
            _idx += 1
        end
    end
    return model
end

function hydraulic_conductivity(node::delta_agent, neighbour::delta_agent, c_sigma::Float64)
    V_avg = 0.5 * (node.w_lev + neighbour.w_lev)
    H_avg = 0.5 * (node.s_ele + neighbour.s_ele)
    if V_avg > H_avg
        return c_sigma * (V_avg - H_avg)
    else
        return 0.0
    end
end

function pairwise_flux(node::delta_agent, neighbour::delta_agent, c_sigma::Float64)
    hc = hydraulic_conductivity(node::delta_agent, neighbour::delta_agent, c_sigma::Float64)
    dh = (node.w_lev - neighbour.w_lev)
    @assert !isnan(hc)
    @assert !isnan(dh)
    return hc * dh
end

function sum_flux(node::delta_agent, model)
    s_flux::Float64 = 0.0
    for neighbour_c in node_neighbors(node, model)
        neighbour = model[get_node_contents(neighbour_c, model)[1]]
        s_flux += pairwise_flux(node, neighbour, model.c_sigma)
    end
    return s_flux
end

function sum_flux_out(node::delta_agent, model)
    s_flux_out::Float64 = 0.0
    for neighbour_c in node_neighbors(node, model)
        neighbour = model[get_node_contents(neighbour_c, model)[1]]
        s_flux_out += max(0.0, pairwise_flux(node, neighbour, model.c_sigma))
    end
    return s_flux_out
end

function sedimentation_erosion_rate(node::delta_agent, neighbour::delta_agent, model)
    term_1::Float64 = model.c1 * (model.I_ero - abs(pairwise_flux(node, neighbour, model.c_sigma)))
    term_2::Float64 = model.c2 * (model.V_ero - abs(node.w_lev - neighbour.w_lev))
    return clamp(term_1 + term_2, model.ero_max, model.dep_max)
end


function landscape_modify!(node::delta_agent, neighbour::delta_agent, model)
    # if node.w_lev > node.s_ele
    s_e_r_raw::Float64 = sedimentation_erosion_rate(node, neighbour, model)
    @assert !isnan(s_e_r_raw)
    s_e_r::Float64 = min(node.i_sed, s_e_r_raw)
    timestep = model.dt
    #node.i_sed_p -= s_e_r
    #neighbour.i_sed_p += s_e_r
    setfield!(node, :i_sed_p, node.i_sed_p - s_e_r)
    #setfield!(neighbour, :i_sed_p, neighbour.i_sed_p + s_e_r)
    # node.s_ele_p += 0.5 * timestep * s_e_r
    # neighbour.s_ele_p += 0.5 * timestep * s_e_r
    setfield!(node, :s_ele_p, node.s_ele_p + 0.5 * timestep * s_e_r)
    setfield!(neighbour, :s_ele_p, neighbour.s_ele_p + 0.5 * timestep * s_e_r)
    # end
end

function water_depth_modify!(node::delta_agent, model)
    s_flux = sum_flux(node, model)
    timestep = model.dt
    @assert !isnan(s_flux)
    node.w_flo = node.w_flo + s_flux
    #node.w_lev_p = node.w_lev - timestep * s_flux
    @assert !isnan(node.w_flo)
    setfield!(node, :w_lev_p, node.w_lev - timestep * node.w_flo)
end

function sediment_pairwise(node::delta_agent, neighbour::delta_agent, model)
    s_flux_out = sum_flux_out(node, model)
    pw_flux = pairwise_flux(node, neighbour, model.c_sigma)
    inflowing_sediment = node.i_sed_p
    @assert !isnan(s_flux_out)
    @assert !isnan(pw_flux)
    @assert !isnan(inflowing_sediment)
    if s_flux_out > 0
        outflowing_flux_pairwise = max(0.0, pw_flux) * max(inflowing_sediment, 0.0) / s_flux_out
    else
        outflowing_flux_pairwise = 0.0
    end
    @assert !isnan(outflowing_flux_pairwise) "$s_flux_out, $pw_flux, $inflowing_sediment"
    return outflowing_flux_pairwise
end

function distribute_sediment(node::delta_agent, model)
    for neighbour_c in node_neighbors(node, model)
        neighbour = model[get_node_contents(neighbour_c, model)[1]]
        of_sed_n::Float64 = sediment_pairwise(node, neighbour, model)
        @assert !isnan(of_sed_n)
        setfield!(node, :i_sed_p, node.i_sed_p - of_sed_n)
        setfield!(neighbour, :i_sed_p, neighbour.i_sed_p + of_sed_n)
    end
end

function smooth_node!(node::delta_agent, model)
    if (node.s_ele > 0.0) && (!node.is_sea_boundary) && (!node.is_inlet) && (node.w_lev > node.s_ele)
        h_neigh::Float64 = 0.0
        n_neigh::UInt8 = 0
        w_neigh::Float64 = 0.0
        s_neigh::Float64 = 0.0
        for neighbour_c in node_neighbors(node, model)
            @inbounds neighbour = model[get_node_contents(neighbour_c, model)[1]]
            if neighbour.s_ele < neighbour.w_lev
                h_neigh += neighbour.s_ele
                w_neigh += neighbour.w_lev
                s_neigh += neighbour.i_sed
                n_neigh += 1
            end
        end
        if n_neigh > 1
            term_1::Float64 = (1.0 - model.ϵ) * node.s_ele
            term_2::Float64 = (1.0 / (n_neigh)) * model.ϵ * h_neigh
            node.s_ele = term_1 + term_2

            term_3::Float64 = (1.0 -  50.0 * model.ϵ) * node.w_lev
            term_4::Float64 = 50.0 / (n_neigh) * model.ϵ * w_neigh
            node.w_lev = term_3 + term_4

            term_5::Float64 = (1.0 - 10.0 * model.ϵ) * node.i_sed_p
            term_6::Float64 = 10.0 / (n_neigh) * model.ϵ * s_neigh
            node.i_sed_p = term_5 + term_6
        end
    end
end

function post_pro!(node::delta_agent, model)
    if (node.s_ele > 0.0) && (!node.is_sea_boundary) && (!node.is_inlet)
        node.w_lev = max(node.w_lev_p, 0.0)
    end
    node.s_ele = node.s_ele_p
    node.w_flo = node.w_flo_p
    node.i_sed = clamp(node.i_sed_p, 0.0, model.saturate)
    # node.i_sed = node.i_sed_p
    node.w_lev_p = 0.0
    # node.s_ele_p = 0.0
    node.w_flo_p = 0.0
    node.i_sed_p = 0.0
    @assert !isnan(node.s_ele)
    @assert !isnan(node.w_lev)
    @assert !isnan(node.i_sed)
    @assert !isnan(node.w_flo)
end

function flow_new_delta!(model)
    for node_c in nodes(model)
        @inbounds node = model[get_node_contents(node_c, model)[1]]

        #set boundary conditions and skip
        if node.is_sea_boundary
            node.w_lev = 0.0
            node.w_flo = 0.0
            node.i_sed = 0.0
            node.i_sed_p = 0.0
            continue
        end

        #check if inlet and if so add flow and sediment
        if node.is_inlet
            node.w_flo = -model.I_0
            node.i_sed = model.s_0
            node.i_sed_p = model.s_0
        end

        #steps 1-3
        water_depth_modify!(node, model)

        #steps 4-7
        for neighbour_c in node_neighbors(node, model)
            if node.s_ele > node.w_lev
                @inbounds neighbour = model[get_node_contents(neighbour_c, model)[1]]
                landscape_modify!(node, neighbour, model)
            end
        end

        #step 8

    end
    for node_c in nodes(model)
        @inbounds node = model[get_node_contents(node_c, model)[1]]
        distribute_sediment(node, model)
        post_pro!(node, model)
        smooth_node!(node, model)
    end
end
