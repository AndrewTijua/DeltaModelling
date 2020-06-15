using StaticArrays
using LinearAlgebra
using Agents

#there is no god

mutable struct cell <: AbstractAgent
    #=
    generic cell structure, all nodes of domain contain one (and only one) of these
    id: cell id in abm
    pos: position of cell in grid as integer tuple, also used for flow direction calculations
    w_level: water level
    s_elev: surface elevation
    in_sed: sediment carried at start of step
    in_sed_new: sediment carried at end of step
    is_sea_boundary: used to check if height should be held at sea level
    is_inlet_boundary: used to define inlet for water parcels and sediment parcels
    =#
    id::Int64
    pos::Tuple{Int64,Int64}
    w_level::Float64
    s_elev::Float64
    w_flow::Float64
    in_sed::Float64
    in_sed_new::Float64
    # visited::Tuple{Bool, Bool, Bool}
    is_sea_boundary::Bool
    is_inlet::Bool
end

# mutable struct water_parcel
#     id::Int64
#     pos::Tuple{Int64, Int64}
#     volume::Float64
# end
#
# #=
# type:
#     fine: sand, noncohesive, carried as bedload
#     coarse: mud, cohesive, carried as suspended load
# =#
# mutable struct fine_sed_parcel
#     id::Int64
#     pos::Tuple{Int64, Int64}
#     volume::Float64
# end
#
# mutable struct coarse_sed_parcel
#     id::Int64
#     pos::Tuple{Int64, Int64}
#     volume::Float64
# end

function delta_model_abm(
    height_matrix,
    depth_matrix,
    sbx,
    sby;
    timestep = 0.1,
    eps_smooth = 1e-4,
    I_0 = 1.7e-4,
    J_0 = 2.5e-4,
    c_sigma = 8.5,
    I_ero = 4e-6,
    V_ero = 4e-6,
    theta_ero = -5e-7,
    c1 = 0.1,
    c2 = 0.1,
    moore = false,
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
    )

    grid_dims = size(height_matrix)
    space = GridSpace(grid_dims, moore = moore)

    model = AgentBasedModel(cell, space; properties = properties)

    _idx = 1
    for x = 1:grid_dims[1]
        for y = 1:grid_dims[2]
            if (x == 1 && y == 1)
                add_agent_pos!(cell(_idx, (x, y), depth_matrix[x, y], height_matrix[x, y], I_0, J_0, 0.0, false, true), model)
            elseif ((x < grid_dims[1] && y < grid_dims[2]) || (x == grid_dims[1] && y < sby) || (y == grid_dims[2] && x < sbx))
                add_agent_pos!(cell(_idx, (x, y), depth_matrix[x, y], height_matrix[x, y], 0.0, 0.0, 0.0, false, false), model)
            else
                add_agent_pos!(cell(_idx, (x, y), depth_matrix[x, y], height_matrix[x, y], 0.0, 0.0, 0.0, true, false), model)
            end
            _idx += 1
        end
    end
    return model
end

function h_conduct(Vi::Float64, Vj::Float64, Hi::Float64, Hj::Float64, c_sigma::Float64)
    V_avg = 0.5 * (Vi + Vj)
    H_avg = 0.5 * (Hi + Hj)
    if V_avg > H_avg
        return c_sigma * (V_avg - H_avg)
    else
        return 0.0
    end
end

function sed_ero_rate(c1::Float64, c2::Float64, I_ero::Float64, V_ero::Float64, I_ij::Float64, Vi::Float64, Vj::Float64, ero_max::Float64)
    dsij = c1 * (I_ero - abs(I_ij)) + c2 * (V_ero - abs(Vi - Vj))
    if dsij > ero_max
        return max(dsij, -ero_max)
    else
        return ero_max
    end
end

function J_ij_out_help(J_ik_in::Float64, I_out::Float64, I_ij_out::Float64)
    if I_out == 0.0
        return 0.0
    else
        return (J_ik_in * I_ij_out / I_out)
    end
end

function setup_wf!(model)
    for node_c in nodes(model)
        node = model[get_node_contents(node_c, model)[1]]
        node.w_flow = 0.0
        in_sediment = node.in_sed
        node.in_sed_new = in_sediment
    end
end

function flow_sediment_balance!(model)
    I_ij::Float64 = 0.0
    sigma_ij::Float64 = 0.0
    J_ij_out::Float64 = 0.0
    dS_ij::Float64 = 0.0
    out_water::Float64 = 0.0
    for node_c in nodes(model)
        node = model[get_node_contents(node_c, model)[1]]

        #set boundary conditions and skip
        if node.is_sea_boundary
            node.w_level = 0.0
            node.w_flow = 0.0
            node.in_sed = 0.0
            node.in_sed_new = 0.0
        end

        #check if inlet and if so add flow and sediment
        if node.is_inlet
            node.w_flow = model.I_0
            node.in_sed = model.s_0
            node.in_sed_new = model.s_0
        end

        #calculate sum(abs(I_ik^out))
        # out_water = 0.0
        # for neighbor_c in node_neighbors(node_c, model)
        #     neighbor = model[get_node_contents(neighbor_c, model)[1]]
        #     I_ij = (node.w_level - neighbor.w_level) * h_conduct(node.w_level, neighbor.w_level, node.s_elev, neighbor.s_elev, model.c_sigma)
        #     if I_ij <= 0
        #         out_water += abs(I_ij)
        #     end
        # end


        #calculate erosive flow and distribute
        for neighbor_c in node_neighbors(node_c, model)
            neighbor = model[get_node_contents(neighbor_c, model)[1]]

            sigma_ij = h_conduct(node.w_level, neighbor.w_level, node.s_elev, neighbor.s_elev, model.c_sigma)
            I_ij = sigma_ij * (node.w_level - neighbor.w_level)

            node.w_flow = node.w_flow + I_ij
            # setfield!(node, :w_flow, node.w_flow + I_ij)

            dS_ij = sed_ero_rate(model.c1, model.c2, model.I_ero, model.V_ero, I_ij, node.w_level, neighbor.w_level, model.ero_max)
            setfield!(node, :s_elev, node.s_elev + 0.5 * model.dt * dS_ij)
            setfield!(neighbor, :s_elev, neighbor.s_elev + 0.5 * model.dt * dS_ij)
            #nH_i = node.s_elev + 0.5 * model.dt * dS_ij
            #nH_j = neighbor.s_elev + 0.5 * model.dt * dS_ij

            #J_ij_out = J_ij_out_help(in_sediment, out_water, I_ij)
            # node.in_sed -= J_ij_out
            # neighbor.in_sed_new += J_ij_out
            setfield!(node, :in_sed_new, node.in_sed_new - dS_ij)
            setfield!(neighbor, :in_sed_new, neighbor.in_sed + dS_ij)

            # node.s_elev = nH_i
            # neighbor.s_elev = nH_j
            #setfield!(node, :s_elev, node.s_elev + 0.5 * model.dt * dS_ij)
            #setfield!(neighbor, :s_elev, neighbor.s_elev + 0.5 * model.dt * dS_ij)
        end

    end
end

function propagate_sediment!(model)
    for node_c in nodes(model)
        node = model[get_node_contents(node_c, model)[1]]
        #change water level for non inlet nodes
        insed = node.in_sed_new
        node.in_sed = insed
        node.in_sed_new = 0.0
        if !((node.is_inlet) || (node.is_sea_boundary))
            if node.s_elev > 0
                node.w_level = node.w_level - node.w_flow * model.dt
            end
        end
    end
end

function finalise_smooth!(model)
    h_neigh::Float64 = 0.0
    num_neigh::UInt8 = 0
    for node_c in nodes(model)
        node = model[get_node_contents(node_c, model)[1]]
        h_neigh = 0.0
        w_neigh = 0.0
        num_neigh = 0
        for neighbor_c in node_neighbors(node_c, model)
            neighbor = model[get_node_contents(neighbor_c, model)[1]]
            num_neigh += 1
            h_neigh += neighbor.s_elev
            w_neigh += neighbor.w_level
        end
        if node.w_level > node.s_elev
            node.s_elev = (1. - model.ϵ) * node.s_elev + model.ϵ / num_neigh * h_neigh
        end
    end
end

function flow_step!(model)
    setup_wf!(model)
    flow_sediment_balance!(model)
    propagate_sediment!(model)
    finalise_smooth!(model)
end
