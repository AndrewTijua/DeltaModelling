using StaticArrays
using LinearAlgebra
using Agents

mutable struct cell
    #=
    generic cell structure, all nodes of domain contain one (and only one) of these
    id: cell id in abm
    pos: position of cell in grid as integer tuple, also used for flow direction calculations
    w_level: water level
    s_elev: surface elevation
    visited: boolean tuple indicating visitation by
        index 1: water parcel
        index 2: sediment, sand parcel
        index 3: sediment, mud parcel
    is_sea_boundary: used to check if height should be held at sea level
    is_wall_boundary: used to check if passable to parcels
    is_inlet_boundary: used to define inlet for water parcels and sediment parcels
    =#
    id::Int64
    pos::Tuple{Int64, Int64}
    w_level::Float64
    s_elev::Float64
    # visited::Tuple{Bool, Bool, Bool}
    is_sea_boundary::Bool
    is_wall_boundary::Bool
    is_inlet_boundary::Bool
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
