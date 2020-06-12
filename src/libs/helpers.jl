# using StaticArrays
# using LinearAlgebra
#
# mutable struct cell
#     #=
#     id: cell id in abm
#     pos: position of cell in grid as integer tuple, also used for flow direction calculations
#     H: sum of depth of water and bed
#         H_old: H at prev timestep
#         H_smooth: H after diffusion applied
#         H_new: weighted sum of H_old and H_smooth
#     η: depth of bed
#     q_w: water unit discharge vector
#     sediments: volume of carried sediments:
#         index 1: sand (coarse, noncohesive)
#         index 2: mud (fine, cohesive)
#     visited: boolean tuple indicating visitation by
#         index 1: water parcel
#         index 2: sediment, sand parcel
#         index 3: sediment, mud parcel
#     =#
#     id::Int64
#     pos::Tuple{Int64, Int64}
#     H_old::Float64
#     H_smooth::Float64
#     H_new::Float64
#     η::Float64
#     q_w::MVector{2,Float64}
#     sediments::MVector{2,Float64}
#     visited::Tuple{Bool, Bool, Bool}
#     is
# end
