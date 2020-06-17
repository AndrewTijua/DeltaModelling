include("../src/libs/structs.jl")
#include("../src/libs/helpers.jl")

using Plots
using Random
using BandedMatrices
plotlyjs()

Random.seed!(1234)

xm = 75
ym = 75

x = 1:xm
y = 1:ym

xs = collect(Float64, x)
ys = collect(Float64, y)


h_mat = zeros(xm, ym)
for col = y, row = x
    h_mat[row, col] = 100. - 2.85*row - 2.85*col + 2.6 * rand()
end
h_mat ./= 1000.

#i hate this
# h_mat = convert(Array{Float64, 2}, sqrt.(x) * sqrt.(y)')
# h_mat ./= 16 * maximum(h_mat)
# h_mat .*= -1.
# h_mat .+= 0.02


w_mat = zeros(size(h_mat))
#w_mat[1,1] = h_mat[1,1] * 1.2
δ = 0.0025
# + 0.0005 * rand(100, 100)


w_mat = max.(w_mat, h_mat .- δ)
# nblock = 5
# for col = 1:nblock+4, row = 1:nblock+4
#     w_mat[row,col] = 1.04 * h_mat[row,col]
# end

w_mat[1,1] = h_mat[1,1] * 1.04

#h_mat -= 0.007 * BandedMatrix(Ones(xm, ym), (1,1))
plot(xs, ys, w_mat, st = :surface, c = :blues)
plot!(xs, ys, h_mat, st = :surface, c = :greens)

sbx = sby = 10

d_model = delta_model_abm_f(h_mat, w_mat, sbx, sby; timestep = 5., moore = true, c_sigma = 8.5, c1 = 0.5, c2 = 0.0)

w_lev(a) = a.w_lev
s_ele(a) = a.s_ele
w_flo(a) = a.w_flo
i_sed(a) = a.i_sed

adata = [w_lev, s_ele, w_flo, i_sed]

flow_new_delta!(d_model)

nsteps = 500

results, _ = run!(d_model, dummystep, flow_new_delta!, nsteps; adata = adata, when = [nsteps])
k = Array(results.w_lev)
m = Array(results.s_ele)
f = Array(results.w_flo)
s = Array(results.i_sed)
wl = reshape(k, size(w_mat))
se = reshape(m, size(w_mat))
wf = reshape(f, size(w_mat))
si = reshape(s, size(w_mat))
plot(xs, ys, wl, st = :surface, c = :blues)
plot!(xs, ys, se, st = :surface, c = :greens)

# plot(xs, ys, wf, st = :surface, c = :reds)
# plot!(xs, ys, si*2000., st = :surface, c = :cividis)
#
# plot(xs, ys, wl - w_mat, st = :surface, c = :blues)
# plot(xs, ys, se - h_mat, st = :surface, c = :greens)

# plot(xs, ys, wl - se .+ δ, st = :surface)

# plot(xs, ys, wf, st = :surface)

# using ProfileView
# @profview run!(d_model, dummystep, flow_step!, nsteps; adata = adata, when = [nsteps])
# @profview run!(d_model, dummystep, flow_new_delta!, 10; adata = adata, when = [10])
# @time flow_step!(d_model)
