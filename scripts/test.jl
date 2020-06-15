include("../src/libs/structs.jl")
#include("../src/libs/helpers.jl")

using Plots
using Random
using BandedMatrices
plotlyjs()

Random.seed!(1234)

x = 1:100
y = 1:100

xs = collect(Float64, x)
ys = collect(Float64, y)


h_mat = zeros(100, 100)
for col = y, row = x
    h_mat[row, col] = 100. - 0.85*row - 0.85*col + 2.6 * rand()
end
h_mat ./= 1000.

#i hate this
# h_mat = convert(Array{Float64, 2}, sqrt.(x) * sqrt.(y)')
# h_mat ./= 16 * maximum(h_mat)
# h_mat .*= -1.
# h_mat .+= 0.02


w_mat = zeros(size(h_mat))
#w_mat[1,1] = h_mat[1,1] * 1.2
δ = 0.0015
w_mat = max.(w_mat, h_mat .- δ)# + 0.0005 * rand(100, 100)
#h_mat -= 0.01 * BandedMatrix(Ones(100, 100), (15,15))
# for i in x
#     w_mat[i,i] = w_mat[i,i] * 1.2
# end

nblock = 5
for col = 10:nblock+9, row = 10:nblock+9
    w_mat[row,col] = 1.02 * h_mat[row,col]
end

plot(xs, ys, w_mat, st = :surface, c = :blues)
plot!(xs, ys, h_mat, st = :surface, c = :greens)

sbx = sby = 10

d_model = delta_model_abm_f(h_mat, w_mat, sbx, sby; timestep = 5., moore = true, c_sigma = 8.5, c1 = 0.0, c2 = 0.0)

w_lev(a) = a.w_lev
s_ele(a) = a.s_ele
w_flo(a) = a.w_flo
i_sed(a) = a.i_sed

adata = [w_lev, s_ele, w_flo, i_sed]

flow_new_delta!(d_model)

nsteps = 400

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
#plot(xs, ys, wf, st = :surface, c = :reds)

# plot(xs, ys, si, st = :surface, c = :cividis)
#
# plot(xs, ys, wl - w_mat, st = :surface, c = :blues)
# plot!(xs, ys, se - h_mat, st = :surface, c = :greens)

# plot(xs, ys, wl - se .+ δ, st = :surface)

# plot(xs, ys, wf, st = :surface)

# using ProfileView
# @profview run!(d_model, dummystep, flow_step!, nsteps; adata = adata, when = [nsteps])
# @profview run!(d_model, dummystep, flow_step!, 100; adata = adata, when = [100])
# @time flow_step!(d_model)
