include("../src/libs/helpers.jl")

using Plots
plotlyjs()

x = 1:100
y = 1:100

xs = collect(Float64, x)
ys = collect(Float64, y)


h_mat = convert(Array{Float64, 2}, sqrt.(x) * sqrt.(y)')
h_mat ./= 16 * maximum(h_mat)
h_mat .*= -1.
h_mat .+= 0.02

w_mat = zeros(size(h_mat))
#w_mat[1,1] = h_mat[1,1] * 1.2
δ = 0.0005
w_mat = max.(w_mat, h_mat .- δ)
for i in x
    w_mat[i,i] = w_mat[i,i] * 1.2
end
plot(xs, ys, w_mat, st = :surface)
plot!(xs, ys, h_mat, st = :surface)

sbx = sby = 65

d_model = delta_model_abm(h_mat, w_mat, sbx, sby; timestep = 0.5, moore = false, c1 = 0.0)

w_lev(a) = a.w_level
s_elev(a) = a.s_elev
w_flow(a) = a.w_flow

adata = [w_lev, s_elev, w_flow]

flow_step!(d_model)

nsteps = 1000

results, _ = run!(d_model, dummystep, flow_step!, nsteps; adata = adata, when = [nsteps])
k = Array(results.w_lev)
m = Array(results.s_elev)
f = Array(results.w_flow)
wl = reshape(k, size(w_mat))
se = reshape(m, size(w_mat))
wf = reshape(f, size(w_mat))
plot(xs, ys, wl, st = :surface, c = :blues)
plot!(xs, ys, se, st = :surface, c = :greens)

plot(xs, ys, wl - se .+ δ, st = :surface)

plot(xs, ys, wf, st = :surface)

using ProfileView
# @profview run!(d_model, dummystep, flow_step!, nsteps; adata = adata, when = [nsteps])
@profview run!(d_model, dummystep, flow_step!, 100; adata = adata, when = [100])
@time flow_step!(d_model)
