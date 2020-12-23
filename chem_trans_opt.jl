using JuMP
using Cbc
using Distances
using Distributions
using Random
using Plots
using TravelingSalesmanExact

## Data Gen
N = 10

Random.seed!(42)
x_low_high = rand(Uniform(0.0, 10.0), N) #  this is our cost 
y_low_high = rand(Uniform(0.0, 10.0), N)
x_low_high = [0.1, 0.2, 0.12, 2.1, 3.2, 6.0, 7.4, 8.1, 9.9, 8.5]
y_low_high = [10.3, 7.9, 6.3, 4.9, 2.0, 1.1, 1.2, 0.9, 0.7, 0.1]
plot(x_low_high, y_low_high, seriestype = :scatter, title = "'Distances'")

x_perterbation = rand(Uniform(1.0, 1.5), N)
y_perterbation = rand(Uniform(1.0, 1.5), N)

x_high_low = x_low_high .* x_perterbation
y_high_low = y_low_high .* y_perterbation

plot!(x_high_low, y_high_low, seriestype = :scatter)
# Making the matrix
# it's easeier thinking of these as cities in 2D so we wil lgenerate the matrix 
# of costs from these starting points

up_points = [[x_low_high[i], y_low_high[i]] for i ∈ 1:N]
down_points = [[x_high_low[i], y_high_low[i]] for i ∈ 1:N]

cost_up = zeros(N, N)
for i ∈ 1:N, j ∈ 1:N
    @inbounds cost_up[i, j] = euclidean(up_points[i], up_points[j])
end

cost_down = zeros(N, N)
for i ∈ 1:N, j ∈ 1:N
    @inbounds cost_down[i, j] = euclidean(down_points[i], down_points[j])
end

cost_down

costs = similar(cost_up)

for i ∈ 1:N, j ∈ 1:N
    if i == j
        costs[i, j] = 0.0
    elseif i < j
        costs[i, j] = cost_up[i, j]
    else
        costs[i, j] = cost_down[i, j]
    end
end

# so we have related but different costs on the way up and down.
# We need find the best up path and down path, then add those together


## Make trip 1 way
max_up_cost = sum(cost_up)

_, idx = findmax(cost_up) # Most expensive transition in single
id1, id2 = idx[1], idx[2]

oneway = zeros(N + 1)
oneway .= max_up_cost
oneway[[id1, id2, end]] .= 0.0
oneway

cost_up = hcat(cost_up, oneway[1:end-1])
cost_up = vcat(cost_up, oneway')

## model
m = Model(Cbc.Optimizer)

@variable(m, UpGrades[1:N+1, 1:N+1], Bin)
@objective(m, Min, sum(UpGrades .* cost_up))

grade_iterator = 1:N+1

for i ∈ grade_iterator
    @constraint(m, UpGrades[i, i] == 0) # Noself linking
    @constraint(m, sum(UpGrades[i, :]) == 1)
end

for j ∈ 1:N
    @constraint(m, sum(UpGrades[:, j]) == 1)
end

for f ∈ 1:N, t ∈ 1:N
    @constraint(m, UpGrades[f, t] + UpGrades[t, f] ≤ 1)
end

optimize!(m)

function is_tsp_solved(m, x)
    N = size(x)[1]
    x_val = JuMP.value.(x)

    # find cycle
    cycle_idx = Int[]
    push!(cycle_idx, 1)
    while true
        v, idx = findmax(x_val[cycle_idx[end], 1:N])
        if idx == cycle_idx[1]
            break
        else
            push!(cycle_idx, idx)
        end
    end
    println("cycle_idx: ", cycle_idx)
    println("Length: ", length(cycle_idx))
    if length(cycle_idx) < N
        @constraint(m, sum(x[cycle_idx, cycle_idx]) <= length(cycle_idx) - 1)
        return false
    end
    return true
end

while !is_tsp_solved(m, UpGrades)
    optimize!(m)
end

## Eval
tran_idxs = Int.(value.(UpGrades))
(cost_up .+ 0.01) .* tran_idxs
## TSP
get_optimal_tour(cost_up, Cbc.Optimizer, verbose = true, symmetric = true)


## Plotting
plot(x_low_high, y_low_high, seriestype = :scatter, title = "'Distances'")

plot_mask =[1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
plot!(
    x_low_high[plot_mask],
    y_low_high[plot_mask],
)