using JuMP
using Cbc
using Distances
using Distributions 
using Random
using Plots

## Data Gen
n_grades = 10

Random.seed!(42)
x_low_high = rand(Uniform(0.0, 10.0), n_grades) #  this is our cost 
y_low_high = rand(Uniform(0.0, 10.0), n_grades)

plot(x_low_high, y_low_high, seriestype = :scatter, title = "'Distances'")

x_perterbation = rand(Uniform(1.0,1.5),n_grades)
y_perterbation = rand(Uniform(1.0, 1.5), n_grades) 

x_high_low = x_low_high .* x_perterbation
y_high_low = y_low_high .* y_perterbation

plot!(x_high_low, y_high_low, seriestype = :scatter)
# Making the matrix
# it's easeier thinking of these as cities in 2D so we wil lgenerate the matrix 
# of costs from these starting points

up_points = [[x_low_high[i], y_low_high[i]] for i ∈ 1:n_grades]
down_points = [[x_high_low[i], y_high_low[i]] for i ∈ 1:n_grades]

cost_up = zeros(n_grades, n_grades)
for i ∈ 1:n_grades, j ∈ 1:n_grades
    @inbounds cost_up[i, j] = euclidean(up_points[i], up_points[j])
end

cost_down = zeros(n_grades, n_grades)
for i ∈ 1:n_grades, j ∈ 1:n_grades
    @inbounds cost_down[i, j] = euclidean(down_points[i], down_points[j])
end

cost_down

costs = similar(cost_up)

for i ∈ 1:n_grades, j ∈ 1:n_grades
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

m = Model(Cbc.Optimizer)

@variable(m, UpGrades[1:n_grades, 1:n_grades], Bin)

for i in 1:n_grades
    
end