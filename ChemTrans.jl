### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ fc94b3fe-43d6-11eb-1fb8-2b0187a77c47
begin
	using MultivariateStats
    using JuMP
    using Cbc
    using Distances
    using Distributions
    using Random
    using Plots
	gr()
end

# ╔═╡ 603a21c6-43d4-11eb-231c-f331a607203c
md"
# Least Path Grade Transitions 
Polypropylene has many potenital grades each with unique charatersitics. Each grade is
is created by a certain reaction environment. Chgangeing this environment takes time
and causes `off-spec` or `wide-spec` product to be produced. 

To minimize this transition, it is best to move to similar reactions instead of very 
different reactions. If we look at a simple example, a planner would want to go from
low melt flow grades to high ones and then back down in a sin wave. However, there are more complex atributes like if a grade is random of homo-polymer etc. To solve this we write a simple linear program. This is a real life application of the classic [Traveling Salesperson Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem) (TSP).

I would like to thank Evan Fields and [this](https://opensourc.es/blog/mip-tsp/) blog post for getting me started.

---

## Loading Libraries
"

# ╔═╡ 55ee3af6-43d7-11eb-3975-d1105ff24f04
md"
## Generating Data
First we need to set the transition costs. In the real world, we multi dimensial
output space to calculate the transition matrix. In order to visualize this, we are going to use a 2D vector space to more easily visualize.

This difference is a bit accademic as multi factored reasoned for costs are projected
down onto a matrix."

# ╔═╡ bee109ec-43d9-11eb-054d-3bd62a0eed31
begin
    len = 10
    Random.seed!(42)
    x_low_high = rand(Uniform(0.0, 10.0), len)
    y_low_high = rand(Uniform(0.0, 10.0), len)
    

end

# ╔═╡ 5c6132b6-43da-11eb-08bc-a13fd48939e7
begin
    plot(x_low_high, y_low_high, seriestype = :scatter, title = "Distances",
		 lab = "Grade")
    plot!(x_low_high, y_low_high, lab = "Random Order Path")
end

# ╔═╡ e483ebc2-43db-11eb-0b2b-935d84aa6152
md"This is obviously not the best paht through the production wheel.

## Randomizing It
It's also not the only problem we face. In polypropylen reactions, the transition
from high to low is different from the transition from low to high even for the same
two grades. To simulate this we add some random noise."

# ╔═╡ 27324bee-43dc-11eb-2a09-ef9978862694
begin
    Random.seed!(42)


    N = length(x_low_high)
    x_perterbation = rand(Uniform(1.0, 1.5), N)
    y_perterbation = rand(Uniform(1.0, 1.5), N)

    x_high_low = x_low_high .* x_perterbation
    y_high_low = y_low_high .* y_perterbation
end

# ╔═╡ e4f39958-43dc-11eb-2d46-71a9943e514f
begin
    plot(x_low_high, y_low_high, seriestype = :scatter, title = "'Distances'", lab = "Up Costs")
    plot!(x_high_low, y_high_low, seriestype = :scatter, lab = "Down Costs")
end

# ╔═╡ 3fecca6e-43dd-11eb-0adb-5b98cc12801e
md"These perturbations exibit the increased time of transitions when going in a 
different direction.

## Product Wheel

Our product wheel should go from low to high and then back down. If our costs were the
same on the way up and down, we could just use the same optimization but we are not guaranteeded that it will be optimal. 

Solving asymetric problems like this is difficult so I am just going to break them up 
into two. If you had a traditional cost matrix that was asymetric, you need to split 
it into two symetric ones, one with the up costs on both sides of the diagonal and one with down costs. Since ours are artificial, I'm just going to keep them in their seperate forms. Just remember that most of the time we would need to split them.

So while we technically have a 'asymetric' problem with different values on the lower and upper triagnles of the cost matrix, we split them into two and solve each individually."

# ╔═╡ a5f881f6-449c-11eb-2e8f-a186c1e61c6d
function cost_generator(d1, d2)
	N = length(d1)
	cost = zeros(N, N)
	points = [[d1[i], d2[i]] for i ∈ 1:N]
	for i ∈ 1:N, j ∈ 1:N
		@inbounds cost[i, j] = euclidean(points[i], points[j])
	end
	return cost
end;

# ╔═╡ 2bfa4b84-44a1-11eb-3b10-8bd19c88bc72
cost_up = cost_generator(x_low_high, y_low_high)

# ╔═╡ d031364c-44a9-11eb-1eb5-177b2a0c4661
heatmap(cost_up[end:-1:1, :]) # make it reverse b/c plots is annoying

# ╔═╡ 9664ee8e-44a1-11eb-12c8-57a5a490ab20
cost_down = cost_generator(x_high_low, y_high_low)

# ╔═╡ 1243216a-44ac-11eb-0cbc-3b8dbec30520
heatmap(cost_down[end:-1:1, :]) # make it reverse b/c plots is annoying

# ╔═╡ c56e5bec-44ac-11eb-28f8-7fdd7289e01b
md"Classic Traveling Salesperson (TSP) problems are designed work for round trips.
While we kind of want to a round trip with different costs assoiated with the up and down portions of the product wheel, we don't want to make our constraints too complex.

As a work around I just calculate one single up trip and one down trip. We still need to change the trips to be one way. To do this we add a trip that has astronomically
high costs except for the start and end nodes which are zero. This way the end nodes are always 'next to each other' because they have the cheapest transition. The downside of this is that we have to pick the start and end nodes. To do this I picked the ones with the greatest cost (the maximum value of the cost matrix) to be the start and end.

We can then ignore this final transition when it comes to our analysis and the trip
will look like it is one way.

---

"

# ╔═╡ e5ab35d8-446b-11eb-2c46-a76a99c13983
"""
	 make_matrix_one_way(A::AbstractMatrix)
Takes a Matrix, calculates the max possible transition cost and pads it with a vector
of that cost. The highest transition cost values (farthest cities) are set to zero to allow for the MIP to complete a tour while maintaining uni-directional logic.
"""
function make_matrix_one_way(A::AbstractMatrix; idx = nothing)
    max_possible = sum(A)
	if idx == nothing
		idx = argmax(A)
	end
    id1, id2 = idx[1], idx[2]

    shortcut_vec = fill(max_possible, size(A, 1) + 1)
    shortcut_vec[[id1, id2, end]] .= 0.0
    cost = hcat(A, shortcut_vec[1:end-1])
    cost = vcat(cost, shortcut_vec')

    return cost, idx
end;


# ╔═╡ e2c3946e-446e-11eb-259d-23eeb00e9a0b
cost1, _ = make_matrix_one_way(cost_up)

# ╔═╡ 237a07d6-4471-11eb-1962-393b7e8f80ed
md"This allows us to simulate a one way trip if we remove the final index.

---

## TSP Solution
I am using the [Dantzig–Fulkerson–Johnson formulation](https://en.wikipedia.org/wiki/Travelling_salesman_problem#Dantzig%E2%80%93Fulkerson%E2%80%93Johnson_formulation)
to solve my problem. 

It is easiest understood as thining of it in two parts. First, I wrote the `dfj` 
function to initialize the model `m` the `Transition` variable, and the objective 
function.  

Secondly, the `dfj` function adds constraints to make sure that there are no zero cost 
joins back to themselves and then to make sure that each row and column have an exact 
value of 1. This makes it so each transition will happen once and only once. 

While these constraints seem complete, there is actually one aditional step that, most
of the time, needs to be taken.
"

# ╔═╡ 37710a8c-4471-11eb-0771-a3d09a1a4ddf
"
	dfj(A::AbstractArray, optimizer)
Takes cost matrix and computes the first step of the optimization.
"
function dfj(A::AbstractArray, optimizer)
    N, M = size(A)

    m = Model(optimizer)
    @assert N == M "Matrix must be square"
    @variable(m, Transition[1:N, 1:M], Bin)

    @objective(m, Min, sum(Transition .* A))

    trans_iterator = 1:N
    for i ∈ trans_iterator
        @constraint(m, Transition[i, i] == 0)
        # ^ No self linking
        @constraint(m, sum(Transition[i, :]) == 1)
        # ^ All must be used
    end

    for j ∈ trans_iterator
        @constraint(m, sum(Transition[:, j]) == 1)
    end

    for i ∈ trans_iterator, j ∈ trans_iterator
        @constraint(m, Transition[i, j] + Transition[j, i] ≤ 1)
        # ^ makes it impossible to have to "same" tranistions
        # on the diagonal or thought of another way, cities self joining with
        # eachother
    end

    optimize!(m)
    return m
end;

# ╔═╡ 03dc832a-44af-11eb-3b0d-bfca4a6a6e0f
md"

Unfortunetly, the constraints in the `djf` function  cause it to often return sub tours instead of a tour that hits every path. Think:

A -> B -> C -> A  and D -> E -> D

This satisfies mini TSP problems but not the global one. We could exhaustivly limit
these issues but MIPs get very difficult to solve if we have too many constraints.

A way around this is to [Delayed Column Generation](https://en.wikipedia.org/wiki/Column_generation)
We run the model, check if it has tours that are too small, limit those tours with a
new constraint, and then run it again. We repeate this until we get a satisfactory result.

This provides the optimal solution with mimimal constraints.
"

# ╔═╡ e60fa9f2-4473-11eb-311a-37be5c9b725f
"	delayed_column_gen(m)
Takes a JuMP model (genereated from the dfj() and adds constraints as necessary to add complete tours.
"
function delayed_column_gen(m)
    var = m.obj_dict[:Transition]
    vals = value.(var)

    N = size(vals, 1)
    Q = [1]
    while true
        idx = argmax(vals[Q[end], :])
        if idx == Q[1] # checks if tour is completed and short
            break
        else
            push!(Q, idx)
        end
    end

    if length(Q) < N
        @show Q
        @constraint(m, sum(var[Q, Q]) ≤ length(Q) - 1) # 'Bans' the current solution
        optimize!(m)
        delayed_column_gen(m)
    else
        return Q
    end

end;

# ╔═╡ 960452c2-447e-11eb-1ff3-2f96a4007940
m = dfj(cost1, Cbc.Optimizer);

# ╔═╡ db0945f0-4481-11eb-3115-03642055cd17
begin
    best_order_with_skip = delayed_column_gen(m)
    skip_idx = argmax(best_order_with_skip)
    best_order =
        vcat(best_order_with_skip[1+skip_idx:end], best_order_with_skip[1:skip_idx-1])
end

# ╔═╡ 35425012-448b-11eb-21c9-eb96657103e5
begin
    plot(
        x_low_high,
        y_low_high,
        seriestype = :scatter,
        title = "Distances",
        series_annotations = [string(i) for i ∈ 1:10],
		lab = "Grades"    
    )
    plot!(x_low_high[best_order], y_low_high[best_order], lab= "Path")
end

# ╔═╡ 599eec44-4559-11eb-16f6-4d1f19d32b45
md"Next, lets put these into convenient functions to find the up and down tours."

# ╔═╡ 0aa28e0a-448c-11eb-0b26-6900659e5015
begin
	function find_path(costs, optimizer)
		costs, idx = make_matrix_one_way(costs)
		m = dfj(costs, optimizer)
		
		best_order = delayed_column_gen(m)
	    skip_idx = argmax(best_order)
	    if skip_idx == length(best_order)
			best_order = best_order[1:end-1]
		elseif skip_idx == 1
			best_order = best_order[2:end]
		else
	    	best_order = vcat(best_order_with_skip[1+skip_idx:end], 
						 	 best_order_with_skip[1:skip_idx-1])
		end
		return best_order, idx
		
	end
	function find_path(costs, optimizer, start_end)
		
		costs, idx = make_matrix_one_way(costs, idx=start_end)
		m = dfj(costs, optimizer)
		
		best_order = delayed_column_gen(m)
		# return best_order
	    skip_idx = argmax(best_order)
		if skip_idx == length(best_order)
			best_order = best_order[1:end-1]
		elseif skip_idx == 1
			best_order = best_order[2:end]
		else
	    	best_order = vcat(best_order_with_skip[1+skip_idx:end], 
						 	 best_order_with_skip[1:skip_idx-1])
		end
		return best_order
	end
end;

# ╔═╡ 2173dfb0-449a-11eb-3a27-937c4337583f
up_order, idx = find_path(cost_up, Cbc.Optimizer)

# ╔═╡ 834be1f8-4559-11eb-20c2-b749e60f142f
md"Note in the bellow function call we are using `idx` which substitutes the maximum 
cost of the down matrix with  the start and end nodes of the previous up path. 

If instead we wanted to arbitrarily set the start and end points we could supply them directly for both the up and down paths."

# ╔═╡ c454902a-449b-11eb-1c7f-27307898dd2d
down_order = find_path(cost_down, Cbc.Optimizer, idx)

# ╔═╡ 7ab6424c-44a5-11eb-0b47-57a7bcd85ae8
begin
    plot(
        x_low_high,
        y_low_high,
        seriestype = :scatter,
        title = "Distances",
        series_annotations = [string(i) * "_up" for i ∈ 1:10],
		lab = "Up Direction"
    )
    plot!(x_low_high[up_order], y_low_high[up_order], lab = "Up Path")
	
    plot!(
        x_high_low,
        y_high_low,
        seriestype = :scatter,
        title = "Distances",
        series_annotations = [string(i) * "_dwn" for i ∈ 1:10],
		lab = "Down Direction"
    )
    plot!(x_high_low[down_order], y_high_low[down_order], lab = "Down Path")

end

# ╔═╡ 4eee53fa-44a7-11eb-2d86-1d917f5c706c
(up_order, down_order)

# ╔═╡ 97ba8db8-455a-11eb-1f01-c1180ca69ccb
md"
With a seed of 42, we get the same values for up and down. This will often be the case but they can potentially be different.

# Final Steps
Often we will just have a cost matrix with no corresponding x, y vectors that we can plot the path on. (Remember above we calculated the costs with the Euclideian function.)

We use a technique called [Classical Multidimensional Scaling (MDS)](https://multivariatestatsjl.readthedocs.io/en/stable/cmds.html)
to project the distance matrix into two vectors. 

The orientation and scaling are different but reletive positions are maintained. 

This captures the non-linear transform of the Euclidean formula very nicely. In real world grade transitions it may not perform exactly as well but it should at least maintain rough differences.
"

# ╔═╡ 5f872722-4539-11eb-3539-a9c58ad4f2b4
est = classical_mds(cost_up, 2)

# ╔═╡ 80399e6e-4539-11eb-2ee1-01ec5c22326a
begin
	d1, d2 = (2, 1)
	l = @layout [a b]
	p1 = plot(est[d1, :], est[d2, :], 
			  seriestype = :scatter, 
			  series_annotations = [string(i) for i ∈ 1:10])
	plot!(est[d1, :][up_order], est[d2, :][up_order])
	p2 = plot(x_low_high, y_low_high,
				seriestype = :scatter, 
		series_annotations = [string(i) for i ∈ 1:10])
	plot!(x_low_high[up_order], y_low_high[up_order])
	
	plot(p1, p2, layout = l, legend = false)
end


# ╔═╡ Cell order:
# ╟─603a21c6-43d4-11eb-231c-f331a607203c
# ╠═fc94b3fe-43d6-11eb-1fb8-2b0187a77c47
# ╟─55ee3af6-43d7-11eb-3975-d1105ff24f04
# ╠═bee109ec-43d9-11eb-054d-3bd62a0eed31
# ╠═5c6132b6-43da-11eb-08bc-a13fd48939e7
# ╟─e483ebc2-43db-11eb-0b2b-935d84aa6152
# ╠═27324bee-43dc-11eb-2a09-ef9978862694
# ╠═e4f39958-43dc-11eb-2d46-71a9943e514f
# ╟─3fecca6e-43dd-11eb-0adb-5b98cc12801e
# ╠═a5f881f6-449c-11eb-2e8f-a186c1e61c6d
# ╠═2bfa4b84-44a1-11eb-3b10-8bd19c88bc72
# ╠═d031364c-44a9-11eb-1eb5-177b2a0c4661
# ╠═9664ee8e-44a1-11eb-12c8-57a5a490ab20
# ╠═1243216a-44ac-11eb-0cbc-3b8dbec30520
# ╟─c56e5bec-44ac-11eb-28f8-7fdd7289e01b
# ╠═e5ab35d8-446b-11eb-2c46-a76a99c13983
# ╠═e2c3946e-446e-11eb-259d-23eeb00e9a0b
# ╟─237a07d6-4471-11eb-1962-393b7e8f80ed
# ╠═37710a8c-4471-11eb-0771-a3d09a1a4ddf
# ╟─03dc832a-44af-11eb-3b0d-bfca4a6a6e0f
# ╠═e60fa9f2-4473-11eb-311a-37be5c9b725f
# ╠═960452c2-447e-11eb-1ff3-2f96a4007940
# ╠═db0945f0-4481-11eb-3115-03642055cd17
# ╠═35425012-448b-11eb-21c9-eb96657103e5
# ╟─599eec44-4559-11eb-16f6-4d1f19d32b45
# ╠═0aa28e0a-448c-11eb-0b26-6900659e5015
# ╠═2173dfb0-449a-11eb-3a27-937c4337583f
# ╟─834be1f8-4559-11eb-20c2-b749e60f142f
# ╠═c454902a-449b-11eb-1c7f-27307898dd2d
# ╠═7ab6424c-44a5-11eb-0b47-57a7bcd85ae8
# ╠═4eee53fa-44a7-11eb-2d86-1d917f5c706c
# ╟─97ba8db8-455a-11eb-1f01-c1180ca69ccb
# ╠═5f872722-4539-11eb-3539-a9c58ad4f2b4
# ╠═80399e6e-4539-11eb-2ee1-01ec5c22326a
