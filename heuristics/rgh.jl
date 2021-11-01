using LightGraphs, Random


# mutable struct heuristic_params{F<:Float64, S<:String, I<:Int64}
#     spanTree::SimpleGraph{I}
#     spanCost::F
#     tree_diameter::I
#     which_heuristic::S
# end

include("../src/read.jl")
include("params.jl")
include("../src/utils.jl")


ins = read_from_files(g_params.file_name)
ins.B = 252


function rgh(ins::model_params, hp::heuristic_params)
    hp.spanTree = readGraph(ins)
    num_nodes = collect(1:ins.n)
    reverse_num_nodes = reverse(2:ins.n)
    best_diameter_viable = MAX_ITER
    best_tree = SimpleGraph{Int64}

    for D in reverse_num_nodes
        flag = true
        T = []
        Q = collect(1:ins.n)
        p = fill(MAX_ITER, ins.n)
        num_nodes = collect(1:ins.n)
        optimal_tree = SimpleGraph(ins.n)

        if D % 2 == 0
            v = rand(num_nodes)
            p[v] = 0
            push!(T, v)
            filter!(e->e≠v, Q)
        else 
            v = rand(num_nodes)
            p[v] = 0
            push!(T, v)
            filter!(e->e≠v, Q)
            node_neighbors = neighbors(hp.spanTree, v)
            w = rand(node_neighbors)
            p[w] = 0
            push!(T, w)
            filter!(e->e≠w, Q)
            add_edge!(optimal_tree, v, w)
        end

        while Q ≠ []
            v = rand(Q)
            melhor = MAX_ITER
            u = -1
            
            for i in T 
                if has_edge(hp.spanTree, i, v) && ins.c[i,v] < melhor && (p[i]+1) <= D/2
                    u = i
                    melhor = ins.c[i,v]
                end 
            end 

            if u == -1 
                flag = false
                break 
            end       
            
            p[v] = p[u] + 1
            push!(T, v)
            filter!(e->e≠v, Q)
            add_edge!(optimal_tree, u, v)
        end

        if flag && LightGraphs.diameter(optimal_tree) < best_diameter_viable && check_cost(ins, optimal_tree) <= ins.B        
            best_diameter_viable = LightGraphs.diameter(optimal_tree)
            best_tree = copy(optimal_tree)
        end
    end 
    return best_diameter_viable
end


