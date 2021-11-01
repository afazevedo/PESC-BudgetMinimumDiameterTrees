using LightGraphs, Random, GraphPlot

include("../read.jl")
include("params.jl")
include("func.jl")

ins = read_from_files(g_params.file_name)


function OTT(ins, initial_node)
    num_nodes = collect(1:ins.n)
    reverse_num_nodes = reverse(num_nodes)

    for D in reverse_num_nodes
        flag = true
        V = 1:ins.n
        T = []
        Q = collect(V)
        chave = zeros(ins.n)
        # rnd_node = rand(V)
        push!(T, initial_node)
        filter!(e->e≠initial_node, Q)
        pt.spanTree = SimpleGraph(ins.n)
        while Q ≠ []
            c_min = 1000000000
            v = -1 
            u = -1 
            for j in Q
                chave[j] = -1 
                melhor = 1000000000000
                
                for i in T
                    # println("Distancia de i", distance(i)+1, " Diametro:", D)
                    # println("Custo:", check_cost(ins, pt.spanTree) + ins.c[i,j])
                    if distance(i) + 1 <= D && check_cost(ins, pt.spanTree) + ins.c[i,j] <= ins.B 
                        melhor = ins.c[i,j] 
                        chave[j] = i 
                    end 
                end

                k = chave[j]
                k = Int(k)

                if k ≠ -1 && ins.c[k,j] < c_min  
                    c_min = ins.c[k,j] 
                    u = k 
                    v = j 
                end    
            end

            if u == -1 || v == -1
                flag = false
                break 
            end            
            push!(T, v)
            filter!(e->e≠v, Q)
            add_edge!(pt.spanTree, u, v)
        end
        if !flag 
            return D-1
            break 
        end 
    end 
end

function distance(node)
    aux = gdistances(pt.spanTree, node)
    filter!(e->e≠9223372036854775807, aux)
    return maximum(aux)
end

function plotGraph(graph)
    nodelabel = 1:nv(graph)
    p = gplot(graph, nodelabel = nodelabel)
    display(p)
end

