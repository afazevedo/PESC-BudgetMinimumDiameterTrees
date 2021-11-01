using GraphPlot 

function MST(ins::model_params, hp::heuristic_params)
    g, cost = readGraph(ins)
    data = prim_mst(g, ins.c)
    
    for k in 1:ins.n-1
        i = src(data[k])
        j = dst(data[k])
        add_edge!(hp.spanTree, i, j)
        hp.spanCost = hp.spanCost + ins.c[i,j]
    end
    return hp.spanTree, hp.spanCost
end 

function budgetCalculator(ins::model_params, percent::Float64)
    sum = 0
    aux = copy(ins.c)
    K = vec(aux)
    sort!(K, rev=true)

    for i in 1:ins.n-1
        sum = sum + K[i]
    end
    return percent*sum
end

function readGraph(ins::model_params)
    g = SimpleGraph(ins.n)
    cost = zeros(ins.n, ins.n)

    for i in 1:ins.n 
        for j in (i+1):ins.n
            if i != j && ins.c[i,j] > g_params.eps
                add_edge!(g, i, j)
                cost[i,j] = ins.c[i,j] 
            end 
        end 
    end
    return g, cost
end

function plotGraph(graph::SimpleGraph)
    nodelabel = 1:nv(graph)
    p = gplot(graph, nodelabel = nodelabel)
    display(p)
end

function getDiameter(ins::model_params, graph::SimpleGraph)
    nodeCenter = LightGraphs.center(graph)
    
    if length(nodeCenter) < 2 
        p = bellman_ford_shortest_paths(graph, nodeCenter[1])

        diam = 0
        for i in 1:ins.n
            if p.dists[i] < 10000
                if p.dists[i] > diam
                    diam = p.dists[i]
                end
            end
        end
        return 2*diam
    else 
        p1 = bellman_ford_shortest_paths(graph, nodeCenter[1])
        p2 = bellman_ford_shortest_paths(graph, nodeCenter[2])
        diam = 0

        for i in 1:ins.n
            if p1.dists[i] < 10000 && p2.dists[i] < 10000
                if p1.dists[i] > diam
                    diam = p1.dists[i]
                elseif p2.dists[i] > diam 
                    diam = p2.dists[i]
                end 
            end
        end
        return 2*diam-1 
    end
end 

function check_cost(ins::model_params, G::SimpleGraph)
    cost = 0
    for i in 1:nv(G)
        for j in (i+1):nv(G)
            if has_edge(G,i,j)
                cost = cost + ins.c[i,j]
            end 
        end 
    end 
    return cost
end 


