using LightGraphs, Random, GraphPlot

include("../src/read.jl")
include("params.jl")
include("../src/utils.jl")
include("rgh.jl")

# g_params.file_name = pwd()*"\\instances\\$type_of_tree\\c_v10_a45_d4.txt"
# ins = read_from_files(g_params.file_name)
# ins.B = 252

function labelNode(ins::model_params, hp::heuristic_params)
    labelNodes = zeros(ins.n)
    
    # Mapear os nós através das distâncias entre ele e o centro
    if length(hp.nodeCenter) < 2
        # Nó central 
        minPath = LightGraphs.bellman_ford_shortest_paths(hp.spanTree, hp.nodeCenter[1])
        for node in 1:ins.n
            labelNodes[node] = minPath.dists[node]
        end 
    else
        # Aresta central 
        minPathLeft = LightGraphs.bellman_ford_shortest_paths(hp.spanTree, hp.nodeCenter[1])
        minPathRight = LightGraphs.bellman_ford_shortest_paths(hp.spanTree, hp.nodeCenter[2])

        for node in 1:ins.n
            if minPathLeft.dists[node] < minPathRight.dists[node]
                labelNodes[node] = minPathLeft.dists[node]
            else 
                labelNodes[node] = minPathRight.dists[node]
            end 
        end 
    end

    return labelNodes
end

function decreaseDiameter(ins::model_params, hp::heuristic_params)
    hp.nodeCenter = LightGraphs.center(hp.spanTree)
    labelNodes = labelNode(ins, hp)
    LCR = sortperm(labelNodes, rev=true)

    for i in LCR # Para cada nó na lista de candidatos
        nodeChose = -1

        #Calcular os vizinhos do nó i que possuem distância menor ou igual a 1
        vizinhos = neighbors(hp.spanTree, i)

        for j in vizinhos # Para cada vizinho do nó i
            if labelNodes[j] <= labelNodes[i]
                nodeChose = j
            end 
        end 


        minBenefit = 10000
        minInd = -1
        for k in 1:ins.n 
            if !has_edge(hp.spanTree, i, k) && nodeChose != k && labelNodes[k] < labelNodes[i] && (ins.c[k,i]-ins.c[nodeChose, i]) < minBenefit
                minBenefit = ins.c[k,i] - ins.c[nodeChose, i]
                minInd = k
            end 
        end

        minInd

        if minInd != -1 
            if hp.spanCost + (ins.c[minInd, i] - ins.c[nodeChose, i]) <= ins.B
                hp.spanCost = hp.spanCost + (ins.c[i, minInd] - ins.c[nodeChose, i])
                add_edge!(hp.spanTree, minInd, i)
                rem_edge!(hp.spanTree, nodeChose, i)
            end
        end
    end

    hp.tree_diameter = LightGraphs.diameter(hp.spanTree)
    return hp.spanTree, hp.spanCost, hp.tree_diameter 
end


