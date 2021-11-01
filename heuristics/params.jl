using LightGraphs

ins = read_from_files(g_params.file_name)
mutable struct heuristic_params{F<:Float64, S<:String, I<:Int64}
    spanTree::SimpleGraph{I}
    spanCost::F
    tree_diameter::I
    which_heuristic::S
end

grafo_inicial = SimpleGraph(ins.n)
custo_inicial = 0.0
which_heuristic = "ott"

hp = heuristic_params{Float64, String, Bool, Int64}(grafo_inicial, custo_inicial, which_heuristic)