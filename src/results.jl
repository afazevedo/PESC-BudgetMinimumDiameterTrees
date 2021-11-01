using CSV, JuliaDB


function results(obj, cost, tempo, path, B, L)

    array_obj = [obj]
    array_cost = [cost]
    array_tem = [tempo]
    array_budget = [B]
    array_L = [L]
    aux = split(path, '\\')
    aux2 = split(aux[length(aux)], '.')
    array_path = [aux2[1]]

    if g_params.heuristic
        heuristica = "heuristica"
    else
        heuristica = "sem-heuristica"
    end 

    model = g_params.type_of_model
    ttree = g_params.type_of_tree
    
    if g_params.warm_start
        ws = "warm-start"
    else
        ws = "sem-warm-start"
    end 
    
    if aux[length(aux)] == "c_v10_a45_d4.txt"
        t = table(array_path, array_obj, array_cost, array_budget, array_tem, array_L, names = [:Instancia, :Objetivo, :CustoÁrvore, :Budget, :Tempo, :L])
        CSV.write("results\\teste_$model-$ttree-$heuristica-$ws.csv", t, newline='\n')
    else
        t = table(array_path, array_obj, array_cost, array_budget, array_tem, array_L, names = [:Instancia, :Objetivo, :CustoÁrvore, :Budget, :Tempo, :L])
        CSV.write("results\\teste_$model-$ttree-$heuristica-$ws.csv", t, newline='\n', append = true)
    end
end
