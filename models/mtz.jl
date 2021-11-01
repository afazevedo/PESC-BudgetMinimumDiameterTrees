include("../params.jl")
include("../read.jl")
include("../heuristics/func.jl")
include("../heuristics/params.jl")
include("../heuristics/ott.jl")

using Gurobi, JuMP, MathOptInterface, LightGraphs, LightGraphsFlows

# g_params.file_name = pwd()*"\\instances\\$type_of_tree\\c_v10_a45_d4.txt"
ins = read_from_files(g_params.file_name)
ins.B = budgetCalculator(ins, 0.30)
min_d = 1000000
for i in 1:ins.n
    d = OTT(ins, i)
    if d < min_d 
        min_d = d
    end 
end 


function isodd(x)
    if x % 2 == 1
        return true
    end
    return false
end

function my_callback_function(cb_data)
    # Inicializar a relaxação linear
    linear_rel = zeros(ins.n+1, ins.n)

    # Pegar a relaxação linear de cada par (i,j)
    for i in 1:ins.n+1
        for j in 1:ins.n
            if i != j
                linear_rel[i,j] = callback_value(cb_data, x[i,j])
            end
        end
    end

    #Criação do grafo suporte direcionado
    flow_graph = LightGraphs.DiGraph(ins.n+1)
    capacity_matrix = fill(10000.0, ins.n+1, ins.n+1)

    for j in 1:ins.n
        for i in 1:ins.n+1
            if linear_rel[i,j] > g_params.eps && i != ins.n+1
                capacity_matrix[i,j] = linear_rel[i,j]
                add_edge!(flow_graph, i, j)
            end
        end
    end

    #Criação do grafo suporte direcionado ao contrário
    flow_graph_inv = LightGraphs.reverse(flow_graph)
    capacity_matrix_inv = transpose(capacity_matrix)

    list_of_vertices = collect(1:ins.n)
    shuffle!(list_of_vertices)

    number_of_cuts = 0
    # Para cada nó t dos nós terminais:
    for t in list_of_vertices
        
        #Rodando o algoritmo de min-cut saindo do vértice artificial para cada vértice terminal
        subtour, complementary, flow_value = LightGraphsFlows.mincut(flow_graph, ins.n+1, t, capacity_matrix, DinicAlgorithm())

        # Verificar desigualdade de cutset
        if flow_value < 1 + g_params.eps
            con = @build_constraint(sum(x[i,j] for i in subtour, j in complementary if i != j) >= 1)
            MOI.submit(model, MOI.UserCut(cb_data), con)
            number_of_cuts += 1
        end

        #Rodando o algoritmo de min-cut saindo de cada vértice terminal para o vértice artifical no grafo invertido
        subtour, complementary, flow_value = LightGraphsFlows.mincut(flow_graph_inv, t, ins.n+1, capacity_matrix_inv, DinicAlgorithm())

        # Verificar desigualdade de cutset
        if flow_value < 1 + g_params.eps
            con = @build_constraint(sum(x[i,j] for i in complementary, j in subtour if i != j ) >= 1)
            MOI.submit(model, MOI.UserCut(cb_data), con)
            number_of_cuts += 1
        end
    
        if number_of_cuts >= g_params.max_cuts
            break
        end
    end
end

# function solve_spanning_mtz(ins::model_params, spanTree)


if isodd(min_d)
    ins.L = (min_d+1)/2
else 
    ins.L = min_d/2
end 

# Inicialização de um modelo 
model = Model(Gurobi.Optimizer)

# Atributos do modelo
MOI.set(model, MOI.RawParameter("PreCrush"), 1)
set_optimizer_attribute(model, "OutputFlag", 0)
set_optimizer_attribute(model, "Threads", 1)
set_optimizer_attribute(model, "TimeLimit", g_params.time_limit)
MOI.set(model, MOI.RawParameter("Cuts"), 0) # Desabilitar cortes
MOI.set(model, MOI.RawParameter("Presolve"), 0) # Desabilitar presolve
MOI.set(model, MOI.RawParameter("Heuristics"), 0) # Desabilitar heurísticas

#Conjuntos
V = 1:ins.n
V0 = 1:ins.n+1

# Definição das variáveis de decisão
@variables(model, begin
    x[i in 1:ins.n+1, j in 1:ins.n; i != j], (base_name = "x[$i,$j]"), Bin
    0 <= u[i in V0] <= ins.L+1, (base_name = "u[$i]"), Int
    z, Bin
    0 <= w[i in 1:ins.n, j in 1:ins.n] <= 1, Int
    s >= 1
end)

# if g_params.warm_start
#     pt.nodeCenter = LightGraphs.center(spanTree)
#     JuMP.set_start_value(s, ins.L+1)
#     JuMP.set_start_value(u[ins.n+1], 0.0)

#     if length(pt.nodeCenter) < 2
#         JuMP.set_start_value(x[ins.n+1, pt.nodeCenter[1]], 1.0)
#         JuMP.set_start_value(z, 0.0)
#     else 
#         JuMP.set_start_value(z, 1.0)
#         JuMP.set_start_value(x[ins.n+1, pt.nodeCenter[1]], 1.0)
#         JuMP.set_start_value(x[ins.n+1, pt.nodeCenter[2]], 1.0)
#         JuMP.set_start_value(w[pt.nodeCenter[1], pt.nodeCenter[2]], 1.0)
#     end 

#     pt.spanTree = spanTree
#     labelNode(ins, pt)

#     for i in 1:nv(spanTree)
#         JuMP.set_start_value(u[i], pt.labelNodes[i]+1)
#         for j in 1:nv(spanTree)
#             if i != j
#                 if has_edge(spanTree, i, j) && pt.labelNodes[i] < pt.labelNodes[j]
#                     JuMP.set_start_value(x[i,j], 1.0)
#                 else
#                     JuMP.set_start_value(x[i,j], 0.0) 
#                 end
#             end 
#             if length(pt.nodeCenter) < 2
#                 JuMP.set_start_value(w[i,j], 0.0)
#             end 
#         end 
#     end 
# end 

# Definição da função objetivo
@objective(model, Min, 2*(s-1) + z)

# Definição das restrições
@constraints(model, begin
        root, sum(x[ins.n+1,j] for j in V) == z+1
        in_deg[j in V], sum(x[i,j] for i in V0 if i != j) == 1
    end)

@constraints(model, begin
        mtz[j in V], u[ins.n+1] - u[j] + (ins.L+1)*x[ins.n+1, j] <= ins.L
        mtz1[i in V, j in V; i != j], u[i] - u[j] + (ins.L+1)*x[i,j] + (ins.L-1)*x[j,i] <= ins.L
        # mtz2[i in V], u[i] - u[ins.n+1] + (ins.L-1)*x[ins.n+1,i] <= ins.L
    end)

@constraints(model, begin
        budget, sum(ins.c[i,j]*x[i,j] for i in V for j in V if i != j) + sum(ins.c[i,j]*w[i,j] for i in V for j in V if i != j) <= ins.B
        art[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, i]
        art2[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, j]
        central, sum(w[i,j] for i in V, j in V if i != j) == z
        mag[i in V], s >= u[i]
    end)   

for i in 1:ins.n
    for j in 1:ins.n
        if ins.c[i,j] < g_params.eps && i != j
            set_upper_bound(x[i,j], 0.0)
            set_lower_bound(x[i,j], 0.0)
            set_upper_bound(w[i,j], 0.0)
            set_lower_bound(w[i,j], 0.0)
        end
    end
end

set_upper_bound(u[ins.n+1], 0.0)
set_lower_bound(u[ins.n+1], 0.0)

# MOI.set(model, MOI.UserCutCallback(), my_callback_function)
optimize!(model)

# Parâmetros do modelo
status = termination_status(model)
tempo = solve_time(model)
nodes = node_count(model)
obj = objective_bound(model)

# println("Valor da função objetivo: ", objective_value(model))
# println("O tempo foi: ", tempo)
# println("O valor de B foi: ", ins.B)
# println("O valor de L foi: ", ins.L)

# return obj, sum(ins.c[i,j]*value(x[i,j]) for i in V for j in V if i != j) + sum(ins.c[i,j]*value(w[i,j]) for i in V for j in V if i != j), tempo
# # Relaxação Linear
# relax_integrality(model)
# optimize!(model)
# relax = objective_bound(model)
# end 

function solve_steiner_mtz(ins::model_params)
    function my_callback_function(cb_data)
        # Inicializar a relaxação linear
        linear_rel = zeros(ins.n+1, ins.n)
    
        # Pegar a relaxação linear de cada par (i,j)
        for i in 1:ins.n+1
            for j in 1:ins.n
                if i != j
                    linear_rel[i,j] = callback_value(cb_data, x[i,j])
                end
            end
        end
    
        #Criação do grafo suporte direcionado
        flow_graph = LightGraphs.DiGraph(ins.n+1)
        capacity_matrix = fill(10000.0, ins.n+1, ins.n+1)
    
        for j in 1:ins.n
            for i in 1:ins.n+1
                if linear_rel[i,j] > g_params.eps && i != ins.n+1
                    capacity_matrix[i,j] = linear_rel[i,j]
                    add_edge!(flow_graph, i, j)
                end
            end
        end
    
        #Criação do grafo suporte direcionado ao contrário
        flow_graph_inv = LightGraphs.reverse(flow_graph)
        capacity_matrix_inv = transpose(capacity_matrix)
    
        # Para cada nó t dos nós terminais:
        for t in ins.terminal_nodes
            
            #Rodando o algoritmo de min-cut saindo do vértice artificial para cada vértice terminal
            subtour, complementary, flow_value = LightGraphsFlows.mincut(flow_graph, ins.n+1, t, capacity_matrix, DinicAlgorithm())
    
            # Verificar desigualdade de cutset
            if flow_value < 1 + g_params.eps
                con = @build_constraint(sum(x[i,j] for i in subtour, j in complementary if i != j) >= 1)
                MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            end
    
            #Rodando o algoritmo de min-cut saindo de cada vértice terminal para o vértice artifical no grafo invertido
            subtour, complementary, flow_value = LightGraphsFlows.mincut(flow_graph_inv, t, ins.n+1, capacity_matrix_inv, DinicAlgorithm())
    
            # Verificar desigualdade de cutset
            if flow_value < 1 + g_params.eps
                con = @build_constraint(sum(x[i,j] for i in complementary, j in subtour if i != j ) >= 1)
                MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            end
        end
    end

    # Inicialização de um modelo 
    model = Model(Gurobi.Optimizer)
    
    # Atributos do modelo
    MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "TimeLimit", g_params.time_limit)
    
    #Conjuntos
    V = 1:ins.n
    V0 = 1:ins.n+1
    
    # Definição das variáveis de decisão
    @variables(model, begin
        x[i in 1:ins.n+1, j in 1:ins.n; i != j], (base_name = "x[$i,$j]"), Bin
        0 <= u[i in V0] <= ins.L+1, (base_name = "u[$i]"), Int
        z, Bin
        0 <= w[i in 1:ins.n, j in 1:ins.n] <= 1, Int
        s >= 1
    end)
    
    # Definição da função objetivo
    @objective(model, Min, 2*(s-1) + z)
    
    # Definição das restrições
    @constraints(model, begin
        root, sum(x[ins.n+1,j] for j in V) == z+1
        in_deg[j in ins.terminal_nodes], sum(x[i,j] for i in V0 if i != j) == 1
        st[j in ins.steiner_nodes], sum(x[i,j] for i in V0 if i != j) <= sum(x[j,k] for k in V if j != k)
    end)

    @constraints(model, begin
        mtz[j in V], u[ins.n+1] - u[j] + (ins.L+1)*x[ins.n+1, j] <= ins.L
        mtz1[i in V, j in V; i != j], u[i] - u[j] + (ins.L+1)*x[i,j] + (ins.L-1)*x[j,i] <= ins.L
        mtz2[i in V], u[i] - u[ins.n+1] + (ins.L-1)*x[ins.n+1,i] <= ins.L
    end)
    
    @constraints(model, begin
        budget, sum(ins.c[i,j]*x[i,j] for i in V for j in V if i != j) + sum(ins.c[i,j]*w[i,j] for i in V for j in V if i != j) <= ins.B
        art[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, i]
        art2[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, j]
        central, sum(w[i,j] for i in V, j in V if i != j) == z
        mag[i in V], s >= u[i]
    end)   
    
    
    for i in 1:ins.n
        for j in 1:ins.n
            if ins.c[i,j] < g_params.eps && i != j
                set_upper_bound(x[i,j], 0.0)
                set_lower_bound(x[i,j], 0.0)
                set_upper_bound(w[i,j], 0.0)
                set_lower_bound(w[i,j], 0.0)
            end
        end
    end
    
    set_upper_bound(u[ins.n+1], 0.0)
    set_lower_bound(u[ins.n+1], 0.0)
    
    if g_params.type_of_tree == "fst"
        for j in ins.terminal_nodes
            set_upper_bound(x[ins.n+1,j], 0.0)
            set_lower_bound(x[ins.n+1,j], 0.0)
        end
        
        for i in ins.terminal_nodes
            for j in V
                if i != j
                    set_upper_bound(x[i,j], 0.0)
                    set_lower_bound(x[i,j], 0.0)
                end
            end
        end
    end
    
    MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function)
    optimize!(model)
    
    # Parâmetros do 
    status = termination_status(model)
    nodes = node_count(model)
    obj = objective_bound(model)
    tempo = solve_time(model)
    
    # Relaxação Linear
    relax_integrality(model)
    optimize!(model)
    relax = objective_bound(model)
    
end

