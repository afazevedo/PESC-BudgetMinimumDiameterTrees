using Gurobi, JuMP, MathOptInterface

function solve_spanning_mcf(ins::model_params)
    # Inicialização de um modelo 
    model = Model(Gurobi.Optimizer)

    # Atributos do modelo
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "TimeLimit", g_params.time_limit)

    #Conjuntos
    V = 1:ins.n
    V0 = 1:ins.n+1

    # Definição das variáveis de decisão
    @variables(model, begin
            x[i in 1:ins.n+1, j in 1:ins.n; i != j], (base_name = "x[$i,$j]"), Bin
            y[k in V, i in V0, j in V; i != j] >= 0
            z, Bin
            0 <= w[i in 1:ins.n, j in 1:ins.n] <= 1, Int
            s >= 1
        end)


    # Definição da função objetivo
    @objective(model, Min, 2*(s-1) + z)

    # Definição das restrições
    @constraint(model, root, sum(x[ins.n+1,j] for j in V) == z+1)

    @constraints(model, begin
                flux[k in V], sum(y[k,ins.n+1,j] for j in V) == 1
                flux1[i in V, k in V; i != k], sum(y[k,i,j] for j in V if i != j) - sum(y[k,j,i] for j in V0 if i != j) == 0
                flux3[k in V, i in V0], sum(y[k,k,j] for j in V if k != j) - sum(y[k,j,k] for j in V0 if j != k) == -1
                diam[k in V], sum(y[k,i,j] for i in V0, j in V if i != j) <= ins.L+1
                gonga[i in V0, j in V, k in V; i != j], y[k,i,j] <= x[i,j]
                tst[k in V, j in V; k != j], y[k,ins.n+1,j] <= x[ins.n+1,j] - w[j,k]
            end)

    @constraints(model, begin
            d[k in V], s >= sum(y[k,i,j] for i in V0, j in V if i != j)
            soma, sum(w[i,j] for i in V, j in V if i != j) == z
            art[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, i]
            art2[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, j]
            budget, sum(ins.c[i,j]*x[i,j] + ins.c[i,j]*w[i,j] for i in V, j in V if i != j) <= ins.B
        end)

    #Esparso
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

    optimize!(model)

    peso_total = sum(ins.c[i,j]*value(x[i,j]) for i in V for j in V if i != j) + sum(ins.c[i,j]*value(w[i,j]) for i in V for j in V if i != j)

    # Parâmetros do 
    # status = termination_status(model)
    # nodes = node_count(model)
    obj = objective_bound(model)
    tempo = solve_time(model)

    # Relaxação Linear
    # relax_integrality(model)
    # optimize!(model)
    # relax = objective_bound(model)

    println(obj)
    println(value.(x))
    println(value.(w))

    return obj, peso_total, tempo
end 

function solve_steiner_mcf(ins::model_params)
    # Inicialização de um modelo 
    model = Model(Gurobi.Optimizer)

    # Atributos do modelo
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "TimeLimit", g_params.time_limit)

    #Conjuntos
    V = 1:ins.n
    V0 = 1:ins.n+1

    # Definição das variáveis de decisão
    @variables(model, begin
        x[i in 1:ins.n+1, j in 1:ins.n; i != j], (base_name = "x[$i,$j]"), Bin
        y[k in V, i in V0, j in V; i != j] >= 0
        z, Bin
        0 <= w[i in 1:ins.n, j in 1:ins.n] <= 1, Int
        s >= 1
    end)

    # Definição da função objetivo
    @objective(model, Min, 2*(s-1) + z)

    # Definição das restrições
    @constraint(model, root, sum(x[ins.n+1,j] for j in V) == z+1)

    @constraints(model, begin
                flux[k in ins.terminal_nodes], sum(y[k,ins.n+1,j] for j in V) == 1
                flux1[i in V, k in ins.terminal_nodes; i != k], sum(y[k,i,j] for j in V if i != j) - sum(y[k,j,i] for j in V0 if i != j) == 0
                flux3[k in ins.terminal_nodes, i in V0], sum(y[k,k,j] for j in V if k != j) - sum(y[k,j,k] for j in V0 if j != k) == -1
                diam[k in ins.terminal_nodes], sum(y[k,i,j] for i in V0, j in V if i != j) <= ins.L+1
                gonga[i in V0, j in V, k in ins.terminal_nodes; i != j], y[k,i,j] <= x[i,j]
                tst[k in ins.terminal_nodes, j in V; k != j], y[k,ins.n+1,j] <= x[ins.n+1,j] - w[j,k]
                fluxzero[i in V0, j in V; i != j], x[i,j] <= sum(y[k,i,j] for k in ins.terminal_nodes)
                passagem[j in ins.steiner_nodes], sum(x[i,j] for i in V0 if i != j) <= sum(x[j,k] for k in V if j != k)
            end)

    @constraints(model, begin
            d[k in V], s >= sum(y[k,i,j] for i in V0, j in V if i != j)
            soma, sum(w[i,j] for i in V, j in V if i != j) == z
            art[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, i]
            art2[i in V, j in V; i != j], w[i,j] <= x[ins.n+1, j]
            tssstst[i in V, j in V; i != j], x[i,j] + x[j,i] <= 1
            st[j in ins.steiner_nodes], sum(x[i,j] for i in V0 if i != j) <= sum(x[j,k] for k in V if j != k)
            budget, sum(ins.c[i,j]*x[i,j] + ins.c[i,j]*w[i,j] for i in V, j in V if i != j) <= ins.B
        end)

    #Esparso
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

    optimize!(model)

    peso_total = sum(ins.c[i,j]*value(x[i,j]) for i in V for j in V if i != j) + sum(ins.c[i,j]*value(w[i,j]) for i in V for j in V if i != j)

    # Parâmetros do 
    # status = termination_status(model)
    # nodes = node_count(model)
    obj = objective_bound(model)
    tempo = solve_time(model)

    # Relaxação Linear
    # relax_integrality(model)
    # optimize!(model)
    # relax = objective_bound(model)
    return obj, peso_total, tempo
end 