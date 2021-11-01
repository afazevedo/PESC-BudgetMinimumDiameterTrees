include("../src/params.jl")
include("params.jl")
include("../src/read.jl")
include("../src/utils.jl")

g_params.file_name = pwd()*"\\instances\\$type_of_tree\\c_v10_a45_d4.txt"
ins = read_from_files(g_params.file_name)


# function subproblema(c, w, num_fac, num_cli)
#     x = zeros(num_fac, num_cli)
#     y = zeros(num_fac) 

#     melhor = MAX_ITER
#     best_fac = 0
    
#     if 1-w <= 0
#         for i in 1:num_fac
#             if sum(ins.c[i,j] for j in 1:num_cli) <= melhor
#                 melhor = sum(ins.c[i,j] for j in 1:num_cli)
#                 best_fac = i
#             end
#         end 

#         y[best_fac] = 1
    
#         for j in 1:num_cli 
#             x[best_fac,j] = 1
#         end
#         lb = -(w*ins.B)
#         for j in 1:num_cli 
#             if x[best_fac, j] == 1 && j != best_fac
#                 lb = lb + ins.c[best_fac,j]*(1-w)
#             end 
#         end 
#     else 
#         lb = -(w*ins.B)
#     end 

#     return x, y, lb/2
# end

function subproblema(c, u, K, num_fac, num_cli)
    x = zeros(num_fac, num_cli)
    y = zeros(num_fac) #abrir ou não a facilidade
    v = zeros(num_fac) #custo lagrangiano (ajuda no subproblema)

    for i in 1:num_fac
        v[i] = 0.0
        for j in 1:num_cli
            v[i] = v[i] + min(0, c[i,j] - u[j])
        end 
    end

    idx = sortperm(v) #lista dos índices
    for i in idx[1:K]
        y[i] = 1
    end 

    for i in 1:num_fac
        for j in 1:num_cli
            if y[i] == 1 && c[i,j] - u[j] < 0
                x[i,j] = 1
            end
        end 
    end

    lb = 0.0 #lower bound
    for j in 1:num_cli
        lb = lb + u[j]
        for i in 1:num_fac
            if x[i,j] == 1
                lb = lb + c[i,j] - u[j]
            end 
        end
    end
   
    return x, y, lb
end

function upper_bound(y, num_fac, num_cli, c)
    x = zeros(num_fac, num_cli)
    
    for j in 1:num_cli
        idx = argmin(c[:,j] + (1 .- y) .* maximum(c))
        x[idx, j] = 1
    end

    ub = 0.0
    for i in 1:num_fac
        for j in 1:num_cli
            if x[i,j] == 1
                ub = ub + c[i,j]*x[i,j]
            end 
        end
    end

    return ub, x
end


maxIter = 100
p_i = 2
pi_min = 0.0001

best_lim_inf = -Inf
best_lim_sup = Inf

num_fac = size(ins.c, 1)
num_cli = size(ins.c, 2)

x_best = zeros(num_fac, num_cli)
y_best = zeros(num_fac)

u = zeros(num_cli)
improve = 0
ins.B = 252

for k in 1:maxIter    
    # x_sub, y_sub, z = subproblema(ins.c, w, num_fac, num_cli)
    x_sub, y_sub, z = subproblema(ins.c, u, 1, num_fac, num_cli)
    println("z: ", z)

    if z > best_lim_inf
        best_lim_inf = z
        improve = 0
        # x_best = x_sub
    else
        improve += 1
    end

    ub, x_up = upper_bound(y_sub, num_fac, num_cli, ins.c)
    # ub = 252

    if ub < best_lim_sup
        best_lim_sup = ub 
        x_best = x_sub
        y_best = y_sub
    end 

    if best_lim_sup - best_lim_inf < 1
        println("Parando por otimalidade (z_up == z_low) - iteração ", k)
        break
    end

    if improve >= maxIter/20
        p_i = p_i/2
        improve = 0
        if p_i < pi_min
            println("Parando por pi pequeno (iteração ", k, ")")
            break
        end
    end

    s = zeros(num_cli)
    norm = 0
    for j in 1:num_cli
        s[j] = 1
        for i in 1:num_fac
            s[j] -= x_sub[i,j]
        end
        norm += s[j]
    end 

    mi = p_i*((ub-z)/(norm)^2)

    for j in 1:num_cli
        u[j] = u[j] + mi*s[j]
    end 
end

sum(ins.c[i,j]*x_best[i,j] for i in 1:num_fac for j in 1:num_cli)

println(best_lim_inf)
println(best_lim_sup)

