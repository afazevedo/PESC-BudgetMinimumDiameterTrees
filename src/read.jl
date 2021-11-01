using DelimitedFiles, Parsers

include("params.jl")

function read_from_files(arq)
    n = Int(readdlm(arq)[1,1])
    m = Int(readdlm(arq)[1,2])
    # L = n//2
    c = zeros(n,n)
    terminal_nodes = []
    steiner_nodes = []
    
    if g_params.type_of_tree == "steiner"
        
        for k in 2:m+1
            i = Int(readdlm(arq)[k, 1])
            j = Int(readdlm(arq)[k, 2])
            c[i,j] = readdlm(arq)[k, 3]
            c[j, i] = c[i, j]
        end

        linhas = countlines(arq)
        number_of_terminals = readdlm(arq)[linhas-1,1]
        terminal = zeros(Int64, number_of_terminals)
        
        for i in 1:number_of_terminals
            terminal[i] = Int(readdlm(arq)[linhas, i])
        end
        
        count_1 = 1
        count_2 = 1
        for i in 1:n
            if i in terminal
                append!(terminal_nodes,i)
                count_1 += 1
            else
                append!(steiner_nodes,i)
                count_2 += 1
            end
        end
    
        ins = model_params{Float64, Int64}(n, m, 0, terminal_nodes, steiner_nodes, c, 0)
        return ins
    end

    for k in 2:m+1
        i = Int(readdlm(arq)[k, 1])
        j = Int(readdlm(arq)[k, 2])
        c[i+1,j+1] = readdlm(arq)[k, 3]
        c[j+1, i+1] = c[i+1, j+1]
    end

    ins = model_params{Float64, Int64}(n, m, 0, terminal_nodes, steiner_nodes, c, 0)
    return ins
end 
