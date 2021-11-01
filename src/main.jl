
using Profile

include("params.jl")
include("read.jl")
include("results.jl")
include("utils.jl")
include("../models/mtz.jl")
include("../models/mcf.jl")
include("../models/hop-cut.jl")

function main(L, spanTree)
    if g_params.type_of_model == "mtz"
        if g_params.type_of_tree == "spanning"
            ins = read_from_files(g_params.file_name)
            ins.B = budgetCalculator(ins, 0.20)
            if g_params.heuristic
                ins.L = L
            else
                if ins.n % 2 == 0
                    ins.L = ins.n/2
                else 
                    ins.L = (ins.n-1)/2
                end 
            end 
            obj, cost, tempo = solve_spanning_mtz(ins, spanTree)
            results(obj, cost, tempo, g_params.file_name, ins.B, ins.L)
        end  
    elseif g_params.type_of_model == "mcf"
        if g_params.type_of_tree == "spanning"
            ins = read_from_files(g_params.file_name)
            ins.B = budgetCalculator(ins, 0.20)
            if g_params.heuristic
                ins.L = L
            else
                if ins.n % 2 == 0
                    ins.L = ins.n/2
                else 
                    ins.L = (ins.n-1)/2
                end 
            end 
            obj, cost, tempo = solve_spanning_mcf(ins)
            results(obj, cost, tempo, g_params.file_name, ins.B, ins.L)
        end  
    end 
    # elseif g_params.type_of_model == "hopcut"
    #     if g_params.type_of_tree == "spanning"
    #         #spanning
    #     elseif g_params.type_of_tree == "steiner"
    #         #steiner
    #     else
    #         #full steiner
    #     end 
    # end 
end

p = pwd()*"\\instances\\$type_of_tree\\"
files = [p*"c_v10_a45_d4.txt" p*"c_v10_a45_d5.txt" p*"c_v10_a45_d6.txt" p*"c_v10_a45_d7.txt" p*"c_v10_a45_d8.txt" p*"c_v10_a45_d10.txt" p*"c_v15_a105_d4.txt" p*"c_v15_a105_d8.txt" p*"c_v20_a190_d4.txt" p*"c_v20_a190_d5.txt" p*"c_v20_a190_d6.txt" p*"c_v20_a190_d7.txt" p*"c_v20_a190_d8.txt" p*"c_v20_a190_d10.txt" p*"c_v25_a300_d4.txt" p*"c_v25_a300_d5.txt" p*"c_v25_a300_d6.txt" p*"c_v25_a300_d8.txt" p*"c_v25_a300_d9.txt"]
# L_optimal = [0, 0, 0, 0, 0, 0, 3, 2, 3, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2]

for i in 7:length(files)
    # @profile main(ins.L, pt.spanTree)
    @profile main(L_optimal[i], pt.spanTree)
end




