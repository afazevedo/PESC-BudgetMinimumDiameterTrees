using LightGraphs, GraphPlot, Random, Plots, Profile, Cairo

include("../src/read.jl")
include("params.jl")
include("../src/utils.jl")
include("ott.jl")
include("rgh.jl")
include("busca-local.jl")

p = pwd()*"\\instances\\spanning\\"
files = [p*"c_v10_a45_d4.txt" p*"c_v10_a45_d5.txt" p*"c_v10_a45_d6.txt" p*"c_v10_a45_d7.txt" p*"c_v10_a45_d8.txt" p*"c_v10_a45_d10.txt" p*"c_v15_a105_d4.txt" p*"c_v15_a105_d8.txt" p*"c_v20_a190_d4.txt" p*"c_v20_a190_d5.txt" p*"c_v20_a190_d6.txt" p*"c_v20_a190_d7.txt" p*"c_v20_a190_d8.txt" p*"c_v20_a190_d10.txt" p*"c_v25_a300_d4.txt" p*"c_v25_a300_d5.txt" p*"c_v25_a300_d6.txt" p*"c_v25_a300_d8.txt" p*"c_v25_a300_d9.txt"]

ott_results = []
rgh_results = []
# busca_results_ott = []
busca_results_rgh = []

ott_time_results = []
rgh_time_results = []

best_diameter = 10000
soma_tempo = 0

for i in 1:length(files)
    g_params.file_name = files[i]
    ins = read_from_files(g_params.file_name)
    ins.B = budgetCalculator(ins, 0.30)
    hp = heuristic_params{Float64, String, Int64}(SimpleGraph(ins.n), 0.0, 0, [], "ott")
    for j in 1:ins.n
        D = OTT(ins, hp, j)
        tempo = @elapsed OTT(ins, hp, j)
        soma_tempo += tempo
        if D < best_diameter 
            best_diameter = D
        end 
    end 
    push!(ott_time_results, soma_tempo)
    push!(ott_results, best_diameter)
end


best_diameter = 10000
soma_tempo = 0
for i in 1:length(files)
    g_params.file_name = files[i]
    ins = read_from_files(g_params.file_name)
    ins.B = budgetCalculator(ins, 0.30)
    for j in 1:ins.n
        hp = heuristic_params{Float64, String, Int64}(SimpleGraph(ins.n), 0.0, 0, [], "ott")
        D2, tree = rgh(ins, hp)
        tempo = @elapsed rgh(ins, hp)
        soma_tempo += tempo
        if D2 < best_diameter 
            best_diameter = D2
        end 
    end 
    push!(rgh_results, best_diameter)
    push!(rgh_time_results, soma_tempo)
end


best_diameter = 10000
for i in 1:length(files)
    g_params.file_name = files[i]
    ins = read_from_files(g_params.file_name)
    ins.B = budgetCalculator(ins, 0.30)
    for j in 1:ins.n
        hp = heuristic_params{Float64, String, Int64}(SimpleGraph(ins.n), 0.0, 0, [], "ott")
        D, t = rgh(ins, hp)
        if ne(t) == ins.n-1
            hp.spanCost = check_cost(ins, t)
            hp.spanTree = t
            tr, cost, D = decreaseDiameter(ins, hp)
            if hp.tree_diameter < best_diameter 
                best_diameter = hp.tree_diameter
            end 
        end 
    end 
    push!(busca_results_rgh, best_diameter)
end


optimal = [4, 4, 4, 5, 4, 5, 4, 3, 3, 4, 4, 4, 4, 3, 4, 3, 3, 4, 4]


plot(optimal, xticks = 1:1:19, yticks = 2:1:10, xlim = [1,19], ylim=[2, 10], xlabel = "Instância", ylabel = "Tempo(s)",  lw = 2, title = "Teste com B sendo 30% do custo total", label = "Solução Ótima")
plot!(ott_results, lw = 2, label = "Solução Heurística OTT")
plot!(rgh_results, lw = 2, label = "Solução Heurística RGH")
# plot!(rgh_time_results, lw = 2, label = "Heurística RGH")

# plot!(busca_results_rgh, lw = 2, label = "Busca Local RGH")
savefig("teste_busca.png")
 

# todo: botar grafico dos tempos