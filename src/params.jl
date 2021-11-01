mutable struct general_params{F<:Float64, S<:String, B<:Bool}
    file_name::S        # Nome do arquivo
    eps::F          # Controle de erro
    time_limit::F       # Tempo limite
    type_of_model::S    # Tipo de modelo
    type_of_tree::S     # Tipo de árvore
    warm_start::B       # Usar ou não um warm start
    heuristic::B        # Usar ou não a heurística
    max_cuts::F        # Número máximo de cortes
end

mutable struct model_params{F<:Float64, I<:Int64}
    n::I
    m::I
    L::I
    terminal_nodes::Array{I,1}
    steiner_nodes::Array{I,1}
    c::Array{F, 2}
    B::F
end

type_of_tree = "spanning"
# file = pwd()*"\\instances\\$type_of_tree\\c_v10_a45_d4.txt"
file = pwd()*"\\instances\\spanning\\c_v25_a300_d4.txt"
eps = 0.0001
time_limit = 1200
type_of_model = "mcf"
warm_start = 0
heuristic = 1
max_cuts = 2

g_params = general_params{Float64, String, Bool}(file, eps, time_limit, type_of_model, type_of_tree, warm_start, heuristic, max_cuts)
