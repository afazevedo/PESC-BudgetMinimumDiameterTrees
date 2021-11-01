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