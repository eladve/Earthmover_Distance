import numpy as np
import networkx as nx

def ell2(a, b):
    return np.linalg.norm(a-b)
    
S = list(map(np.array, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]))
T = list(map(np.array, [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]]))

def compute_optimal_transport(S, T, underlying_distance = ell2, S_weights = None, T_weights = None):
    """
    input:
    S and T should each be a list of numpy vectors, all unique (no repeats! If you want a point to appear more than once, just increase its weight)
    S_weights should be a list of positive integers of length len(S). Ditto T_weights. If you're not sure, just use the default. The weight is the "number of times" that a point appears, i.e. the mass of the distribution at that point.
    underlying_distance should be a function taking two numpy vectors and returning their distance -- this should be the underlying distance for the purposes of the EMD definition you'd like to use. If you're not sure, just use the default ell2
    Idea for making this more efficient: instead of applying the underlying_distance function to each pair by itself, run something like scipy.spatial.distance.pdist to compute them all at the same time. Might be more efficient.
    """
    dim = len(S[0])
    for v in S: assert len(v) == dim
    for u in T: assert len(u) == dim
    
    if S_weights is None:
        S_weights = [1.0] * len(S)
    total_S_weight = sum(S_weights)
    if T_weights is None:
        T_weights = [1.0] * len(T)
    total_T_weight = sum(T_weights)
    
    # now let's build the distance matrix:
    dist_matrix = []
    for v in S:
        temp = []
        for u in T:
            temp.append(underlying_distance(v, u))
        dist_matrix.append(temp)
    
    # now let's define the graph:
    # (documentation: https://networkx.github.io/documentation/networkx-1.9/tutorial/tutorial.html )
    # Note that we set all capacities to be integers in order to make the algorithms happy,
    # but we'll then need to divide the final cost by total_S_weight*total_T_weight so that the result is indifferent
    # to the cardinalities and weights of the point-sets.
    
    G=nx.DiGraph()
    G.add_node("src") 
    G.add_node("sink")
    for i in range(len(S)):
        vtx_name = ("S", i)
        G.add_node(vtx_name)
        G.add_edge("src", vtx_name, capacity = total_T_weight * S_weights[i], weight = 0)
    for j in range(len(T)):
        vtx_name = ("T", j)
        G.add_node(vtx_name)
        G.add_edge(vtx_name, "sink", capacity = total_S_weight * T_weights[j], weight = 0)
    for i in range(len(S)):
        for j in range(len(T)):
            G.add_edge(("S", i), ("T", j), weight = dist_matrix[i][j]) # infinite capacity edge
            
    # now run the flow algorithm:
    # (documentation: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.flow.max_flow_min_cost.html )
    
    mincostFlow_unnormalized = nx.max_flow_min_cost(G, "src", "sink")
    # sanity check -- check that indeed total_S_weight*total_T_weight units of flow have gone through
    amount_of_flow_unnormalized = sum([mincostFlow_unnormalized[v]["sink"] for v in G.predecessors("sink")])
    assert amount_of_flow_unnormalized == total_S_weight * total_T_weight

    # calculate normalized cost of flow:
    mincost_unnormalized = nx.cost_of_flow(G, mincostFlow_unnormalized)
    mincost_normalized = mincost_unnormalized / (total_S_weight * total_T_weight)
    
    # Now remove the source and sink, replace vertex names by their corresponding vectors, and renormalize
    resulting_flow = dict()
    for i in range(len(S)):
        source_vec = tuple(S[i])
        source_vec_weight = S_weights[i]
        outgoing_from_source_vec = mincostFlow_unnormalized[("S", i)]
        result_for_source_vec = dict()
        for u in outgoing_from_source_vec.keys():
            target_vec = tuple(T[u[1]])
            target_vec_weight = T_weights[u[1]]
            flow_amount = outgoing_from_source_vec[u]
            if flow_amount > 0:
                result_for_source_vec[target_vec] = flow_amount / total_T_weight
        resulting_flow[source_vec] = result_for_source_vec
        
    return resulting_flow, mincost_normalized

resulting_flow, EMD_value = compute_optimal_transport(S, T)

print("the EMD value is", EMD_value)

print ("here are the details of the optimal transport:")
print ("==============================================")
for v in resulting_flow.keys():
    print(v, "gets transported to the following:")
    outgoing_from_v = resulting_flow[v]
    for u in outgoing_from_v.keys():
        flow_amount = outgoing_from_v[u]
        if flow_amount > 0:
            print (round(flow_amount, 2), "units get transported to", u)
    print()
