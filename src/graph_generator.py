"""
    Waxman model - graph generator for TOC paper

"""
import sys
import networkx as nx
import random
import numpy
import math
from math import ceil

RANGE = 300
G = nx.Graph()
vis = {}

def dfs(node):
    global vis
    global G
    if vis[node] == True:
        return
    vis[node] = True
    for v in G.neighbors(node):
        if vis[v] == False:
            dfs(v)
    return

def is_connect():
    global vis
    global G
    # print(type(G.nodes()))
    vis = {node : False for node in G.nodes()}
    for node in G.nodes():
        dfs(node)
        break
    for node in G.nodes():
        if(vis[node] == False):
            return False
    return True

def dist(p1, p2):
    (x1, y1) = p1
    (x2, y2) = p2
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1 / 2)

def link_prob(entangle_lambda, dis, times):
    one_prob = math.exp(-entangle_lambda * dis)
    print("one_prob =", 1 - ((1 - one_prob) ** times), file=sys.stderr)
    return 1 - ((1 - one_prob) ** times)

if len(sys.argv) <= 2:
    print("missing argv")
    sys.exit()

filename = sys.argv[1]
num_of_node = int(sys.argv[2])
# Paper Eq.2: Fe(u,v) = 1/4 + 3/4 * exp(-gamma * l(u,v))
#   fid_ratio = exp(-gamma * l)，之後在 Graph.cpp 中做 min/max 縮放
# Paper Eq.13: Pr(u,v) = 1 - (1 - exp(-lambda * l(u,v)))^xi
gamma = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
entangle_lambda_param = float(sys.argv[4]) if len(sys.argv) > 4 else 0.045
xi = int(sys.argv[5]) if len(sys.argv) > 5 else 8

print("======== generating graph ========", file=sys.stderr)
print("filename =", filename, file=sys.stderr)
print("num_of_node =", num_of_node, file=sys.stderr)
# print("min_fidelity =", min_fidelity, ", max_fidelity =", max_fidelity, file=sys.stderr)
# print("min_memory_cnt =", min_memory_cnt, ", max_memory_cnt =", max_memory_cnt, file=sys.stderr)

while True:
    G = nx.waxman_graph(num_of_node, beta=0.85, alpha=0.08, domain=(0, 0, 0.5, 1))
    # G = nx.waxman_graph(num_of_node, beta=0.85, alpha=10, domain=(0, 0, 0.5, 1))
    positions = nx.get_node_attributes(G, 'pos')
    add_edge = []
    # connect the u's nearest vertex with u
    for u in range(G.order()-1):
        mi_dist = dist(positions[u], positions[G.order()-1])
        mi_idx = G.order()-1
        for v in range(u+1, G.order()):
            if(G.has_edge(u, v)):
                continue
            if(mi_dist > dist(positions[u], positions[v])):
                mi_dist = dist(positions[u], positions[v])
                mi_idx = v
        if(G.has_edge(u, mi_idx)):
            continue
        add_edge.append((u, mi_idx))

    for e in add_edge:
        (u, v) = e
        G.add_edge(u, v)

    if is_connect():
        break 
    else:
        print("topo is not connected", file=sys.stderr)

path = filename
with open(path, 'w') as f:
    positions = nx.get_node_attributes(G, 'pos')
    # write node
    print(num_of_node, file=f)
    for n in G.nodes():
        (x, y) = positions[n]
        pos_x = str(x*RANGE)
        pos_y = str(y*RANGE)
        num_of_memory = random.randint(-1, 1)
        print(num_of_memory, file = f)
    
    # write edge
    num_of_edge = 0
    for e in G.edges():
        if e[0] != e[1]:
            num_of_edge += 1
    print(num_of_edge, file=f)
    avg_l = 0
    for e in G.edges():
        if e[0] != e[1]:
            e0 = str(e[0])
            e1 = str(e[1])
            dis = RANGE * dist(positions[e[0]], positions[e[1]])  # distance
            # F = random.random()*(max_fidelity-min_fidelity) + min_fidelity  # fidelity
            # ratio_list=[0.98,0.7]
            # ratio=numpy.random.choice(ratio_list,p=[0.3,0.7])
            # ratio = numpy.random.normal(0.7, 0.1)
            #dif = abs(1 - ratio)
            #ratio = 1 - dif
            #if ratio > 0.95:
            #    ratio = 0.95
            #if ratio < 0.55:
            #    ratio = 0.55
            # Paper Eq.2: fid_ratio = exp(-gamma * l), 寫入後由 Graph.cpp 做 min/max 縮放
            fid_ratio = math.exp(-gamma * dis)
            fid_ratio = max(0.01, min(0.99, fid_ratio))

            # Paper Eq.13: Pr(u,v) = 1 - (1 - exp(-lambda * l))^xi
            one_prob = math.exp(-entangle_lambda_param * dis)
            entangle_prob = 1.0 - ((1.0 - one_prob) ** xi)
            entangle_prob = max(0.001, min(0.999, entangle_prob))

            print(e0 + " " + e1 + " " + str(entangle_prob) + " " + str(fid_ratio), file=f)
            avg_l += dis
    avg_l /= num_of_edge

print("num_of_edge =", num_of_edge, file=sys.stderr)
print("avg_edge_len =", avg_l, file=sys.stderr)
print("\n======== graph generate finished ! ========", file=sys.stderr)


# print(prob(entangle_lambda, 150))