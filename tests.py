import sys
import mapfalgo
import tateo

import json
import subprocess
import numpy as np
import copy
import igraph
import time 
import os
sys.path.append(os.path.normpath(os.path.join(os.getcwd(),"../cmapf-gui")))
import create_graph_from_png

def create_instance(G_Mname, G_Cname, G_size):
    G_C = igraph.read(G_Cname)
    G_M = igraph.read(G_Mname)
    nb_agent = np.random.randint(10,20,1)[0]-1
    sources = choose_config(G_C,np.random.randint(0,G_M.vcount(),1)[0], nb_agent)
    #create an execution 
    print("sources are", sources)
    #nb_steps = 50
    #last_seen=[[sources[a],sources[a], sources[a],sources[a]] for a in range(nb_agent)]
    fake_first_chosen = sources[0]
    nb_find = 20
    dist_max = 0
    while nb_find>0:
        fake_first = np.random.randint(0,G_M.vcount(),1)[0]
        if mapfalgo.get_distance(G_M, sources[0], fake_first)>dist_max:
            fake_first_chosen=fake_first
            dist_max = mapfalgo.get_distance(G_M, sources[0], fake_first)
        nb_find-=1
    fake_targets= choose_config(G_C, fake_first_chosen, nb_agent)
    targets = tateo_construct_targets(G_M, G_C, sources, fake_targets)
    print("targets are", targets)
    return sources, targets

def choose_config(G_C, first_vertex, nb_agent):
    config = [first_vertex]
    neighbours = G_C.neighbors(config[0])
    for n in range(nb_agent):
        config.append(np.random.choice(neighbours,1)[0])
        for i in G_C.neighbors(config[-1]):
            if not(i in neighbours):
                neighbours.append(i)
    return config

def tateo_construct_targets(G_M, G_C, sources, targets):
    begin = tateo.configuration.Configuration(sources)
    end = tateo.configuration.Configuration(targets)
    first_instance = tateo.instance.Instance(G_M, G_C, begin, end, "astar")
    stack = [begin]
    path = [begin]
    closed = tateo.closed_Tree.Closed_tree()
    nb_steps = 20
    while nb_steps>=0:
        print("nb_steps left", nb_steps)
        current = path[-1]
        closed.add_configuration(current)
        if nb_steps ==0 or current.same(end):
            return current.l_config
        best_child = tateo.find_best_child(first_instance, current, end, closed)
        if not best_child.is_empty():
            stack.insert(0, best_child)
            path.append(best_child)
        else:
            stack.pop(0)
            path.pop()
        nb_steps-=1
    return None

def construct_targets(G_M, G_C, sources, targets):
    return 0

def makeTests(png):
    nb_instance = 1
    if png=="coast.png":
        size = [755, 770]
    elif png == "open.png":
        size = [800,600]
    elif png=="offices.png":
        size = [80,60]
    elif png=="opensmall.png":
        size = [70,59]
    else :
        size = [944,944]
    results = []
    for n in range(nb_instance):
        radius = np.random.randint(3, size[1]//10)
        physFileName, commFileName = create_graph_from_png.cgfpng(radius, png)
        start = time.time()
        init, target = create_instance(physFileName, commFileName, size)
        sol_divide_and_conquer = mapfalgo.mapf_algo(physFileName, commFileName, init, target)
        time_divide_and_conquer = time.time()
        sol_tateo = tateo.tateo(physFileName, commFileName, init, target)
        end = time.time()
        results.append([len(sol_divide_and_conquer[0]), time_divide_and_conquer-start, len(sol_tateo[0]), end-time_divide_and_conquer])
    return results

if __name__ == "__main__":
    png = "opensmall.png"
    print(makeTests(png))
    # physFileName, commFileName = create_graph_from_png.cgfpng(10, "test2.png")
    # # #print(create_instance(physFileName, commFileName, [80,60]))
    # start = time.time()
    # G_M = igraph.read(physFileName)
    # G_C = igraph.read(commFileName)
    # init = [85,83]
    # targets =  [8,6]
    # sol_divide_and_conquer = mapfalgo.mapf_algo(physFileName, commFileName, init, targets)
    # time_div_and_conquer = time.time()
    # sol_tateo = tateo.tateo(physFileName, commFileName, init, targets)
    # end = time.time()
    # print([time_div_and_conquer-start, end-time_div_and_conquer, len(sol_divide_and_conquer[0]),len(sol_tateo[0])])

''' Tests results

Offices
sources = [2092, 1954, 1762, 1713, 1679, 1763, 1586, 1746, 1757, 1689, 1876]
targets = [1765, 1760, 1775, 1487, 1319, 1785, 1254, 1790, 1652, 1780, 1770]
d&c too long

'''