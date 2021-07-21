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
    nb_agent = np.random.randint(10,20,1)[0]
    sources = [np.random.randint(0,G_M.vcount(),1)[0]]
    neighbours = G_C.neighbors(sources[0])
    for n in range(nb_agent):
        sources.append(np.random.choice(neighbours,1)[0])
        for i in G_C.neighbors(sources[-1]):
            if not(i in neighbours):
                neighbours.append(i)
    #create an execution 
    print("sources are", sources)
    nb_steps = 50
    last_seen=[[sources[a],sources[a], sources[a],sources[a]] for a in range(nb_agent+1)]
    targets = [sources[i] for i in range(nb_agent+1)]
    for step in range(nb_steps):
        i=0
        while i < nb_agent+1:
            targets[i] = np.random.choice(G_M.neighbors(sources[i]))
            if (mapfalgo.is_connected(targets[:i+1], G_C)) and (not(targets[i] in last_seen[i])):
                last_seen[i]=last_seen[i][1:]+[targets[i]]
                i+=1
    print("targets are", targets)
    return sources, targets

def makeTests(png):
    nb_instance = 1
    if png=="coast.png":
        size = [755, 770]
    elif png == "open.png":
        size = [800,600]
    elif png=="offices.png":
        size = [80,60]
    else :
        size = [944,944]
    results = []
    for n in range(nb_instance):
        radius = np.random.randint(3, size[1]//10)
        physFileName, commFileName = create_graph_from_png.cgfpng(radius, png)
        start = time.time()
        init, target = create_instance(physFileName, commFileName, size)
        sol_divide_and_conquer = mapfalgo.mapf_algo(physFileName, commFileName, init, target)
        end_first = time.time()
        sol_tateo = tateo.tateo(physFileName, commFileName, init, target)
        end = time.time()
        results.append([sol_divide_and_conquer, end_first-start, sol_tateo, end-end_first])
    return results

if __name__ == "__main__":
    png = "offices.png"
    print(makeTests(png))