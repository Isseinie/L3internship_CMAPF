import numpy as np
import math
import igraph
from copy import deepcopy
import heapq
from numpy import random
import mapfalgo

import time 
import os
import sys
sys.path.append(os.path.normpath(os.path.join(os.getcwd(),"../InternshipM1_CMAPF")))
import configuration
import closed_Tree
import instance
import heap_item
import dfs_tateo

def find_best_child(instance, current, targets, closed):
    nb_total_agents = targets.nb_agent
    heap = []
    partial_config = configuration.Configuration([])
    h_cost = dfs_tateo.compute_h(instance, current, partial_config)
    #heap item structure is : ((h + g, -g, config), (g, h, config))
    #Used to store g and h and compare only with g + h
    item = heap_item.Heap_item((h_cost, 0, partial_config), 
                     (0, h_cost, partial_config.copy())
                     )
    heap_item.heappush(heap, item)
    while not heap == [] :
        next_partial = heap_item.heappop(heap)
        (current_g, _, partial_config) = next_partial.item
        num_agent = partial_config.nb_agent
        if num_agent == nb_total_agents:
            if mapfalgo.is_connected(partial_config.l_config, instance.comm_graph) \
                and not closed.is_in(partial_config) \
                and not partial_config.same(current):
                    return partial_config
            continue
        else:
            successors = dfs_tateo.get_successors(instance, current.get_agent_pos(num_agent))
            for node in successors:
                new_config = partial_config.copy()
                new_config.add_agent(node)
                g_cost = current_g + dfs_tateo.compute_distance(
                    current.get_agent_pos(num_agent), node)
                h_cost = dfs_tateo.compute_h(instance, current, new_config)
                if not h_cost == dfs_tateo.m.inf and not g_cost == dfs_tateo.m.inf:
                    item = heap_item.Heap_item((h_cost + g_cost, - g_cost, 
                                      new_config.copy()), 
                                     (g_cost, h_cost, new_config.copy()))
                    heap_item.heappush(heap, item)
    return configuration.Configuration([])

def tateo(G_Mname, G_Cname, sources, targets):
    '''Tateo algorithm
    Input: graphs names, lists of sources and targets
    Output: execution'''
    G_C = igraph.read(G_Cname)
    G_M = igraph.read(G_Mname)
    begin = configuration.Configuration(sources)
    end = configuration.Configuration(targets)
    first_instance = instance.Instance(G_M, G_C, begin, end, "astar")
    stack = [begin]
    path = [begin]
    closed = closed_Tree.Closed_tree()
    while len(stack)>0:#??
        current = path[-1]
        closed.add_configuration(current)
        if current.same(end):
            path_agents = []
            for a in range(len(path[0].l_config)):
                path_a = []
                for e in path :
                    path_a.append(e.l_config[a])
                path_agents.append(path_a)
            return path_agents
        best_child = find_best_child(first_instance, current, end, closed)
        if not best_child.is_empty():
            stack.insert(0, best_child)
            path.append(best_child)
        else:
            stack.pop(0)
            path.pop()
    print("Cannot find a path")
    return None


if __name__== "__main__":
    0
    
