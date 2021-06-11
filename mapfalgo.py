#from pysat.solvers import Solver
import numpy as np
import igraph
from copy import deepcopy

#2 graphs : movement & connection
#A : list of sources and targets (agent a goes from sources[a] to targets[a])
#exec : list of configurations [[a0_t,...] for each t] (or list of paths)

def decoupled_exec(G_M, sources, targets) :
    '''algo shortest path for each agent : A*
    for now, we use a less efficient algorithm
    Output : execution P'''
    paths = []
    for a in range(0, len(sources)) :
       paths.append(extract_path_from_pred(get_pred(G_M, sources[a], targets[a]),sources[a], targets[a])) 
    configs = [targets for t in max(paths.map(len))]
    for i_a in range(0, len(paths)) :
        for t in len(paths[i_a]):
            configs[t][i_a] = paths[i_a][t]
    return configs

def extract_path_from_pred(pred, source, dest) :
    if source == dest :
        return [source]
    else :
        extract_path_from_pred(pred, source, pred[dest])+ [dest]


def get_pred(G_M, source, dest) :
    pred = [-1 for x in len(G_M.vcount())]
    queue = [source]
    visited = [False for x in len(G_M.vcount())]
    while len(queue) > 0 :
        x = queue.pop()
        if x == dest :
            return pred
        neighbours = G_M.neighbors(x)
        for n in neighbours :
            if not(visited[n]) :
                queue.append(n)
                pred[n] = x
    return None




def nb_conflicts(exec, G_C) :
    '''Use an algorithm already created or use A_ordered to see for each configuration if there is a disconnected agent 
    Input : execution exec
    Output : int number of conflicts'''
    return len(filter(lambda x : not x, map(lambda config : is_connected(config, G_C), exec)))
    

def is_connected(config, G_C):
    '''return True if the configuration is connected '''
    queue = [0]
    visited = [False for a in config]
    while len(queue) > 0 :
        x = queue.pop()
        for i_a in range(0, len(config)):
            if G_C.are_connected(config[x], config[i_a]) :
                if not(visited[i_a]):
                    queue.append(i_a)
    return all(visited)




def is_ordered_connected(G_C, A_ordered_id, i, t, exec): 
    '''true if a_i is connected to agents a_0... a_i-1 at time t 
    look in the list of neighbours of a_i at time t if there is a_j, j<i'''
    if i == 0 :
        return True
    neighbours = G_C.neighbors(exec[t][i], mode = "all")
    for j in range(0, i):
        if exec[t][A_ordered_id[j]] in neighbours:
            return True
    return False

def pick_time_with_conflict(exec, G_C) : 
    '''Choose t around the middle of the execution, with conflicts '''
    for i in range(0,len(exec)//2):
        if nb_conflicts(exec[len(exec)//2 + i], G_C) > 0 :
            return len(exec)//2 +i 
        elif nb_conflicts(exec[len(exec)//2 - i], G_C) > 0 :
            return len(exec)//2 -i 
    return len(exec)//2

def choose_order(G_C, A) :
    '''Choose an order of agents
    Output : list of index '''
    i = np.random(len(A))
    A_ordered_id = [i]
    #BFS on A depending on the initial configuration
    queue = G_C.neighbors(A[i], mode = "all")
    while len(A_ordered_id)<len(A):
        v = queue.pop(0)
        for j in range(0, len(A)):
            if A[j]==v :
                A_ordered_id.append(j)
                queue += G_C.neighbors(A[j], mode = "all")
    return A_ordered_id




def path_of_length(G_M, start, arrival, t):
    if t == 0 :
        if start == arrival : 
            return [arrival] 
        else : 
            return []
    else :
        Neighbours = G_M.neighbors(start, mode = "all")
        return [([n]+x for x in path_of_length(G_M, n, arrival, t-1)) for n in Neighbours]

def pick_path_of_length(G_M, start, arrival, t):
    paths = path_of_length(G_M, start, arrival, t)
    res = []
    for p in paths :
        if len(p) == t :
            res.append(p)
    return res


def choose_best_neighbour(G_M, G_C, sources, targets, A_ordered_id, i, t, exec):
    '''Choose a neighbour u of a_0...a_i-1 with path of length t from s_i to u, minimize d(u, g_i) and nb of conflicts '''
    Neighbours = []
    for j in range(0,i):
        Neighbours+= G_C.neighbors(exec[t][A_ordered_id[j]], mode = "all")
    min_dist = len(exec)
    min_nb_conflicts = nb_conflicts(exec, G_C)
    best = Neighbours[0]
    for u in Neighbours :
        p_final = decoupled_exec(G_M, G_C, [u], [targets[i]])
        paths = pick_path_of_length(G_M, sources[i], u, t)
        for p in paths:
            exec_copy = deepcopy(exec)
            for time in range(t)):
                exec_copy[time][i] = p[time]
            for time in range(len(p_final)):
                exec_copy[time+t][i] = p_final[time]
            if len(p_final) < min_dist :
                if nb_conflicts(exec_copy, G_C) < min_nb_conflicts:
                    best = u
                    min_dist = len(p_final)
                    min_nb_conflicts = nb_conflicts(exec_copy, G_C)
        #if there is one, then compute the distance d(u, g_i) and the nb of conflicts: 
        #if it's less than min_dist then and min_nb_conflicts then replace best by u
    return best




def update(G_M, exec, u, i, t):
    '''update exec with a_i going through u at t'''
    sources_first = [exec[0][j] for j in range(0, len(exec[0]))]
    targets_first = [exec[t][j] for j in range(0, len(exec[0]))]
    targets_first[i] = u
    sources_second = [exec[t][j] for j in range(0, len(exec[0]))]
    targets_second = [exec[len(exec)][j] for j in range(0, len(exec[0]))]
    sources_second[i] = u
    return decoupled_exec(G_M, sources_first, targets_first) + decoupled_exec(G_M, sources_second, targets_second)




def mapfdivideandconquer(G_M, G_C, sources, targets):
    exec = decoupled_exec(G_M, sources, targets)
    nb_it = 0
    while nb_it < 5 : #number of attempts to find a better P
        A_ordered_id = choose_order(G_C, sources) 
        exec_changed = aux_divide(exec, sources, targets, A_ordered_id, G_C, G_M, 0)
        if nb_conflicts(exec_changed, G_C) == 0 :
            return exec_changed
        else :
            nb_it+=1 #if exec_changed has less conflicts, exec is replaced by exec_changed
    return exec

def aux_divide(exec, sources, targets, A_ordered_id, G_C, G_M, n):
    if n < 5: #number of recursive calls = 5
        t = pick_time_with_conflict(exec, G_C)
        exec_changed = deepcopy(exec)
        for i in A_ordered_id:
            if not(is_ordered_connected(G_C, A_ordered_id, i, t, exec_changed)):
                u = choose_best_neighbour(G_M,G_C, sources, A_ordered_id, i, t, exec_changed)
                update(G_M, exec_changed, u, i, t) #update of exec_i
            if nb_conflicts(exec_changed, G_C) < nb_conflicts(exec):
                exec = exec_changed
        return aux_divide(exec[::len(exec)//2],sources, targets, A_ordered_id, G_C, G_M, n+1) + aux_divide(exec[len(exec)//2::], sources, targets, A_ordered_id, G_C, G_M, n+1)
    else :
        return exec







'''def last_is_disconnected(G_M, G_C, A, P):
    return False, 0

def find_path(G_M, G_C, A, n, constraints):
    return Solver()

def extract_exec(G_M, G_C, A, constraints):
    return 0

def pick_config(G_M, G_C, A, t):
    return [[0] for a in A]

def MAPF1(G_M, G_C, A) :
    constraints = []
    P = decoupled_exec(G_M, G_C, A)
    n = max([len(Pai) for Pai in P])
    while has_conflict(P):
        agent, t = pick_disconnected(G_M, G_C, A)
        constraints.append(is_connected(G_M, G_C, A, agent, t))
        while(find_path(G_M, G_C, A, n, constraints).solve() == False) :
            n+= 1
        P = extract_exec(G_M, G_C, A, constraints)
    return P

def MAPFmalin(G_M, G_C, A, final_state):
    constraints = []
    P = decoupled_exec(G_M, G_C, A)
    n = max([len(Pai) for Pai in P])
    if not(has_conflict(P)) :
        return P
    else :
        agent, t = pick_disconnected(G_M, G_C, A)
        constraints.append(is_connected(G_M, G_C, A, agent, t))
        while n < NMAX :
            while(find_path(G_M, G_C, A, n, constraints).solve() == False) :
                n+= 1
            P_extr = extract_exec(G_M, G_C, A, n, constraints)
            P_final = MAPFmalin(G_M, G_C, A, P_extr[-1]) + MAPFmalin(G_M, G_C, P_extr[-1], final_state)
            if not(has_conflict(P_final)) :
                return P_final
            else :
                constraints.append(pick_config(G_M, G_C, A, t) != P_extr[t])'''

