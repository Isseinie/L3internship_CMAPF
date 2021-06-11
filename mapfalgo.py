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
    Output : execution configs'''
    paths = []
    for a in range(0, len(sources)) :
       paths.append(extract_path_from_pred(get_pred(G_M, sources[a], targets[a]),sources[a], targets[a])) 
    configs = [deepcopy(targets) for t in range(max(map(len, paths)))]
    for i_a in range(0, len(paths)) :
        for t in range(len(paths[i_a])):
            configs[t][i_a] = paths[i_a][t]
    return configs

def extract_path_from_pred(pred, source, dest) :
    if source == dest :
        return [source]
    else :
        return extract_path_from_pred(pred, source, pred[dest])+ [dest]


def get_pred(G_M, source, dest) :
    pred = [-1 for x in range(G_M.vcount())]
    queue = [source]
    visited = [False for x in range(G_M.vcount())]
    while len(queue) > 0 :
        x = queue.pop(0)
        visited[x] = True
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
    print(exec)
    return 0
    is_connected_array = map(lambda config : is_connected(config, G_C), exec)
    return len(list(filter(lambda x : not x, is_connected_array)))
    

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
    i = np.random.randint(0, len(A), 1)[0]
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




'''def path_of_length(G_M, start, arrival, t):
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
    return res'''


def execution_with_best_neighbour(G_M, G_C, sources, targets, A_ordered_id, i, t, exec):
    '''Choose a neighbour u of a_0...a_i-1 with path of length t from s_i to u, minimize d(u, g_i) and nb of conflicts 
    Output : execution with a_i going through u at t'''
    Neighbours = []
    for j in range(0,i):
        Neighbours+= G_C.neighbors(exec[t][A_ordered_id[j]], mode = "all")
    min_dist = len(exec)
    min_nb_conflicts = nb_conflicts(exec, G_C)
    best = Neighbours[0]
    best_exec = deepcopy(exec)
    for u in Neighbours :
        dist_u_gi = len(decoupled_exec(G_M, [u], [targets[i]])) #d(u, g_i)
        '''paths = pick_path_of_length(G_M, sources[i], u, t) #get possible paths for a_i through u at t
        for p in paths:
            exec_first = [deepcopy(exec[time]) for time in range(t)] #execution from beginning to time t
            for time in range(t):
                exec_first[time][i] = p[time]#replace path of a_i by the found one'''
        sources_first = [sources[j] for j in range(len(sources))]
        targets_first = [exec[t][j] for j in range(len(sources))]
        targets_first[i] = u
        sources_second = [exec[t][j] for j in range(len(sources))]
        targets_second = [targets[j] for j in range(len(sources))]
        sources_second[i] = u
        exec_tested = decoupled_exec(G_M, sources_first, targets_first)+decoupled_exec(G_M, sources_second, targets_second)
        if dist_u_gi < min_dist :
            if nb_conflicts(exec_tested, G_C) < min_nb_conflicts:
                best = u
                min_dist = dist_u_gi
                min_nb_conflicts = nb_conflicts(exec_tested, G_C)
                best_exec = exec_tested
        #if there is one, then compute the distance d(u, g_i) and the nb of conflicts: 
        #if it's less than min_dist then and min_nb_conflicts then replace best by u
    return best_exec


###Algorithm must return list of paths

def mapfdivideandconquer(G_Mname, G_Cname, sources, targets):
    '''This algorithm's method is divide and conquer '''
    G_C = igraph.read(G_Cname)
    G_M = igraph.read(G_Mname)
    #exec = decoupled_exec(G_M, sources, targets)
    nb_it = 0
    while nb_it < 5 : #number of attempts to find a better P
        A_ordered_id = choose_order(G_C, sources) 
        exec_changed = aux_divide(sources, targets, A_ordered_id, G_C, G_M, 5) 
        if nb_conflicts(exec_changed, G_C) == 0 :
            return exec_changed
        else :
            nb_it+=1 #if exec_changed has less conflicts, exec is replaced by exec_changed
    return exec

def aux_divide(sources, targets, A_ordered_id, G_C, G_M, n):
    '''This function fix the connection problem around the middle of the execution, then redo it for each part'''
    exec = decoupled_exec(G_M, sources, targets)
    if nb_conflicts(exec, G_C) == 0 :
        return exec
    if n >0: #number of recursive calls = 5
        t = pick_time_with_conflict(exec, G_C)
        for i in A_ordered_id:
            if not(is_ordered_connected(G_C, A_ordered_id, i, t, exec)):
                exec_changed = execution_with_best_neighbour(G_M,G_C, sources, targets, A_ordered_id, i, t, exec) #update of exec_i
                if nb_conflicts(exec_changed, G_C) < nb_conflicts(exec, G_C):
                    exec = exec_changed
        return aux_divide(sources, exec[t], A_ordered_id, G_C, G_M, n-1) + aux_divide(exec[t], targets, A_ordered_id, G_C, G_M, n-1)
    else :
        return exec


### Tests

'''G_comm = igraph.read("map1.png_comm_uniform_grid_1_range_6.graphml")
G_mov = igraph.read("map1.png_phys_uniform_grid_1_range_6.graphml")

sources = [47, 14]
targets = [55, 20]

pred_ex = get_pred(G_mov, 47, 55)
path_ex = extract_path_from_pred(pred_ex, 47, 55)

exec_ex = decoupled_exec(G_mov, sources, targets)
print(exec_ex)'''





### Abandonned functions, algorithms


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


'''def update(G_M, exec, u, i, t):
    ''update exec with a_i going through u at t''
    sources_first = [exec[0][j] for j in range(0, len(exec[0]))]
    targets_first = [exec[t][j] for j in range(0, len(exec[0]))]
    targets_first[i] = u
    sources_second = [exec[t][j] for j in range(0, len(exec[0]))]
    targets_second = [exec[len(exec)][j] for j in range(0, len(exec[0]))]
    sources_second[i] = u
    return decoupled_exec(G_M, sources_first, targets_first) + decoupled_exec(G_M, sources_second, targets_second)'''