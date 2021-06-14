#from pysat.solvers import Solver
import numpy as np
import igraph
from copy import deepcopy


#2 graphs : movement & connection
#A : list of sources and targets (agent a goes from sources[a] to targets[a])
#exec : list of paths for each agent

def decoupled_exec(G_M, sources, targets) :
    '''algo shortest path for each agent : A*
    for now, we use a less efficient algorithm
    Output : execution configs'''
    paths = []
    for a in range(0, len(sources)) :
       paths.append(extract_path_from_pred(get_pred(G_M, sources[a], targets[a]),sources[a], targets[a]))
    max_t = max(map(len, paths))
    for p in paths :
        while len(p) < max_t :
            p.append(p[len(p)-1]) #we want same-length paths for each agent, so we make them wait at their arrival to complete their path
    return paths

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
    list_config = [[0 for i_a in range(len(exec))] for t in range(len(exec[0]))]
    for i_a in range(len(exec)) :
        for t in range(len(exec[i_a])):
            list_config[t][i_a] = 1
            list_config[t][i_a] = exec[i_a][t]
    is_connected_array = map(lambda config : is_connected(config, G_C), list_config)
    return len(list(filter(lambda x : not x, is_connected_array)))
    

def is_connected(config, G_C):
    '''return True if the configuration is connected '''
    queue = [0]
    visited = [False for a in config]
    visited[0] = True
    while len(queue) > 0 :
        x = queue.pop()
        for i_a in range(0, len(config)):
            if G_C.are_connected(config[x], config[i_a]) :
                if not(visited[i_a]):
                    queue.append(i_a)
                    visited[i_a] = True
    return all(visited)




def is_ordered_connected(G_C, A_ordered_id, i, t, exec): 
    '''true if a_i is connected to agents a_0... a_i-1 at time t 
    look in the list of neighbours of a_i at time t if there is a_j, j<i'''
    if i == 0 :
        return True
    neighbours = G_C.neighbors(exec[i][t], mode = "all")
    for j in range(0, i):
        if exec[A_ordered_id[j]][t] in neighbours:
            return True
    return False

def pick_time_with_conflict(exec, G_C) : 
    '''Choose t around the middle of the execution, with conflicts '''
    max_len = max(map(len, exec))
    for i in range(0,max_len//2):
        if nb_conflicts([[exec_i[max_len//2 + i] for exec_i in exec]], G_C) > 0 :
            return max_len//2 +i 
        elif nb_conflicts([[exec_i[max_len//2 - i] for exec_i in exec]], G_C) > 0 :
            return max_len//2 -i 
    return max_len//2

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



def execution_with_best_neighbour(G_M, G_C, sources, targets, A_ordered_id, i, t, exec):
    '''Choose a neighbour u of a_0...a_i-1 with path of length t from s_i to u, minimize d(u, g_i) and nb of conflicts 
    Output : execution with a_i going through u at t'''
    Neighbours = []
    for j in range(0,i):
        Neighbours+= G_C.neighbors(exec[A_ordered_id[j]][t], mode = "all")
    min_dist = len(exec)
    min_nb_conflicts = nb_conflicts(exec, G_C)
    best = Neighbours[0]
    best_exec = deepcopy(exec)
    for u in Neighbours:
        dist_u_gi = len(decoupled_exec(G_M, [u], [targets[i]])) #d(u, g_i)
        sources_first = [sources[j] for j in range(len(sources))]
        targets_first = [exec[j][t] for j in range(len(sources))]
        targets_first[i] = u
        sources_second = [exec[j][t] for j in range(len(sources))]
        targets_second = [targets[j] for j in range(len(sources))]
        sources_second[i] = u
        exec_tested = concatanate_executions(decoupled_exec(G_M, sources_first, targets_first),decoupled_exec(G_M, sources_second, targets_second))
        if dist_u_gi < min_dist :
            if nb_conflicts(exec_tested, G_C) < min_nb_conflicts:
                best = u
                min_dist = dist_u_gi
                min_nb_conflicts = nb_conflicts(exec_tested, G_C)
                best_exec = exec_tested
        #if there is one, then compute the distance d(u, g_i) and the nb of conflicts: 
        #if it's less than min_dist then and min_nb_conflicts then replace best by u
    return best_exec


###Algorithm (must return list of paths)

def mapfdivideandconquer(G_Mname, G_Cname, sources, targets):
    '''This algorithm's method is divide and conquer '''
    G_C = igraph.read(G_Cname)
    G_M = igraph.read(G_Mname)
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
    #print(exec)
    if nb_conflicts(exec, G_C) == 0 :
        return exec
    if n >0: #number of recursive calls = 5
        t = pick_time_with_conflict(exec, G_C)
        #print(t)
        for i in A_ordered_id:
            if not(is_ordered_connected(G_C, A_ordered_id, i, t, exec)):
                exec_changed = execution_with_best_neighbour(G_M,G_C, sources, targets, A_ordered_id, i, t, exec) #update of exec_i
                if nb_conflicts(exec_changed, G_C) < nb_conflicts(exec, G_C):
                    exec = exec_changed
        L1 =  aux_divide(sources, [exec_i[t] for exec_i in exec], A_ordered_id, G_C, G_M, n-1) 
        L2 = aux_divide([exec_i[t] for exec_i in exec], targets, A_ordered_id, G_C, G_M, n-1)
        return concatanate_executions(L1, L2)
    else :
        return exec

def concatanate_executions(ex1, ex2):
    '''Input : 2 executions with same number of agents
    Output : the concatenation (in time) of the executions '''
    ex_final = deepcopy(ex1)
    for i_a in range(len(ex2)):
        ex_final[i_a]+= ex2[i_a][1:]
    return ex_final

### Tests

'''G_comm = igraph.read("map1.png_comm_uniform_grid_1_range_6.graphml")
G_mov = igraph.read("map1.png_phys_uniform_grid_1_range_6.graphml")

sources = [47, 14]
targets = [55, 20]
sources2 = [18, 10]
targets2 = [58, 62]
exec2 = decoupled_exec(G_mov, sources2, targets2)

pred_ex = get_pred(G_mov, 47, 55)
path_ex = extract_path_from_pred(pred_ex, 47, 55)

exec_ex = decoupled_exec(G_mov, sources, targets)
#print(exec2)

#print(nb_conflicts(exec2, G_comm))
#print(G_comm.are_connected(18,10))

#print(concatanate_executions([[1,2,3], [4,5,6]], [[3,7,8],[6,9,10]]))

#print(G_mov.vs[3]["x_coord"])'''





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

