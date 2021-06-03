#from pysat.solvers import Solver
import numpy as np
import igraph

#2 graphs : movement & connection
#A : each agent is represented by a couple (current position, goal)
#P : list of configurations [[a0_t,...] for each t] (or list of paths)

def decoupled_exec(G_M, G_C, A) :
    '''algo shortest path for each agent : A*
    for now, we use a less efficient algorithm
    Output : execution P'''
    return 0

def nb_conflicts(P) :
    '''Use an algorithm already created or use A_ordered to see for each configuration if there is a disconnected agent 
    Outpub : int, number of conflicts'''
    return 0

def is_connected(G_C,A, A_ordered_id, i, t, P): 
    '''true if a_i is connected to agents a_0... a_i-1 at time t 
    look in the list of neighbours of a_i at time t if there is a_j, j<i'''
    neighbours = igraph.neighbors(G_C, P[t][i][0], mode = "all")
    for j in range(0, i):
        if P[t][A_ordered_id[j]][0] in neighbours:
            return True
    return False

def pick_time_with_conflict(P) : 
    '''Choose t around the middle of the execution, with conflicts '''
    for i in range(0,len(P)//2):
        if nb_conflicts(P[len(P)//2 + i]) > 0 :
            return len(P)//2 +i 
        elif nb_conflicts(P[len(P)//2 - i]) > 0 :
            return len(P)//2 -i 
    return len(P)//2

def choose_order(G_C, A) :
    '''Choose an order of agents
    Output : list of index '''
    i = np.random(len(A))
    A_ordered_id = [i]
    #BFS on A depending on the initial or final configuration
    queue = igraph.neighbors(G_C, A[i][0], mode = "all")
    while len(A_ordered_id)<len(A):
        v = queue.pop(0)
        for j in range(0, len(A)):
            if A[j][0]==v :
                A_ordered_id.append(j)
                queue += igraph.neighbors(G_C, A[j][0], mode = "all")
    return A_ordered_id

def choose_best_neighbour(G_C, A, A_ordered_id, i, t, P):
    '''Choose a neighbour u of a_0...a_i-1 with path of length t from s_i to u, minimize d(u, g_i) and nb of conflicts '''
    return 0 

def update(G_M, G_C, P, u, i, t):
    '''update P with a_i going through u at t'''
    A_first = [(P[0][j], P[t][j]) for j in range(0, len(P[0]))]
    A_first[i][1] = u
    A_second = [(P[t][j], P[len(P)][j]) for j in range(0, len(P[0]))]
    A_second[i][0] = u
    return decoupled_exec(G_M, G_C, A_first) + decoupled_exec(G_M, G_C, A_second)

def mapfdivideandconquer(G_M, G_C, A):
    P = decoupled_exec(G_M, G_C, A)
    nb_it = 0
    while nb_it < 5 : #number of attempts to find a better P
        A_ordered_id = choose_order(G_C, A) 
        P_changed = aux_divide(P,A, A_ordered_id, G_C, G_M, 0)
        if nb_conflicts(P_changed) == 0 :
            return P_changed
        else :
            nb_it+=1 #if P_changed has less conflicts, P is replaced by P_changed
    return 0

def aux_divide(P,A, A_ordered_id, G_C, G_M, n):
    if n < 5: #number of recursive calls = 5
        t = pick_time_with_conflict(P)
        P_changed = P
        for i in A_ordered_id:
            if not(is_connected(G_C,A, A_ordered_id, i, t, P_changed)):
                u = choose_best_neighbour(G_C, A, A_ordered_id, i, t, P_changed)
                update(G_M, G_C, P_changed, u, i, t) #update of P_i
            if nb_conflicts(P_changed) < nb_conflicts(P):
                P = P_changed
        return aux_divide(P[::len(P)//2],A, A_ordered_id, G_C, G_M, n+1) + aux_divide(P[len(P)//2::],A,A_ordered_id, G_C, G_M, n+1)
    else :
        return P


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

