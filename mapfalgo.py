#from pysat.solvers import Solver
from igraph import *
import numpy as np

#2 graphs : movement & connection
#A : each agent is represented by a couple (current position, goal)
#P : list of configurations [[a0_t,...] for each t] (or list of paths)

def decoupled_exec(G_M, G_C, A) :
    #algo shortest path for each agent : A*
    return 0

def nb_conflicts(P) :
    return 0

def is_connected(G_C, A, A_ordered, i, t, P): #a_i is connected to agents a_0... a_i-1 at time t 
    #look in the list of neighbours of a_i at time t if there is a_j, j<i
    return True

def pick_time_with_conflict(P) : #around the middle
    return len(P)//2

def choose_order(A) : #return index list
    i = np.random(len(A))
    A_ordered = [i]
    #BFS on A depending on the initial or final configuration
    return A_ordered

def choose_best_neighbour(G_C, A, A_ordered, i, t, P):
    return 0 #neighbour u of a_0...a_i-1 with path of length t from s_i to u, minimize d(u, g_i) and nb of conflicts 

def update(P, u, i):
    return P #update P with a_i going through u

def mapfdivideandconquer(G_M, G_C, A):
    P = decoupled_exec(G_M, G_C, A)
    nb_it = 0
    while nb_it < 5 : #number of try to find a better P
        A_ordered = choose_order(A) 
        P_changed = aux_divide(P,A, A_ordered, G_C, G_M, 0)
        if nb_conflicts(P_changed) == 0 :
            return P_changed
        else :
            nb_it+=1 #or if P_changed has less conflicts, replace P by P_changed
    return 0

def aux_divide(P,A, A_ordered, G_C, G_M, n):
    if n < 5: #number of recursive calls = 5
        t = pick_time_with_conflict(P)
        P_changed = P
        for i in range(1, len(A)):
            if not(is_connected(G_C,A, A_ordered, i, t, P_changed)):
                u = choose_best_neighbour(G_C, A, A_ordered, i, t, P_changed)
                update(P_changed, u, i) #update of P_i
            if nb_conflicts(P_changed) < nb_conflicts(P):
                P = P_changed
        return aux_divide(P[::len(P)//2],A, A_ordered, G_C, G_M, n+1) + aux_divide(P[len(P)//2::],A,A_ordered, G_C, G_M, n+1)
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

