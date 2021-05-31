#from pysat.solvers import Solver
from igraph import *
import numpy as np

#2 graphes : mouvement et connexion
#A : chaque agent est représenté par un couple (position, but)
#P : pour chaque t, donne la position de chaque agent ? ou renvoie un A_t

def decoupled_exec(G_M, G_C, A) :
    return 0

def nb_conflicts(P) :
    return 0

def is_connected(G_C, A, i, t, P): #a_i est connecté aux agents a_0... a_i-1 au temps t 
    return True

def pick_time_with_conflict(P) :
    return len(P)//2

def choose_order(A) :
    i = np.random(len(A))
    A_ordered = [i]
    #BFS sur A selon la config initiale ou finale
    return A_ordered

def choose_best_neighbour(G_C, A, i, t, P):
    return 0 #voisin des a_0...a_i-1 tq chemin de longueur t de s_i à u, minimise d(u, g_i) et le nb de conflits 

def update(P, u, i):
    return P #update P avec a_i qui passe par u

def mapfdichotomie(G_M, G_C, A):
    P = decoupled_exec(G_M, G_C, A)
    t = pick_time_with_conflict(P)
    nb_it = 0
    P_changed = []
    while nb_it < 5 :
        A_ordered = choose_order(A)
        P_changed[nb_it] = P
        for i in range(1, len(A)):
            if not(is_connected(G_C, A_ordered, i, t, P_changed[nb_it])):
                u = choose_best_neighbour(G_C, A, i, t, P_changed[nb_it])
                update(P_changed[nb_it], u, i) #màj de P_i
            if nb_conflicts(P_changed[nb_it]) < nb_conflicts(P):
                P = P_changed[nb_it]
        nb_it += 1

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

