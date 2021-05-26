import pysat
from pysat.solvers import Solver

NMAX = 800

def decoupled_exec(V, Ec, Em, A) :
    return [[] for a in A]

def has_conflict(P) :
    return False

def is_connected(V, Ec, Em, A, agent, t):
    return True

def pick_disconnected(V, Ec, Em, A):
    return A[0], 0

def find_path(V, Ec, Em, A, n, constraints):
    return Solver()

def extract_exec(V, Ec, Em, A, constraints):
    return [[] for a in A]

def pick_config(V, Ec, Em, A, t):
    return [[0] for a in A]

def MAPF1(V, Ec, Em, A) :
    constraints = []
    P = decoupled_exec(V, Ec, Em, A)
    n = max([len(Pai) for Pai in P])
    while has_conflict(P):
        agent, t = pick_disconnected(V, Ec, Em, A)
        constraints.append(is_connected(V, Ec, Em, A, agent, t))
        while(find_path(V, Ec, Em, A, n, constraints).solve() == False) :
            n+= 1
        P = extract_exec(V, Ec, Em, A, constraints)
    return P

def MAPFmalin(V, Ec, Em, A, final_state):
    constraints = []
    P = decoupled_exec(V, Ec, Em, A)
    n = max([len(Pai) for Pai in P])
    if not(has_conflict(P)) :
        return P
    else :
        agent, t = pick_disconnected(V, Ec, Em, A)
        constraints.append(is_connected(V, Ec, Em, A, agent, t))
        while n < NMAX :
            while(find_path(V, Ec, Em, A, n, constraints).solve() == False) :
                n+= 1
            P_extr = extract_exec(V, Ec, Em, A, n, constraints)
            P_final = MAPFmalin(V, Ec, Em, A, P_extr[-1]) + MAPFmalin(V, Ec, Em, P_extr[-1], final_state)
            if not(has_conflict(P_final)) :
                return P_final
            else :
                constraints.append(pick_config(V, Ec, Em, A, t) != P_extr[t])


