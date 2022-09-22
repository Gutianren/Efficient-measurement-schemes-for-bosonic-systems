import numpy as np
from sampleClifford import sampleClifford, getr
from stabilizerstate import state 

def lstinv(L,d):
    # Get the inverse list of L (U^dagger)
    L1 = []
    for u in L[::-1]:
        if u[0] == 'F':
            L1.append(('Fd',u[1],u[2]))
        else:
            L1.append((u[0],(-u[1]) % d, u[2]))
    return L1

def applyunitary(state, L):
    applydict = {
        'P': lambda x,k:state.P(x,k),
        'F': lambda x,k:state.F(x,k),
        'Fd': lambda x,k:state.Fd(x,k),
        'X': lambda x,k:state.X(x,k),
        'Z': lambda x,k:state.Z(x,k),
        'C': lambda x,k: state.CNOT(x[0],x[1],k)
    }
    j = 0
    for u in L:
        #for i in range(u[1]):
        applydict[u[0]](u[2],u[1])
        #print('--------{}--------'.format(j))
        #print(state.tostate())
        j += 1


def shadow(input_state, sample_num, r):
    # input_state is a stabilizer (state class)
    # sample_num is the number of classical shadows
    # output the classical shadows U^\dagger |b><b| U 
    d = input_state.d
    n = input_state.n
    shadows = [input_state.copy() for i in range(sample_num)]
    count = 0
    for s in shadows:
        LU = sampleClifford(d,n,r)
        #LU = LU[::-1]
        applyunitary(s,LU)

        for a in range(n):
            s.measure(a)

        LUd = lstinv(LU,d)
        applyunitary(s,LUd)
        count += 1
        #if count % (sample_num // 100) == 0:
        #    print('{}% shadow finished'.format(100 * count/sample_num))

    return shadows

def trace(P,n,d):
    assert len(P) == 2*n
    result = 1
    for i in range(n):
        if P[i] != 0 or P[i+n] != 0:
            return 0
        else:
            result *= d
    return result

def prediction(shadows, N, K, LP, d, n):
    # N * K must be equivalent with the number of samples
    # LP is the list of operators
    assert N * K == len(shadows)
    M = len(LP)
    result = []
    count = 0
    for P in LP:
        trs = np.array([(d**n+1) * s.traceoperator(P) - trace(P,n,d) for s in shadows])
        #trs = np.array([(d**n+1) * s.traceoperator(P) - traceP(np.eye(9),P,d,n) for s in shadows])
        #print(trs)
        aves = np.array([np.average(trs[N*i:N*(i+1)]) for i in range(K)])
        result.append(np.median(aves))
        count += 1
        #if M >= 100:
        #    if count % (M // 100) == 0:
        #        print('{}% prediction finished'.format(100 * count/M))
        #else:
        #    print('{}% prediction finished'.format(100 * count/M))
            
    return result

def shadowtest(shadows, N_samples, d, n):
    arr = np.zeros((d**n,d**n))
    for s in shadows:
        st = s.getstate()
        st = (d**n+1) * (st @ st.conj().T)- np.eye(d**n)
        arr = arr + st
    return arr / N_samples


def main(n, d, observables, initial_state, N_samples, N, K):
    assert N * K == N_samples
    assert n == initial_state.n 
    assert d == initial_state.d
    r = getr(d)
    Sh = shadow(initial_state, N_samples, r)
    results = prediction(Sh, N, K, observables, d, n)
    return results 


def traceP(rho,P,d,n):
    # P is a 2*n opearator
    assert d == 3
    assert n == 2
    r3 = np.sqrt(3)
    w = -1/2 + 1j * r3 / 2
    w2 = -1/2 - 1j * r3 / 2
    I3 = np.eye(3)
    X = np.array([[0,0,1],[1,0,0],[0,1,0]])
    Z = np.array([[1,0,0],[0,w,0],[0,0,w2]])

    S0 = I3.copy()
    for i in range(P[0]):
        S0 = S0 @ X
    for i in range(P[2]):
        S0 = S0 @ Z

    S1 = I3.copy()
    for i in range(P[1]):
        S1 = S1 @ X
    for i in range(P[3]):
        S1 = S1 @ Z
    
    SP = np.kron(S0,S1)
    return np.trace(SP @ rho)

if __name__ == '__main__':

    n = 2
    d = 5
    r = getr(d)

    np.set_printoptions(precision=2,suppress=True,linewidth=np.inf)

    s = state(n,d)
    s.F(0)
    L = sampleClifford(d,n,r)
    #applyunitary(s,L)
    print(s.tostate())
    st = s.getstate()
    st = st @ st.conj().T

    shadows = shadow(s,10000,r)
    print(st)
    print('-'*30)
    st1 = shadowtest(shadows,10000,d,n)
    print(st1)
    print(np.max(np.abs(st-st1)))

    '''
    P = np.array([2,0,0,1])
    P1 = np.array([0,2,1,0])
    ans_shadow = prediction(shadows,10000,1,[P],d,n)
    ans_exact = s.traceoperator(P)
    ans_pred = traceP(st,P1,d,n)
    ans_shadow2 = traceP(st1,P1,d,n)
    print(ans_shadow)
    print(ans_exact)
    print(ans_pred)
    print(ans_shadow2)
    '''


    '''
    P = np.zeros(2*n)   # Z1
    P[0] = 2
    P[3] = 1
    print(s.traceoperator(P))   # Exact result

    time_test = 0
    
    if time_test == 1:
        profile = LineProfiler(shadow)
        profile.runcall(shadow,s,10000,r)
        profile.print_stats()

    else:
        results = main(n,d,[P],s,10000,10000,1)
        print(results)
    '''