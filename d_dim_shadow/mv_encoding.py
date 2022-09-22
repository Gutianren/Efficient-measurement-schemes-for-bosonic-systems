import numpy as np
import numpy.linalg as LA
from clifford_d import shadow, prediction, main
from stabilizerstate import dit,adit 

def decompose_verification(d,n,A,a):
    # Test whether the decomposed form equals the original matrix
    X = np.zeros((d,d))
    for i in range(d-1):
        X[i+1,i] = 1
    X[0,d-1] = 1
    Z = np.zeros((d,d),dtype=np.complex128)
    w = np.exp(1j * 2*np.pi / d)
    for i in range(d):
        Z[i,i] = w ** i

    S0 = np.zeros((d**n,d**n),dtype=np.complex128)
    for i in range(d**n):
        for j in range(d**n):
            x = dit(i,d,n)[::-1]
            z = dit(j,d,n)[::-1]
            S = np.eye(1,dtype=np.complex128)
            for k in range(n):
                s = np.eye(d,dtype=np.complex128)
                for l in range(x[k]):
                    s = s @ X
                for l in range(z[k]):
                    s = s @ Z 
                S = np.kron(S,s)
            S0 += A[i,j] * S
    
    return np.max(np.abs(S0-a))#,np.max(np.abs(S0)),np.max(np.abs(a))

def getdecompose(d,a):
    # decompose matrix a into stabilizers
    X = np.zeros((d,d))
    for i in range(d-1):
        X[i+1,i] = 1
    X[0,d-1] = 1
    Xd = X.conj().T 
    Z = np.zeros((d,d),dtype=np.complex128)
    w = np.exp(1j * 2*np.pi / d)
    for i in range(d):
        Z[i,i] = w ** i
    Zd = Z.conj().T

    A = np.zeros((d,d),dtype = np.complex128)
    for i in range(d):
        for j in range(d):
            S = np.eye(d)
            for k in range(i):
                S = S @ Xd
            for k in range(j):
                S = S @ Zd
            A[i,j] = 1/d * (w**(i*j)) * np.trace(S @ a)

    # Verification
    err = decompose_verification(d,1,A,a)
    assert err < 1e-10
    return A


def geta(d):
    a = np.zeros((d,d))
    for i in range(d-1):
        a[i,i+1] = np.sqrt(i+1)
    ad = a.conj().T
    A = getdecompose(d,a)
    Ad = getdecompose(d,ad)
    # A and Ad are the decomposing coefficients of a and ad
    # the latter 2 terms are the real a and ad matrix
    return A,Ad,a,ad

def getxp(d,M):
    A,Ad,a,ad = geta(d)
    x0 = ad+a
    p0 = ad-a
    x = np.eye(d)
    p = np.eye(d)
    X = [getdecompose(d,x)]
    P = [getdecompose(d,p)]
    for i in range(M):
        x = x @ x0
        p = p @ p0
        X.append(getdecompose(d,x))
        P.append(getdecompose(d,p))
    return x0,p0,X,P



def string_from_product(n,LA):
    # decompose a n-mode product form to Pauli string
    # the decomposed matrixes of n modes are stored in LA
    # return a (d**n) * (d**n) coefficient list
    # eps = 1e-12
    L = np.eye(1,dtype = np.complex128)
    for i in range(n):
        L = np.kron(L,LA[i])
    return L

def encoding(d,n,k,V):
    # k,V are lists, denoting ki = k_[Vi1,Vi2,...]
    # default: idx 0,1,...,d-1 has V = (1,1), (2,2), ..., (d,d)
    # ansreal is the real Hamiltonian
    # ansL returns the decomposition of the Hamiltonian into stabilizers

    assert len(k) == len(V)
    w = np.zeros(n)
    x,p,XX,PP = getxp(d,4)
    
    I0 = np.zeros((d,d))
    I0[0,0] = 1

    ansL = np.zeros((d**n,d**n),dtype=np.complex128)
    ansreal = np.zeros((d**n,d**n),dtype=np.complex128)
    for i in range(n):
        assert len(V[i]) == 2 and V[i][0] == i and V[i][1] == i
        w[i] = np.sqrt(2*k[i])

    # potential
    for i in range(len(k)):

        hp = np.zeros(n,dtype=int)
        hc = np.ones(n) # coefficient
        hreal = [np.eye(d) for j in range(n)]
        temp = np.eye(1,dtype=np.complex128)
        
        v = V[i]
        for j in range(len(v)):
            hp[v[j]] = hp[v[j]] + 1
            hc[v[j]] = hc[v[j]] / np.sqrt(2*w[v[j]])
            hreal[v[j]] = hreal[v[j]] @ x / np.sqrt(2*w[v[j]])

        h_product = []
        for j in range(n):
            h_product.append(hc[j] * XX[hp[j]])

        for j in range(n):
            temp = np.kron(temp, hreal[j]) 
        
        ansreal += temp * k[i]
        ansL += string_from_product(n,h_product) * k[i]    
    
    # kinetic energy
    
    h_product = [I0.copy() for j in range(n)]
    hreal = [np.eye(d) for j in range(n)]

    for i in range(n):
        h_product[i] = - PP[2] * w[i] / 2 / 2
        hreal[i] = - p @ p * w[i] / 2 / 2
        temp = np.eye(1,dtype=np.complex128)

        for j in range(n):
            temp = np.kron(temp, hreal[j])
        
        ansreal += temp
        ansL += string_from_product(n,h_product)
        h_product[i] = I0.copy()
        hreal[i] = np.eye(d)
    
    return ansL, ansreal

def readfile(filename):
    with open(filename,'r') as f:
        data = f.readlines()
    k = []
    V = []
    n = 0
    for da in data:
        dsep = [x.strip() for x in da.split(',')]
        k.append(float(dsep[-1]))
        v = [int(x)-1 for x in dsep[:-1]]
        V.append(v)
        n = max(n,np.max(v)+1)
    return n,k,V

def encode(d,filename):
    n,k,V = readfile(filename)
    ansL,ansreal = encoding(d,n,k,V)
    w,v = LA.eig(ansreal)
    idx = np.argsort(w)
    w = w[idx]
    v = v[:,idx]
    ground_energy = w[0]
    ground_state = v[:,0]

    P = []
    alpha = []
    for i in range(d**n):
        for j in range(d**n):
            if np.abs(ansL[i,j]) > 1e-12:
                P.append(np.concatenate((dit(i,d,n)[::-1], dit(j,d,n)[::-1])))
                alpha.append(ansL[i][j])

    return n, (alpha,P), ground_energy, ground_state


if __name__ == '__main__':
    d = 3
    n,k,V = readfile('H2OCoe0.txt')
    ansL,ansreal = encoding(d,n,k,V)
    print('Encoding Error :', decompose_verification(d,n,ansL,ansreal))
    count = 0
    for i in range(d**n):
        for j in range(d**n):
            if np.abs(ansL[i,j]) > 1e-10:
                count += 1
    w,v = LA.eig(ansreal)
    #print(np.min(np.abs(w)))
    #print(count)
    n,P,ge,gs = encode(3,'H2OCoe0.txt')
    print(ge)
