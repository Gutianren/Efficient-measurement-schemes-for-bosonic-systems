import numpy as np
import numpy.random as nr


def dit(x,d,n):
    # change number x to Dinary
    # inversed order / ex: 5(10) -> 21(3)   Trinary: 5(10) -> 12(3)
    r = np.zeros(n,dtype = np.int32)
    for i in range(n):
        r[i] = x % d
        x = x // d
    return r

def adit(x,d,n):
    # Inverse operation of dit
    j = 0
    for i in range(n):
        j += x[i] * (d ** i)
    return j

def sympletic_product(x1,z1,x2,z2,d):
    # Calculate the sympletic product of two Pauli operators
    # P1P2 = w^sym P2P1
    sym = (x2 @ z1.T - x1 @ z2.T) % d
    return sym

def getinv(d):
    # Solve for x the equation ax = b (mod d)
    inv = np.zeros((d,d),dtype=int)  # ax = b (mod d), a != 0
    for a in range(1,d):
        for x in range(1,d):
            inv[a,(a*x)%d] = x
    return inv

class state(object):

    def __init__(self, n, d):
        self.s = np.eye(2*n+1,dtype=int)
        self.n = n
        self.d = d
        self.s[-1,-1] = 0
        self.rho_d = d%2
        self.inv = getinv(d)  # ax = b (mod d)
    
    def copy(self):
        newState = state(self.n,self.d)
        newState.s = self.s.copy()
        return newState

    def CNOT(self,a,b,k=1):
        # CNOT * X1^a1 Z1^a2 * X2^b1 Z2^b2 * CNOT'
        # = X1^a1 Z1^(a2-b2) * X2^(a1+b1) Z2^b2
        n = self.n 
        s = self.s
        d = self.d
        s[:-1,b] = (s[:-1,b] + s[:-1,a] * k) % d
        s[:-1,a+n] = (s[:-1,a+n] - s[:-1,b+n] * k) % d

    def F(self,a,k=1):  
        # Fourier transformation (~ Hadamard)
        # F X^aZ^b F' = w^(-ab) X^(-b) Z^a
        assert k == 1
        n = self.n
        s = self.s
        d = self.d
        s[:-1,-1] = (s[:-1,-1] - s[:-1,a] * s[:-1,a+n]) % d
        s[:-1,a], s[:-1,a+n] = (-s[:-1,a+n].copy()) % d, s[:-1,a].copy()
    
    def Fd(self,a,k=1): 
        # Inverse Fourier transformation
        # F' X^aZ^b F = w^(-ab) X^b Z^(-a)
        assert k == 1
        n = self.n
        s = self.s
        d = self.d
        s[:-1,-1] = (s[:-1,-1] - s[:-1,a] * s[:-1,a+n]) % d
        s[:-1,a], s[:-1,a+n] = s[:-1,a+n].copy(), (-s[:-1,a].copy()) % d

    def P(self,a,k=1):  
        # Phase Gate
        # P X^aZ^b P' = w^[a(a+rho_d)/2] X^a Z^(a+b)
        n = self.n
        s = self.s
        d = self.d
        s[:-1,-1] = (s[:-1,-1] + s[:-1,a] * (s[:-1,a] + self.rho_d) // 2 * k) % d
        s[:-1,a+n] = (s[:-1,a+n] + s[:-1,a] * k) % d

    def X(self,a,k=1):
        # X X^aZ^b X' = w^(-b) X^aZ^b
        n = self.n
        s = self.s
        d = self.d
        s[:-1,-1] = (s[:-1,-1] - s[:-1,a+n] * k) % d
    
    def Z(self,a,k=1):
        # Z X^aZ^b Z' = w^a X^aZ^b
        s = self.s
        d = self.d
        s[:-1,-1] = (s[:-1,-1] + s[:-1,a] * k) % d

    def g(self,a,b,c,d,m,n):
        # Extra phase when adding(multiplying) two XZ operators
        # (X^a Z^b)^m * (X^c Z^d)^n = w^g X^(ma+nc) Z^(mb+nd)
        D = self.d 
        factor = (m * (m-1) / 2 * a * b + n * (n-1) / 2 * c * d + m * n * b * c) % D
        return factor

    def eliminateX(self,h,i,a):
        # eliminate the X component at position a in row h using row i
        s = self.s
        assert s[i,a] != 0
        self.rowsum(h,i,s[i,a],-s[h,a])

    def rowsum(self,h,i,p,q):
        # multiplying Ri^q to Rh^p and save it to row h
        if  p == 1 and q == 0:
            return
        n = self.n
        s = self.s
        d = self.d
        for j in range(n):
            s[h,-1] = (s[h,-1] + self.g(s[h,j],s[h,j+n],s[i,j],s[i,j+n],p,q)) % d
        s[h,:] = (s[h,:] * p + s[i,:] * q) % d

    def sympletic_product(self, h, i):
        # calculate the sympletic product of row h and row i
        n = self.n
        s = self.s
        d = self.d
        return sympletic_product(s[h,:n],s[h,n:-1],\
                                 s[i,:n],s[i,n:-1],d)

    def measure(self, a):
        n = self.n 
        d = self.d 
        s = self.s
        p0 = np.where(s[n:-1,a]>0)[0] + n 
        if len(p0) > 0: # random case
            
            p = p0[0]
            # Using row p to eliminate the X components in other rows
            for i in range(2*n):
                if i != p and s[i,a] > 0:
                    self.eliminateX(i,p,a)
            
            # Copy row p to the destabilizer
            s[p-n] = s[p]
            s[p] = 0
            result = int(nr.rand() * d)
            
            # New stabilizer should be Z_a
            # The probability of measuring all terms equally distributes
            s[p,-1] = (-result) % d
            s[p,a+n] = 1
            return result

        else: # deterministic case

            s[-1] = 0
            s[-1,a+n] = 1

            # Calculate the coefficients c that
            # c1R1 + ... + cnRn = Za (Ris are the stabilizers)
            # as Za must in the stabilizer space
            c = np.zeros(n)
            for j in range(n):
                a = self.sympletic_product(j,j+n)
                b = self.sympletic_product(j,-1)
                c[j] = self.inv[a,b]
            
            # Summing up the rows according to c
            s[-1] = 0
            for j in range(n):
                self.rowsum(-1,j+n,1,c[j])
            
            return (-s[-1,-1]) % d

    def traceoperator(self,P):
        # Get trace(rho*P) = <phi | P | phi>
        # P is some n-mode d-dim Pauli operator
        # saved by 2n dits numpy array
        n = self.n
        d = self.d
        s = self.s
        w = np.cos(2*np.pi / d) + 1j * np.sin(2*np.pi / d)

        # If P does not commute with all the stabilizers
        # The result must be 0
        for i in range(n):
            sym = sympletic_product(s[i+n,:n],s[i+n,n:-1],P[:n],P[n:],d)
            if sym != 0:
                return 0
        
        # If P lies in the stabilizer space
        # we use the same method as the deterministic measurement
        s[-1,:2*n] = P
        s[-1,-1] = 0
        c = np.zeros(n)
        for j in range(n):
            a = self.sympletic_product(j,j+n)
            b = self.sympletic_product(j,-1)
            c[j] = self.inv[a,b]
        
        s[-1] = 0
        for j in range(n):
            self.rowsum(-1,j+n,1,c[j])
        
        return w ** ((-s[-1,-1]) % d)

    def stab(self):
        # print the stabilizer form of the state
        n = self.n
        s = self.s
        d = self.d
        ss = 'Stabilizer:\n'
        for i in range(n,2*n):
            for j in range(n):
                s1 = ''
                if s[i,j] != 0:
                    s1 += f'X{s[i,j]}'
                if s[i,j+n] != 0:
                    s1 += f'Z{s[i,j+n]}'
                if s[i,j] == 0 and s[i,j+n] == 0:
                    s1 = 'I'
                s1 += ' ' * (4-len(s1))
                ss += s1 + ' '
            ss += f'w{s[i,-1]}\n'
        ss += 'Destabilizer:\n'
        for i in range(n):
            for j in range(n):
                s1 = ''
                if s[i,j] != 0:
                    s1 += f'X{s[i,j]}'
                if s[i,j+n] != 0:
                    s1 += f'Z{s[i,j+n]}'
                if s[i,j] == 0 and s[i,j+n] == 0:
                    s1 = 'I'
                s1 += ' ' * (4-len(s1))
                ss += s1 + ' '
            ss += f'w{s[i,-1]}\n'
        return ss

    def _getstate(self, initial_c):
        n = self.n
        d = self.d
        s = self.s
        w = np.cos(2*np.pi / d) + 1j * np.sin(2*np.pi / d)
        #c = np.zeros(d**n,dtype=np.complex128)
        c = initial_c.copy()
        #c[12] = 1
        for i in range(n): # all stabilizers
            r = np.zeros(d**n,dtype = np.complex128)
            for k in range(d):  # S_i^k
                for j in range(d**n):
                    if np.abs(c[j]) > 1e-5:
                        x = dit(j,d,n)
                        fact = s[i+n,-1] * k
                        for l in range(n): # for all modes:
                            fact = (fact  \
                                    + k*(k-1)/2 * s[i+n,l] * s[i+n, l+n] \
                                    + s[i+n,l+n] * x[l] * k) % d
                            x[l] = (x[l] + k * s[i+n,l]) % d
                        j1 = adit(x,d,n)
                        r[j1] += c[j] * (w**fact)
            c = r.copy()
            # c /= np.sqrt(np.sum(c * c.conj()))
        c = np.reshape(c,(1,d**n))
        return c

    def getstate(self):  
        # Output the vector form of the state
        # Inversed order: c[0] = |00>, c[1] = |10>, ..., c[8] = |22>
        # output shape: (D,1) 
        # very slow
        # For test
        n = self.n
        d = self.d
        for i in range(d**n):
            initial_c = np.zeros(d**n)
            initial_c[i] = 1
            c = self._getstate(initial_c)
            norm = np.abs(np.sum(c*c.conj()))
            if norm > 1e-10:
                c = c / np.sqrt(norm)
                return c.reshape((d**n,1))

    def tostate(self):
        n = self.n
        d = self.d
        c = self.getstate().flatten()
        ss = ''
        for i in range(d**n):
            if np.abs(c[i]) > 1e-5:
                x = dit(i,d,n)
                ss += '{:f}|'.format(c[i])
                for j in range(n):
                    ss += str(x[j])
                ss += '>+\n'
        return ss                 

    def __str__(self):
        return str(self.s)

if __name__ == '__main__':

    d = 3
    n = 1
    s = state(n,d)
    #print(s.stab())
    #print(s.tostate())
    s.Fd(0)
    s.Z(0)
    print(s.tostate())
    print(s.stab()) 
    '''
    s.CNOT(1,2)
    print(s.tostate())
    print(s.stab())
    print(s.measure(1))
    print(s.stab())
    print(s.tostate())
    print(s.measure(1))
    print(s.tostate())
    '''
