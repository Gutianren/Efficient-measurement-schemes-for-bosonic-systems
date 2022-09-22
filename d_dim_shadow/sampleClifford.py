import numpy as np
import numpy.random as nr
from stabilizerstate import dit,sympletic_product,getinv

def getr(d):
    r = np.zeros((d**2,d**2))
    for a in range(1,d**2):
        for b in range(d**2):
            x = dit(a,d,2)
            y = dit(b,d,2)
            r[a,b] = (x[1] * y[0] - x[0] * y[1]) % d
    return r
    
def addF(a,b,d,k,j,m,L,sa,sb):
    assert k == 1
    L.append(('F',k,j))
    a[j], a[j+m] = (-a[j+m]) % d, a[j]
    b[j], b[j+m] = (-b[j+m]) % d, b[j]
    sa1 = (sa + a[j]*a[j+m]) % d 
    sb1 = (sb + b[j]*b[j+m]) % d
    return sa1, sb1

def addP(a,b,d,k,j,m,L,sa,sb):
    rho_d = d%2
    L.append(('P',k,j))
    a[j+m] = (a[j+m] + a[j] * k) % d
    b[j+m] = (b[j+m] + b[j] * k) % d
    sa1 = (sa + a[j] * (a[j] + rho_d) // 2 * k) % d
    sb1 = (sb + b[j] * (b[j] + rho_d) // 2 * k) % d
    return sa1,sb1

def addC(a,b,d,k,i,j,m,L,sa,sb):

    L.append(('C',k,(i,j)))
    a[j] = (a[j] + k * a[i]) % d
    a[i+m] = (a[i+m] - k * a[j+m]) % d
    b[j] = (b[j] + k * b[i]) % d
    b[i+m] = (b[i+m] - k * b[j+m]) % d
    return sa,sb

def sampleClifford(d,n,r):
    # sample an d-dimension n-bit Clifford operator
    # output a list with gates/powers and qudit numbers
    # Refer to arxiv:2008.06011

    inv = getinv(d)
    #rho_d = d%2

    L0 = []

    for l in range(n):
        m = n - l # Sample length
        
        # sample the first row
        ra = nr.randint(1,d**(2*m))
        a = dit(ra,d,2*m)
        sa = nr.randint(0,d)   # sign of the first row

        # sample the second row
        j = np.where(a>0)[0][0]
        if j >= m:
            j = j - m
        
        rb = nr.randint(d**(2*m))
        b = dit(rb,d,2*m)
        
        # b[1-j] = rb
        # b[j] = 0            # sign of the second row
        sb = nr.randint(0,d)

        # confirm that ab = w^{-1} * ba (commutation relations)
        # global phase will not affect
        # b[j] = 0
        # b[j+m] = 0

        sym = sympletic_product(a[:m],a[m:],b[:m],b[m:],d)
        symbj = (a[j+m] * b[j] - a[j] * b[j+m]) % d
        #print('sym = ',sym)
    
        needsym = (-1 - (sym - symbj)) % d

        # Now a[j],a[j+m] is fixed
        # We want to rearrange b[j],b[j+m] to 
        # make the sympletic product be -1
        # This arrangement should be equal.
        # needsym = (-1-sym) % d

        x_a = a[j] + d * a[j+m]
        x_b = b[j] + d * b[j+m]
        idx_need = np.where(r[x_a]==needsym)[0]
        idx_equals_symbj = np.where(r[x_a] == symbj)[0]
        try:
            newidx = np.where(idx_equals_symbj == x_b)[0][0]
        except:
            assert 0
        new_x_b = idx_need[newidx]
        b[j] = new_x_b % d
        b[j+m] = ((new_x_b) // d) % d

        # print('After fixing commutation:')
        # print('a = ', a)
        # print('b = ', b)

        # Sweeping Process
        # Change a and b to X_l^kx and Z_l^kz
        # The unitarys guarantee the commutation relations
        # Get the gate list

        L = []  # Save the gates/powers/index (temp)

        # Clear the Z components of row a with Fourier(F) and Phase(P):
        # F(Z^b) = X^(-b)
        # F(X^a) = Z^a
        # P(X^aZ^b) = w^(1/2a(a+rho)) * X^a Z^(a+b)
        for i in range(m): 

            if a[i+m] > 0:
                
                if a[i] == 0: # X^0, apply F for 1 time
                    sa,sb = addF(a,b,d,1,i,m,L,sa,sb)

                else:   # X^a, apply P for k times
                    k = inv[a[i],(-a[i+m]) % d]
                    sa,sb = addP(a,b,d,k,i,m,L,sa,sb)
        
        # Clear the X components of row a with CNOT
        # CNOT(X^a,X^b) = X^a, X^(a+b)
        J = np.where(a[:m]>0)[0]
        while len(J) > 1:
            for i in range(len(J)):
                if i & 1 == 0 and i+1 < len(J):
                    k = inv[a[J[i]], (-a[J[i+1]]) % d]
                    sa,sb = addC(a,b,d,k,J[i],J[i+1],m,L,sa,sb)
            J = np.where(a[:m] > 0)[0]

        # Now we have only one X^k in row a
        # Change this X^k to qudit 0 if it lies at some other position
        # We do this by using 2 times of CNOT^k

        if J[0] != 0:

            sa,sb = addC(a,b,d,1,J[0],0,m,L,sa,sb)

            k = inv[a[J[0]],(-a[J[0]]) % d]
            sa,sb = addC(a,b,d,k,0,J[0],m,L,sa,sb)

        
        # Clear Row b
        J1 = np.where(b[:2*m]>0)[0]
        if not(len(J1) == 1 and J1[0] == m): # If b is not 'ZI...I'
            
            sa,sb = addF(a,b,d,1,0,m,L,sa,sb)

            # Clear the Z components of row b
            for i in range(m):
                if b[i+m] > 0:
                    if b[i] == 0: # X^0, apply F for 1 time
                        sa,sb = addF(a,b,d,1,i,m,L,sa,sb)

                    else:   # X^a, apply P for k times
                        k = inv[b[i],(-b[i+m]) % d]
                        sa,sb = addP(a,b,d,k,i,m,L,sa,sb)
            

            # Clear the X components of row b
            # Due to the commutation relation we know the remaining element
            # in set J must be 1, so we needn't swap the qudits again
            J = np.where(b[:m]>0)[0]
            while len(J) > 1:
                for i in range(len(J)):
                    if i & 1 == 0 and i+1 < len(J):
                        k = inv[b[J[i]], (-b[J[i+1]]) % d]
                        sa,sb = addC(a,b,d,k,J[i],J[i+1],m,L,sa,sb)

                J = np.where(b[:m] > 0)[0]
            
            sa,sb = addF(a,b,d,1,0,m,L,sa,sb)

        # Now there are only X^a and Z^b left
        # We use P and F to eliminate them
        # PX^aP' = w^(1/2a(a+rho)) X^aZ^a
        # PZ^bP' = Z^b
        # F(X^aZ^b) = w^(ab)X^(-b)Z^a
        
        try:
            assert a[m] == 0
            assert b[0] == 0 
        except:
            pass
            assert 0

        if a[0] != 1:
            
            aa, bb = a[0], b[m] # a and b
            
            # X^a / Z^b
            sa,sb = addF(a,b,d,1,0,m,L,sa,sb)
            # Z^a / X^(-b)  no extra phase

            k = inv[d-bb,d-1]
            sa,sb = addP(a,b,d,k,0,m,L,sa,sb)
            # Z^a / X^(-b)Z^(-1)

            sa,sb = addF(a,b,d,1,0,m,L,sa,sb)
            # X^(-a) / XZ^(-b)

            sa,sb = addP(a,b,d,bb,0,m,L,sa,sb)
            # X^(-a) Z^(-1) / X

            sa,sb = addF(a,b,d,1,0,m,L,sa,sb)
            # XZ^(-a) / Z

            sa,sb = addP(a,b,d,aa,0,m,L,sa,sb)
            # X / Z

        # print('Final:')
        # print('a = ', a)
        # print('b = ', b)

        # Clear the signs with X and Z
        # X(XaZb)X^(-1) = w^(-b) XaZb
        # Z(XaZb)Z^(-1) = w^a XaZb
        # so we use ZXZ^(-1) = wX and XZX^(-1) = w^(-1) Z
        
        if sa != 0: # X
            L.append(('X',(-sa)%d,0))
            sa = 0
        
        if sb != 0:
            L.append(('Z',sb % d,0))
            sb = 0

        # Save all this operations into the final list L0
        # Note that each index should be increased by l
        try:
            for i in range(len(L)):
                if L[i][0] != 'C':
                    L0.append((L[i][0],L[i][1],L[i][2]+l))
                else:
                    L0.append((L[i][0],L[i][1],(L[i][2][0]+l,L[i][2][1]+l)))
        except:
            print('Error Here')
    # print('L0 = ', L0) 
    return L0


    

if __name__ == '__main__':
    pass
