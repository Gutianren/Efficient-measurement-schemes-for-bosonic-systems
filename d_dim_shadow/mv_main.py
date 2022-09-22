import numpy as np
import numpy.random as nr
from mv_encoding import encode
from stabilizerstate import state
from clifford_d import main


def main_process(d,g,R,T,K):
    #d = 3   # must be a prime
    N_samples = T
    Repeat = R
    filename = 'H2OCoe0.txt'
    n, observablesP, ge, gs = encode(d,filename)
    alpha, P = observablesP
    M = len(P)
    print('Encoded observables:',M)
    #print('{}-dim Exact Ground Energy : {}'.format(d,ge))
    #print(dit(5,3,2))
    #state_after_measure = measure_state(gs, sample_num = N_samples)
    #print(state_after_measure)
    #print(alpha)
    #print(P)
    ghz = g
    st = state(n,d)
    #st.F(0)
    #for i in range(1,n):
    #    st.CNOT(0,i)
    
    real_pred = [st.traceoperator(p) for p in P]
    real_ans = 0
    for i in range(M):
        real_ans += alpha[i] * real_pred[i]
    print('Exact_result =',real_ans)

    #K = 10
    results = []
    err = []
    err_total = 0
    filename = 'zdata_ground_d{}R{}T{}K{}.txt'.format(d,Repeat,N_samples,K)
    with open(filename,'w') as f:
        for j in range(Repeat):
            pred = main(n,d,P,st,N_samples,N=N_samples//K,K=K)
            ans = 0
            for i in range(M):
                ans += alpha[i] * pred[i]
            err.append(np.abs(ans-real_ans))
            err_total += err[-1] ** 2
            print('total_variance=',err_total)
            print('Repeat {} of {}: {}, err = {}'.format(j, Repeat, ans, err[-1]))
            results.append(ans)
        results = np.array(results)
        err = np.array(err)

        sample_ans = np.average(results)
        deviation = np.sqrt(np.sum(err*err)/Repeat)
        print('Sample_result =',np.average(results))
        print('Exact_result =',real_ans)
        print('Err = ', deviation)
        print('Err_total=',err_total)
        
        print(M,file=f)
        for i in range(Repeat):
            print(results[i],err[i],file=f)
        print(sample_ans,file=f)
        print(real_ans,file=f)
        print(deviation,file=f)

if __name__ == '__main__':
    #d = 3
    g = 3
    R = 10
    #T = 10000
    K = 10
    
    for T in [100,200,500,2000,5000]:
        for d in [3,5,7]:
            main_process(d,g,R,T,K)
       