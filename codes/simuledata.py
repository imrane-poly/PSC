import numpy as np
kbar = 3
display = True

def simuleData(star, param, M, v, s, T):
    lam, b, m0, beta, kappa, sigma, rho, thetaBar, alpha, b2Star = param
    d = 2**kbar
    rho2 = np.sqrt(1-(rho**2))
    A = transitionMatrix(star, m0, lam, b2Star, b)
    J = jumpsMatrix(lam, b, beta)
    JBar = JBarFun(J, A)
    
    dat = np.zeros(T)
    
    dat[0] = s    
    #s : s_{0}, v : v_{0} 
    
    W = np.random.randn(1) #correspond en fait à (rho*W1 + rho2*W2) (il n'y a pas de couplage à t = 0)
    
    theta = thetaMt(m0, M, thetaBar)
    
    if display:
        vtab = np.zeros(T)
        vtab[0] = v
        thetatab = np.zeros(T)
        thetatab[0] = theta

    
    v += kappa*(theta - v) + sigma*np.sqrt(v)*W    
    v = max(0, v) #ou faut-il laisser v négatif?

    
    #s : s{0}, v : v{1}
    
    for t in range(T-1):
        i = M #correspond à M_t

        M = np.random.choice(d, 1, True, A[i, :])[0]    #correspond à M_{t+1} 
                                            #(on choisit dans range(d) avec les probas données par la ligne de A)
        theta = thetaMt(m0, M, thetaBar)
        (W1, W2) = np.random.randn(2)
        
        #s : s_{t}, v : v_{t+1}
        
        mu = (alpha - 0.5)*v #mu : mu_{t+1}
        s += mu + np.sqrt(v)*W1 + J[i, M] - JBar[i]
        
        
        #s : s_{t+1}, v : v_{t+1}
        
        dat[t+1] = s
        
        if display:
            vtab[t+1] = v
            thetatab[t+1] = theta
        

        v += kappa*(theta - v) + sigma*np.sqrt(v)*(rho*W1 + rho2*W2)
        v = max(0, v) #ou faut-il laisser v négatif?
            
        #s : s_{t+1}, v : v_{t+2}
        
        dat[t+1] = s
    if display:
        return dat, vtab, thetatab
    else:
        return dat
