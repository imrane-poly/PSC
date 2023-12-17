from random import gauss
import numpy as np

def simulategarchdata(T, omega, alpha1,alpha2, beta1, beta2):
    n = T
    series = [gauss(0,1), gauss(0,1)]
    vols = [1, 1]

    for _ in range(n):
        new_vol = np.sqrt(omega + alpha1*series[-1]**2 + alpha2*series[-2]**2 + beta1*vols[-1]**2 + beta2*vols[-2]**2)
        new_val = gauss(0,1) * new_vol

        vols.append(new_vol)
        series.append(new_val)
    return [vols,series]

