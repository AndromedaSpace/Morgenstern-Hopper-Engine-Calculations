import engineSimulator
import scipy
from scipy.optimize import minimize

def score(x,*args):
    print(x)
    P0 = x[0]
    OF0 = x[1]
    eps = x[2]
    L = x[3]
    Pmax , Tmax, Mmax , dt = args
    engine = engineSimulator.engineSimulator()
    engine.setMechanicalLimits(Pmax, Tmax, Mmax)
    res = engine.stateSimulationHandler(P0,OF0,eps,L,dt,breakAtFailure=False)
    print(res['medianIsp'])
    return 1/res['medianIsp']

res = minimize(score,[20*10**5 , 7 , 2 , 0.2], args=(30*10**5,7000,4,0.1),method='Nelder-Mead')
