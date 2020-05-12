from copy import deepcopy
from multiprocessing import Process, Queue
from CEA_DataGenerator import CEADataGenerator
import math
import fluids
from scipy.optimize import fsolve

class engineSimulator():
    
    accentDecentAccel = 0
    Tb = 0
    m0 = 0
    Cd = 0
    A = 0
    rho0 = 1.225
    Patm = 101325
    h0 = 0
    T0 = 288.15
    ox = ""
    fuel = ""
    fuelRho = 0
    a = 0
    n = 0
    expansionHalfAngle = 0
    inefficiencyFactor = 0
    nozzleIneffiencyFactor = 0


    def __init__ (self, y0 = 0, accentDecentAccel = 5 , Tb = 20 , m0 = 7, A = 0.565486678 , oxName = "N2O", fuelName = "paraffin", fuelRho = 924.0, a = 0.472, n = 0.555, expsHalf = 0.261799, Inef = 0.94):
        self.setInitialHeight(h=y0)
        self.setFlightPofile(accentDecentAccel=accentDecentAccel, Tb = Tb)
        self.setInitialRocketMass(m=m0)
        self.setDragArea(A=A)
        self.setOX(oxName=oxName)
        self.setFuel(fuelName=fuelName)
        self.setFuelRho(fuelRho = fuelRho)
        self.setFuelRegressionCoeffs(a=a,n=n)
        self.setExpansionHalfAngel(a=expsHalf)
        self.setInefficienryFactor(n=Inef)

    def setFlightPofile(self , accentDecentAccel, Tb):
        self.accentDecentAccel = accentDecentAccel
        self.Tb = Tb

    def setInitialRocketMass(self,m):
        self.m0 = m

    def setDragArea(self, A):
        self.A = A

    def setOX(self, oxName):
        self.ox = oxName

    def setFuel(self, fuelName):
        self.fuel = fuelName
    
    def setFuelRho(self ,fuelRho):
        self.fuelRho = fuelRho

    def setFuelRegressionCoeffs(self,a,n):
        self.a = a
        self.n = n

    def setExpansionHalfAngel(self,a):
        self.setNozzleIneffiencyFactor(a)
        self.expansionHalfAngle = a

    def setInefficienryFactor(self,n):
        self.inefficiencyFactor = n

    def setNozzleIneffiencyFactor(self,a):
        self.nozzleIneffiencyFactor = 1/2 * (1+math.cos(a))

    def setInitialHeight(self,h):
        self.h0 = h

    def flightProfile(self,t):
        if t <= self.Tb/5:
            return self.accentDecentAccel
        elif t <= 3/10 * self.Tb:
            return -20/self.Tb * self.accentDecentAccel * (t-1/5*self.Tb) + self.accentDecentAccel
        elif t <= 7/10 * self.Tb:
            return -1 * self.accentDecentAccel
        elif t <= 4/5 * self.Tb:
            return  20 / self.Tb * self.accentDecentAccel * (t- 7/10 * self.Tb) - self.accentDecentAccel
        elif t <= self.Tb:
            return self.accentDecentAccel
        return 0

    def getDrag(self,Vy,rho,Cd):
        return 0.5 * rho * (Vy ** 2) * self.A * Cd

    def getCd(self, Vy):
        return 0.052
    
    def getPe(self , h):
        return self.Patm

    def getG(self,h):
        return 9.81

    def getRho(self,y):
        return self.rho0

    def getReqThrust(self , m , D , accel  , g):
        return (accel + g) *  m + D

    def getAb(self,L,r):
        return 2 * math.pi * r * L

    def getAc(self, r):
        return math.pi * (r ** 2)

    def getGox(self,moxdot,Ac):
        return moxdot / Ac

    def getRdot(self,moxdot,r):
        Ac = self.getAc(r)
        return self.a * (self.getGox(moxdot,Ac) ** self.n)/1000
    
    def getMoxdot(self,mdot,OF):
        return mdot/(1+1/OF)
    
    def getOF(self,moxdot,L,r, rdot):
        Ab = self.getAb(L=L , r = r)
        return moxdot/(Ab * rdot)
    
    def getPcFromT(self,T,At,Cf):
        return T/(self.inefficiencyFactor * self.nozzleIneffiencyFactor * Cf * At)

    def getMdotFromPc(self, Pc , At , Cstr):
        return (Pc * At) / Cstr

    def updateZeroOrderIntegral(self, val , valdot, dt):
        return val + valdot * dt

    def updateFirstOrderIntergral(self,val,valdot,valdotdot,dt):
        return val + valdot * dt + 1/2 * valdotdot * (dt ** 2)
    
    def getAt0(self,F,Cf,Pc):
        return F / (self.inefficiencyFactor * self.nozzleIneffiencyFactor * Cf * Pc)

    def getr0(self,mfdot,moxdot,L,fuelRho):
        def funcToMin(r0 ,*data):
            a , n , L, mfueldot, moxdot, rho = data
            print("moxdot", moxdot ,"mfdot",mfueldot,"gox" , (moxdot) / ( math.pi * ( (r0) ** 2 ) ), "rdot in meters" , a * ( ( (moxdot) / ( math.pi * ( (r0) ** 2 ) ) )  ** n )/1000 , 'r0' , r0 , 'error' , ((2 * math.pi * r0 * rho * L * a * ( ( (moxdot ) / ( math.pi * ( (r0) ** 2 ) ) )  ** n )))/1000 - mfueldot ) 
            return ((2 * math.pi * r0 * rho * L * a * ( ( (moxdot) / ( math.pi * ( (r0) ** 2 ) ) )  ** n )))/1000 - mfueldot

        r0 = fsolve(funcToMin , 0.001 , (self.a , self.n , L , mfdot , moxdot , fuelRho) , factor=0.5)
        return r0

    def simInit(self, P0 , OF0, eps, L):
        generator = CEADataGenerator()
        generator.setFuel(fuelName = self.fuel)
        generator.setOX(oxName = self.ox)
        Ivac, Isp , Cstr, Tc, Cf, SeparationState = generator.singleWorker(Pe = self.getPe(self.h0), EPS=eps, P = P0 , OF = OF0)

        At = self.getAt0(F=self.getReqThrust(m = self.m0, D = 0 , accel = self.flightProfile(t=0), g = self.getG(self.h0) ) , Cf = Cf , Pc = P0)

        mdot = self.getMdotFromPc(Pc = P0 , At = At , Cstr = Cstr)
        
        moxdot = self.getMoxdot(mdot = mdot , OF = OF0) 

        mfueldot = mdot - moxdot

        r0 = self.getr0(mfdot = mfueldot , moxdot = moxdot , L = L , fuelRho = self.fuelRho)

        return At , r0

    def runEngineSim(self,Pc,eps,OF,At,r, L , h ,dt = 0.01):
        t = 0
        v = 0
        generator = CEADataGenerator()
        generator.setFuel(fuelName = self.fuel)
        generator.setOX(oxName = self.ox)
        medianIvac = 0
        medianIsp = 0
        medianTc = 0
        medianPc = 0

        m  = self.m0
        Tb = self.Tb
        h = self.h0

        iters  = int(Tb /dt)

        while t <= Tb:
            Pe = self.getPe(h = h)
            Ivac, Isp , Cstr, Tc, Cf, SeparationState = generator.singleWorker(Pe = Pe, EPS=eps, P = Pc , OF = OF)
            D = self.getDrag(Vy = v , rho = self.getRho(h), Cd = self.getCd(Vy = v))
            accel = self.flightProfile(t=t)
            T = self.getReqThrust(m = m , D = D,accel=accel,g = self.getG(h))
            Pc = self.getPcFromT(T = T ,At = At, Cf = Cf)
            mdot = self.getMdotFromPc(Pc = Pc , At = At , Cstr = Cstr)
            moxdot = self.getMoxdot(mdot = mdot , OF = OF)
            rdot = self.getRdot(moxdot = moxdot, r = r)
            r = self.updateZeroOrderIntegral(val = r , valdot = rdot, dt = dt)
            OF = self.getOF(moxdot = moxdot , L = L, r = r, rdot = rdot)
            h = self.updateFirstOrderIntergral(val = h , valdot = v , valdotdot = accel , dt = dt)
            v = self.updateZeroOrderIntegral(val = v , valdot = accel , dt = dt)
            m = self.updateZeroOrderIntegral(val = m , valdot = -mdot ,dt = dt)
            
            medianIvac += Ivac
            medianIsp += Isp
            medianPc += Pc
            medianTc += Tc
            
            if SeparationState['state'] == 'Separated':
                return{
                    'state' : 'Separation'  ,
                    'sepData' : SeparationState['data'] ,
                    't' : t ,
                    'Pc' : Pc
                }

            t += dt

        medianIvac /= iters
        medianIsp /= iters
        medianPc /= iters
        medianTc /= iters

        return {
            'state' : 'success',
            'medianIvac' : medianIvac,
            'medianIsp' : medianIsp,
            'medianPc' : medianPc,
            'medianTc' : medianTc
        }
        
if __name__ == '__main__':
    engine = engineSimulator(accentDecentAccel=200,m0=3,n=0.46,a=0.15)
    At , r0 = engine.simInit(30 * 10 ** 5, OF0 = 4, eps = 4, L = 0.33)
    print('At' ,At)
    print('Throat Radius' , math.sqrt(At / math.pi))
    print(r0)
