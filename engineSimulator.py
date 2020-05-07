from copy import deepcopy
from multiprocessing import Process, Queue
from CEA_DataGenerator import CEADataGenerator
import math

class engineSimulator():
    
    accentDecentAccel = 0
    Tb = 0
    m0 = 0
    Cd = 0
    A = 0
    rho0 = 1.225
    Patm = 101325
    T0 = 288.15
    ox = ""
    fuel = ""
    a = 0
    n = 0
    expansionHalfAngle = 0
    inefficiencyFactor = 0
    nozzleIneffiencyFactor = 0


    def setFlightPofile(self , accentDecentAccel = 5, Tb = 20):
        self.accentDecentAccel = accentDecentAccel
        self.Tb = Tb

    def setInitialRocketMass(self,m):
        self.m = m

    def setDragArea(self, A):
        self.A = A

    def setOX(self, oxName = "N2O"):
        self.ox = oxName

    def setFuel(self, fuelName = "paraffin"):
        self.fuel = fuelName
    
    def setFuelRegressionCoeffs(self,a,n):
        self.a = a
        self.n = n

    def setExpansionHalfAngel(self,a=0.261799):
        self.nozzleIneffiencyFactor(a)
        self.expansionHalfAngle = a

    def setInefficienryFactor(self,n = 0.94):
        self.inefficiencyFactor = n

    def setNozzleIneffiencyFactor(self,a):
        self.nozzleIneffiencyFactor = 1/2 * (1+math.cos(a))

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
    
    def getRho(self,y):
        return self.rho0

    def getReqThrust(self , t, D):
        return self.flightProfile(t) / self.m + D

    def getAb(self,L,r):
        return 2 * math.pi * r * L

    def getAc(self, r):
        return math.pi * (r ** 2)

    def getGox(self,moxdot,Ac):
        return moxdot / Ac

    def getRdot(self,moxdot,r):
        Ac = self.getAc(r)
        return self.a * (self.getGox(moxdot,Ac) ** self.n)
    
    def getMoxdot(self,mdot,OF):
        return mdot/(1+OF)
    
    def getOF(self,moxdot,Ac,rdot):
        return moxdot/(Ac*rdot)
    
    def getPcFromT(self,T,At,Cf):
        return T/(self.inefficiencyFactor * self.nozzleIneffiencyFactor * Cf * At)

    def getMdotFromPc(self, Pc , At , Cstr):
        return (Pc * At) / Cstr