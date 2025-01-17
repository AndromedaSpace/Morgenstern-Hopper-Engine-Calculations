import json
from engineSimulator import engineSimulator
from scipy.optimize import minimize

class engineOptimiser:
    confData = None
    engine = None

    def __init__ (self,configFile="optimiser.conf"):
        self.confData = self.readData(filename=configFile)
        self.engine = self.initEngine(self.confData['engine'])


    def readData(self, filename):
        inFile = open(filename , 'r')
        jsonString = inFile.read()
        inFile.close()
        return json.loads(jsonString)

    def initEngine(self , data):
        engine = engineSimulator(
            y0 = data['y0'],
            accentDecentAccel = data['accentDecentAccel'],
            Tb = data['Tb'],
            Ti = data['Ti'],
            m0 = data['m0'],
            A = data['A'],
            oxName = data['oxName'],
            fuelName = data['fuelName'],
            fuelRho = data['fuelRho'],
            a = data['a'],
            n = data['n'],
            expsHalf = data['expsHalf'],
            Inef = data['Inef']
        )
        
        engine.setMechanicalLimits(
            Pmax = data['Pmax'],
            Tmax = data['Tmax'],
            Mmax = data['Mmax'],
            PtT = data['PtT']
        )

        return engine

    def score(self,x,*args):
        P0 = x[0]
        OF0 = x[1]
        eps = x[2]
        L = x[3]
        ThroatResizeCoeff = x[4]
        CIsp, CP, CT, Cm, dt, engine = args
        simResults = engine.stateSimulationHandler(P0,OF0,eps,L,ThroatResizeCoeff,dt=dt,breakAtFailure=False)
        print(CIsp * simResults['medianIsp'] + CP * simResults['medianPc'] + CT * simResults['medianTc'] + Cm * simResults['mprop'] )
        return - CIsp * simResults['medianIsp'] + CP * simResults['medianPc'] + CT * simResults['medianTc'] + Cm * simResults['mprop'] 

    def rangeGenerator(self,data):
        temp = []
        step = data['step']
        cur = data['min']

        while cur <= data['max']:
            temp.append(cur)
            cur += step

        return temp

    def optimise(self):
        startingPositions = {
            'P0' : self.rangeGenerator(self.confData['startingStates']['P0']),
            'OF0' : self.rangeGenerator(self.confData['startingStates']['OF0']),
            'eps' : self.rangeGenerator(self.confData['startingStates']['eps']),
            'L' : self.rangeGenerator(self.confData['startingStates']['L']),
            'ThroatResizeCoeff' : self.rangeGenerator(self.confData['startingStates']['ThroatResizeCoeff']),
        }
        bounds = (
            (self.confData['startingStates']['P0']['min'],self.confData['engine']['Pmax']),
            (0.01,None),
            (0.01,None),
            (0.01,None),
            (0.01,1)
        )
        localMinima = []
        for P0 in startingPositions['P0']:
            for OF0 in startingPositions['OF0']:
                for L in startingPositions['L']:
                    for eps in startingPositions['eps']:
                        for ThroatResizeCoeff in startingPositions['ThroatResizeCoeff']:
                            print([P0,OF0,eps,L,ThroatResizeCoeff])
                            localMinima.append(minimize(
                                self.score,
                                [P0,OF0,eps,L,ThroatResizeCoeff],
                                (self.confData['optimiser']['CIsp'], self.confData['optimiser']['CP'], self.confData['optimiser']['CT'],self.confData['optimiser']['Cm'],self.confData['optimiser']['dt'],self.engine),
                                method='TNC',
                                bounds=bounds 
                                ))
                            print(localMinima[-1])

        


if __name__ == "__main__":
    optimiser = engineOptimiser()
    optimiser.optimise()