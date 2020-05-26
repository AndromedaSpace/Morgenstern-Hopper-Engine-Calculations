import json
from engineSimulator import engineSimulator

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

    def cost(self,simResults):
        return -simResults['medianIsp'] + simResults['medianPc'] + simResults['medianTc'] + simResults['mprop'] 

    def rangeGenerator(self,data):
        temp = []
        step = data['step']
        cur = data['min']

        while cur <= data['max']:
            temp.append(cur)
            cur += step

        return temp

    def calculatePartialDerivative(self, currentState , derivativeTerm, confData):
        dv = (confData['startingStates'][derivativeTerm]['max'] - confData['startingStates'][derivativeTerm]['min']) * confData['optimiser']['d']
        
        posState = currentState
        posState[derivativeTerm] += dv / 2

        negState = currentState
        negState[derivativeTerm] -= dv / 2

        posCost = self.cost(self.engine.stateSimulationHandler(
            P0 = posState['P0'],
            OF0 = posState['OF0'],
            eps = posState['eps'],
            L = posState['L'],
            ThroatResizeCoeff = posState['ThroatResizeCoeff'],
            dt = confData['optimiser']['dt']
        ))

        negCost = self.cost(self.engine.stateSimulationHandler(
            P0 = negState['P0'],
            OF0 = negState['OF0'],
            eps = negState['eps'],
            L = negState['L'],
            ThroatResizeCoeff = negState['ThroatResizeCoeff'],
            dt = confData['optimiser']['dt']
        ))

        return (posCost - negCost) / dv

if __name__ == "__main__":
    optimiser = engineOptimiser()