import json
from engineSimulator import engineSimulator

class engineOptimiser:
    confData = None
    engine = None

    def __init__ (self,configFile="optimiser.conf", nProcesses = 12):
        self.confData = self.readData(filename=configFile)
        self.engine = self.initEngine(self.confData['engine'])
        self.nProcesses = nProcesses


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

    def score(self,simResults):
        return - 1 * simResults['medianIsp'] + 1 * simResults['medianPc'] + 1 * simResults['medianTc'] + 1 * simResults['mprop'] 

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
            'ThroatResizeCoeff' : self.rangeGenerator(self.confData['startingStates']['ThroatResizeCoeff']),
        }
        


if __name__ == "__main__":
    optimiser = engineOptimiser()
    optimiser.optimise()