import json
from engineSimulator import engineSimulator

class engineOptimiser:
    setupData = None
    engine = None

    def __init__ (self,configFile="optimiser.conf"):
        self.setupData = self.readData(filename=configFile)
        self.engine = self.initEngine(self.setupData['engine'])


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



if __name__ == "__main__":
    optimiser = engineOptimiser()