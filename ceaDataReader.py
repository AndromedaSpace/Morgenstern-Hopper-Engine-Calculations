from copy import deepcopy

class ceaDataReader:
    data = []
    valNames = []

    def readDataFromFile(self,fileName):
        inFile = open(fileName, 'r')

        self.data = inFile.read()

        self.data = self.data.split("\n")
        self.valNames = self.data[0].split(' ')[2:]
        self.data.pop(0)
        self.data.pop(len(self.data)-1)

       

    def packAndGetData(self):
        template = [[],[],[],[],[],[],[],[]]
        lastEPS = ""
        packed = []
        EpsVals = []
        for line in self.data:
            line = line.split(' ')
            if line[-3] == "Separated":
                    line[-2] += line[-1]
                    line.pop(len(line)-1)
            if line[0] != lastEPS:
                lastEPS = line[0]
                EpsVals.append(line[0])
                packed.append(deepcopy(template))
            
            for i , entry in enumerate(line[1:]):
                if i == len(line[1:]) - 2:
                    packed[-1][i].append({
                        'state' : entry,
                        'data' : ""
                    })    
                elif i == len(line[1:]) - 1:
                    packed[-1][i-1][-1]['data'] = entry
                else:
                    packed[-1][i].append(float(entry))


        return packed , EpsVals


    def getDataByType(self):
        dataEps = []
        dataPc = []
        dataOF = []
        dataIvac = []
        dataCstr = []
        dataTc = []
        dataCf = []
        dataSepState = []
        for line in self.data:
            line = line.split(' ')
            dataEps.append(float(line[0]))
            dataPc.append(float(line[1]))
            dataOF.append(float(line[2]))
            dataIvac.append(float(line[3]))
            dataCstr.append(float(line[4]))
            dataTc.append(float(line[5]))
            dataCf.append(float(line[6]))
            dataSepState.append({
                'state' : line[7],
                'data' : line[8]
            })
        return dataEps, dataPc , dataOF , dataIvac , dataCstr , dataTc , dataCf , dataSepState
    
    def getValNames(self):
        return self.valNames