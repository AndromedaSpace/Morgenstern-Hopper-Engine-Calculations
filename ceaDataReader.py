class ceaDataReader:
    dataPc = []
    dataOF = []
    dataEps = []
    dataIvac = []
    dataCstr = []
    dataTc = []
    dataCf = []
    dataSepState = []

    def readData(self,fileName):
        inFile = open(fileName, 'r')

        data = inFile.read()

        data = data.split("\n")

        data.pop(0)
        data.pop(len(data)-1)

        for line in data:
            line = line.split(' ')
            self.dataEps.append(float(line[0]))
            self.dataPc.append(float(line[1]))
            self.dataOF.append(float(line[2]))
            self.dataIvac.append(float(line[3]))
            self.dataCstr.append(float(line[4]))
            self.dataTc.append(float(line[5]))
            self.dataCf.append(float(line[6]))
            self.dataSepState.append({
                'state' : line[7],
                'data' : line[8]
            })

    def getData(self):
        return self.dataEps, self.dataPc , self.dataOF , self.dataIvac , self.dataCstr , self.dataTc , self.dataCf , self.dataSepState