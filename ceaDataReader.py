class ceaDataReader(object):
    dataPc = []
    dataOF = []
    dataEps = []
    dataIvac = []
    dataCstr = []
    dataTc = []
    dataCf = []

    def readData(self,fileName):
        inFile = open(fileName, 'r')

        data = inFile.read()

        data = data.split("\n")

        data.pop(0)
        data.pop(len(data)-1)

        for line in data:
            line.split(' ')
            print(line)
            self.dataPc.append(line[0])
            self.dataOF.append(line[1])
            self.dataEps.append(line[2])
            self.dataIvac.append(line[3])
            self.dataCstr.append(line[4])
            self.dataTc.append(line[5])
            self.dataCf.append(line[6])

    def getData(self):
        return self.dataPc , self.dataOF , self.dataEps , self.dataIvac , self.dataCstr , self.dataTc , self.dataCf 