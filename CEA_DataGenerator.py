from rocketcea.cea_obj import add_new_fuel, add_new_propellant
from rocketcea.cea_obj_w_units import CEA_Obj
import time
from multiprocessing import Process, Queue

class CEADataGenerator:
    card_str = """
        fuel C50H102(S)  C 50.0   H 102.0 wt%=100.00
        h,cal= -343.7     t(k)=298.15   rho=.0333816
    """
    add_new_fuel('paraffin', card_str )


    ox = ""
    fuel = ""

    data = []

    def setOX(self, oxName = "N2O"):
        self.ox = oxName

    def setFuel(self, fuelName = "paraffin"):
        self.fuel = fuelName

    def singleWorker(self,Pe,EPS, P,OF):
        C = CEA_Obj(oxName=self.ox, fuelName=self.fuel, pressure_units='Pa')
        Ivac,Cstr,Tc = C.get_IvacCstrTc(Pc=P, MR=OF, eps=EPS)
        Isp = C.estimate_Ambient_Isp(Pc=P,MR=OF,eps=EPS, Pamb=Pe)[0]
        Cstr /=  3.2808
        Cf = C.get_PambCf(Pamb=Pe, Pc=P, MR=OF, eps=EPS)
        return [
            Ivac,
            Isp,
            Cstr,
            Tc,
            {
                'state' : Cf[1],
                'data' : Cf[2]
            }
        ]


    def rangeWorker(self,q,Pe,EPSmin, EPSmax , dEPS , Pmin , Pmax, dP, OFmin , OFmax, dOF, printProgress):
        C = CEA_Obj(oxName=self.ox, fuelName=self.fuel, pressure_units='Pa')
        cEPS = EPSmin
        if  (printProgress):
            last_time = time.time()
        while cEPS <= EPSmax:
            cP = Pmin
            while cP <= Pmax:
                cOF = OFmin
                while cOF <= OFmax:
                    Ivac,Cstr,Tc = C.get_IvacCstrTc(Pc=cP, MR=cOF, eps=cEPS)
                    Isp = C.estimate_Ambient_Isp(Pc=cP,MR=cOF,eps=cEPS, Pamb=Pe)[0]
                    Cstr /=  3.2808
                    Cf = C.get_PambCf(Pamb=Pe, Pc=cP, MR=cOF, eps=cEPS)

                    q.put([cEPS,cP,cOF,Ivac,Isp,Cstr,Tc,Cf[1],Cf[2]])

                    cOF += dOF
                cP += dP
            cEPS += dEPS
            if  (printProgress):
                cTime = time.time()
                print("Progress: %.2f%% ETA: %.2fs" % ((cEPS-EPSmin)/(EPSmax-EPSmin) *100 , (cTime-last_time) * (EPSmax-cEPS)/dEPS) )
                last_time = cTime
        return

    def runRange(self,nProcesses = 3,Pe=101325,EPSmin=1, EPSmax=1.8 , dEPS=0.05 , Pmin=2 * 10**5  , Pmax= 30 * 10**5, dP = 0.05 * 10**5, OFmin = 2 , OFmax = 10, dOF = 0.05):
        EPS_step = (EPSmax-EPSmin)/nProcesses
        jobs = []
        q = Queue()
        print("Generating Processes")
        for i in range(nProcesses):
            pEPSmin = EPSmin+EPS_step*i
            pEPSmax = EPSmin+EPS_step*(i+1)
            if i == nProcesses - 1:
                printProgress = True
            else:
                printProgress = False
            
            p = Process(target=self.rangeWorker, args=(q,Pe,pEPSmin,pEPSmax,dEPS,Pmin,Pmax,dP,OFmin,OFmax,dOF,printProgress))
            jobs.append(p)
            p.start()
        print("Processes Started.\nGenerating Data.")
        self.data = [[],[],[],[],[],[],[],[],[]]
        for job in jobs:
            while job.is_alive():
                while not q.empty():
                    entry = q.get()
                    for i in range(len(entry)):
                        self.data[i].append(entry[i])
        while not q.empty():
            self.data.append(q.get())

        indices = list(range(len(self.data[0])))
        indices.sort(key=self.data[0].__getitem__)

        for i , sublist in enumerate(self.data):
            self.data[i] = [sublist[j] for j in indices]
        
        print("Data Collected.")

    def saveDataToFile(self,filename = "cea_results.txt"):
        print("Writing data.")
        outFile = open(filename, 'w')
        outFile.write("#Format: eps Pc(Pa) OF Ivac(s) Isp(s) Cstr(m/s) Tc(K) Cf SeparetionState\n")
        for i in range(len(self.data[0])):
            outStr = ""
            for j in range(len(self.data)):
                outStr += str(self.data[j][i]) + ' '
            outStr = outStr[:-1] + '\n'
            outFile.write(outStr)
        outFile.close()
        print("Data written.")
    
    def clearData(self):
        self.data = []


if __name__ == '__main__':
    generator = CEADataGenerator()
    generator.setOX("N2O")
    generator.setFuel()
    generator.runRange(Pmax=2.2 * 10**5)
    generator.saveDataToFile()