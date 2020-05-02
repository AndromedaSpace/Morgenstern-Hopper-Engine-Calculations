from rocketcea.cea_obj import add_new_fuel, add_new_propellant
from rocketcea.cea_obj_w_units import CEA_Obj
import time
from multiprocessing import Process, Queue


card_str = """
    fuel C50H102(S)  C 50.0   H 102.0 wt%=100.00
    h,cal= -343.7     t(k)=298.15   rho=.0333816
"""
add_new_fuel( 'paraffin', card_str )


Pe = 101325 #Pa

EPSmin = 1
EPSmax = 1.8 #3
dEPS = 0.05

Pmin = 2 * 10**5 # Pa
Pmax = 2.2 * 10**5 # Pa 30 
dP = 0.05 * 10**5 # Pa

OFmin = 2
OFmax = 12
dOF = 0.05

ox = "N2O"
fuel = "paraffin"

filename = "cea_results.txt"

nProcesses = 3



def worker(q,EPSmin, EPSmax , dEPS , Pmin , Pmax, dP, OFmin , OFmax, dOF, printProgress):
    C = CEA_Obj(oxName=ox, fuelName=fuel, pressure_units='Pa')
    cEPS = EPSmin
    if  (printProgress):
        last_time = time.time()
    while cEPS <= EPSmax:
        cP = Pmin
        while cP <= Pmax:
            cOF = OFmin
            while cOF <= OFmax:
                Ivac,Cstr,Tc = C.get_IvacCstrTc(Pc=cP, MR=cOF, eps=cEPS)

                Cf = C.get_PambCf(Pamb=Pe, Pc=cP, MR=cOF, eps=cEPS)

                q.put([cEPS,cP,cOF,Ivac,Cstr,Tc,Cf[1],Cf[2]])

                cOF += dOF
            cP += dP
        cEPS += dEPS
        if  (printProgress):
            cTime = time.time()
            print("Progress: %.2f%% ETA: %.2fs" % ((cEPS-EPSmin)/(EPSmax-EPSmin) *100 , (cTime-last_time) * (EPSmax-cEPS)/dEPS) )
            last_time = cTime
    return

if __name__ == '__main__':
    EPS_step = (EPSmax-EPSmin)/nProcesses
    jobs = []
    q = Queue()
    print("Generating Processes")
    for i in range(nProcesses):
        pEPSmin = EPSmin+EPS_step*i
        pEPSmax = EPSmin+EPS_step*(i+1)
        #sprint("EPSmin:" , pEPSmin, "EPSmax", pEPSmax)
        if i == nProcesses - 1:
            printProgress = True
        else:
            printProgress = False
        
        p = Process(target=worker, args=(q,pEPSmin,pEPSmax,dEPS,Pmin,Pmax,dP,OFmin,OFmax,dOF,printProgress))
        jobs.append(p)
        p.start()
    print("Processes Started.\nGenerating Data.")
    data = [[],[],[],[],[],[],[],[]]
    for job in jobs:
        while job.is_alive():
            while not q.empty():
                entry = q.get()
                for i in range(len(entry)):
                    data[i].append(entry[i])
    print("Data Collected.")

    while not q.empty():
        data.append(q.get())
    
    indices = list(range(len(data[0])))
    indices.sort(key=data[0].__getitem__)

    for i , sublist in enumerate(data):
        data[i] = [sublist[j] for j in indices]
    print("Writing data.")
    outFile = open(filename, 'w')
    outFile.write("#Format: eps Pc OF Ivac Cstr Tc Cf SeparetionState")
    for i in range(len(data[0])):
        outStr = ""
        for j in range(len(data)):
            outStr += str(data[j][i]) + ' '
        outStr = outStr[:-1] + '\n'
        outFile.write(outStr)
    outFile.close()
    print("Data written. Exiting.")