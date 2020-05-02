from rocketcea.cea_obj import add_new_fuel, add_new_propellant
from rocketcea.cea_obj_w_units import CEA_Obj
import time


card_str = """
    fuel C50H102(S)  C 50.0   H 102.0 wt%=100.00
    h,cal= -343.7     t(k)=298.15   rho=.0333816
"""
add_new_fuel( 'paraffin', card_str )




Pe = 101325 #Pa

EPSmin = 1
EPSmax = 3
dEPS = 0.05

Pmin = 2 * 10**5 # Pa
Pmax = 30 * 10**5 # Pa
dP = 0.05 * 10**5 # Pa

OFmin = 2
OFmax = 12
dOF = 0.05

ox = "N2O"
fuel = "paraffin"

outFile = open("cea_results.txt", 'w')
C = CEA_Obj(oxName=ox, fuelName=fuel, pressure_units='Pa')

cP = Pmin

totalPSteps = (EPSmax-EPSmin)/dP
cEPS = EPSmin
last_time = time.time()
while cEPS < EPSmax:
    cP = Pmin
    while cP <= Pmax:
        print(cP)
        cOF = OFmin
        while cOF <= OFmax:
            outFile.write(str(cEPS) + ' ')
            outFile.write(str(cP) + ' ')
            outFile.write(str(cOF) + ' ')

            Ivac,Cstr,Tc = C.get_IvacCstrTc(Pc=cP, MR=cOF, eps=cEPS)
            outFile.write(str(Ivac) + ' ')
            outFile.write(str(Cstr) + ' ')
            outFile.write(str(Tc) + ' ')

            Cf = C.get_PambCf(Pamb=Pe, Pc=cP, MR=cOF, eps=cEPS)
            outFile.write(str(Cf[1]) + '\n')

            cOF += dOF
        cP += dP
    cEPS += dEPS
    cTime = time.time()
    print("Progress", str((cEPS-EPSmin)/(EPSmax-EPSmin) *100),'%. ETA: ' , str((cTime-last_time) * (EPSmax-cEPS)/dEPS), 's')
    last_time = cTime

outFile.close()