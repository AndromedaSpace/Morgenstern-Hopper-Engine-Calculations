import matplotlib.pyplot as plt
from rocketcea.cea_obj import add_new_fuel, add_new_propellant
from rocketcea.cea_obj_w_units import CEA_Obj


card_str = """
    fuel C50H102(S)  C 50.0   H 102.0 wt%=100.00
    h,cal= -343.7     t(k)=298.15   rho=.0333816
"""
add_new_fuel( 'paraffin', card_str )




Pe = 101325 #Pa
Pmin = 2 * 10**5 # Pa
Pmax = 3 * 10**5 # Pa
dP = 0.05 * 10**5 # Pa

OFmin = 2
OFmax = 12
dOF = 0.05

ox = "N2O"
fuel = "paraffin"

outFile = open("cea_resutls.txt", 'w')
outFile.write("#Format : Pc OF eps Ivac Cstr Tc Cf\n")
C = CEA_Obj(oxName=ox, fuelName=fuel, pressure_units='Pa')

cP = Pmin
dataPc = []
dataOF = []
dataEps = []
dataIvac = []
dataCstr = []
dataTc = []
dataCf = []
totalPSteps = (Pmax-Pmin)/dP
while cP <= Pmax:
    cOF = OFmin
    print("Progress", str((cP-Pmin)/(Pmax-Pmin) *100),'%')
    while cOF <= OFmax:
        outFile.write(str(cP) + ' ')
        outFile.write(str(cOF) + ' ')
        #dataPc.append(cP)
        #dataOF.append(cOF)


        eps = C.get_eps_at_PcOvPe(Pc=cP, MR=cOF, PcOvPe=cP/Pe)
        outFile.write(str(eps) + ' ')
        #dataEps.append(eps)

        Ivac,Cstr,Tc = C.get_IvacCstrTc(Pc=cP, MR=cOF, eps=eps)
        outFile.write(str(Ivac) + ' ')
        outFile.write(str(Cstr) + ' ')
        outFile.write(str(Tc) + ' ')
        #dataIvac.append(Ivac)
        #dataCstr.append(Cstr)
        #dataTc.append(Tc)

        Cf = C.get_PambCf(Pamb=Pe, Pc=cP, MR=cOF, eps=eps)
        outFile.write(str(Cf[1]) + '\n')
        #dataCf.append(Cf[1])

        cOF += dOF
    cP += dP

outFile.close()