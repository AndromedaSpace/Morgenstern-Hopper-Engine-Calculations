from ceaDataReader import ceaDataReader
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.tri as mtri

plt.style.use('dark_background')
fig = plt.figure()
ax1 = fig.add_subplot(321, projection='3d')
ax2 = fig.add_subplot(322, projection='3d')
ax3 = fig.add_subplot(323, projection='3d')
ax4 = fig.add_subplot(324, projection='3d')
ax5 = fig.add_subplot(325, projection='3d')

dataReader = ceaDataReader()

dataReader.readData("cea_results.txt")

dataPc , dataOF , dataEps , dataIvac , dataCstr , dataTc , dataCf = dataReader.getData()

ax1.set_title("Vaccume Isp vs Pc vs OF")
ax1.plot_trisurf(dataPc, dataOF , dataIvac, cmap=cm.coolwarm)

ax2.set_title("Cstr vs Pc vs OF")
ax2.plot_trisurf(dataPc, dataOF , dataCstr, cmap=cm.coolwarm)

ax3.set_title("Cf vs Pc vs OF")
ax3.plot_trisurf(dataPc, dataOF , dataCf, cmap=cm.coolwarm)

ax4.set_title("Tc vs Pc vs OF")
ax4.plot_trisurf(dataPc, dataOF , dataTc, cmap=cm.coolwarm)

ax5.set_title("Expansion Ration vs Pc vs OF")
ax5.plot_trisurf(dataPc, dataOF , dataEps, cmap=cm.coolwarm)

plt.show()