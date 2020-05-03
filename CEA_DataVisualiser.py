from ceaDataReader import ceaDataReader
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.tri as mtri

plt.style.use('dark_background')
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')


dataReader = ceaDataReader()

dataReader.readData("cea_results.txt")

packeData , EpsVals = dataReader.packAndGetData()
print(len(packeData))
print(EpsVals)
