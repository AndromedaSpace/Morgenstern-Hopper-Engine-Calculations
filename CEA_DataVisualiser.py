from ceaDataReader import ceaDataReader
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.tri as mtri
import keyboard
import time
import matplotlib.animation as animation

interval = 100

plt.style.use('dark_background')
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')


dataReader = ceaDataReader()

dataReader.readData("cea_results.txt")

packedData , EpsVals = dataReader.packAndGetData()
cEPS = 0
dataToShow = 2


def showPlot(ax , x , y , z , state):
    ax.clear()
    ax.plot_trisurf(x,y,z,cmap=cm.coolwarm)
    
def main_animated(i):
    global cEPS , dataToShow, dataCut
    if keyboard.is_pressed('Right'):  
        if cEPS < len(packedData) - 1:
            cEPS += 1
            
    elif keyboard.is_pressed('Left'):
        if cEPS > 0:
            cEPS -= 1
            
    elif keyboard.is_pressed('Up'):
        if dataToShow < len(packedData[cEPS])-2:
            dataToShow += 1
            dataCut = len(packedData[cEPS][dataToShow]) - 1
            
    elif keyboard.is_pressed('Down'):
        if dataToShow > 2:
            dataToShow -= 1
            dataCut = len(packedData[cEPS][dataToShow]) - 1
                    
    
    showPlot(ax1 , packedData[cEPS][0] , packedData[cEPS][1] , packedData[cEPS][dataToShow] , packedData[cEPS][-1])


ani = animation.FuncAnimation(fig, main_animated, interval = interval)
plt.show()