from ceaDataReader import ceaDataReader
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.tri as mtri
import keyboard
import time

plt.style.use('dark_background')
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')


dataReader = ceaDataReader()

dataReader.readData("cea_results.txt")

packedData , EpsVals = dataReader.packAndGetData()
cEPS = 0
dataToShow = 2
dataCut = len(packedData[cEPS][dataToShow]) - 1

while True:
    if keyboard.is_pressed('esc'):
        print('Exiting')
        break
    elif keyboard.is_pressed('Up'):  
        if cEPS < len(packedData) - 1:
            cEPS += 1
            time.sleep(0.2)
    elif keyboard.is_pressed('Down'):
        if cEPS > 0:
            cEPS -= 1
            time.sleep(0.2)
    elif keyboard.is_pressed('PgUp'):
        if dataToShow < len(packedData[cEPS])-2:
            dataToShow += 1
            dataCut = len(packedData[cEPS][dataToShow]) - 1
            time.sleep(0.2)
    elif keyboard.is_pressed('PgDown'):
        if dataToShow > 2:
            dataToShow -= 1
            dataCut = len(packedData[cEPS][dataToShow]) - 1
            time.sleep(0.2)
    elif keyboard.is_pressed('Right'):
        length = len(packedData[cEPS][dataToShow])
        dataCut += length * 0.05
        if dataCut > length - 1:
            dataCut = length
    elif keyboard.is_pressed('Down'):
        length = len(packedData[cEPS][dataToShow])
        dataCut -= length * 0.05
        if dataCut < 0:
            dataCut = 0

def showPlot(ax , x , y , z , state):
    return 