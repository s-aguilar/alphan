from Analyzor.DataFactory import DataFactory
from Analyzor import VertexAnalyzor
import pylab as plt
import matplotlib.path as mplPath
import seaborn as sns
import numpy as np

plt.rcParams['figure.figsize'] = (20, 20)


dp = DataFactory('0085_0000.db','ProtoMap.db')
dp.InitT3()

dp.InitMesh()
dp.mesh_df.loc[6].plot()

image = dp.ConstructImage(7)
image = VertexAnalyzor.FilterBackground(image)
points,(xc,yc) = VertexAnalyzor.GetEventPositions(image,1)
plt.scatter(xc,yc,s=200,c='r')

plt.show()
