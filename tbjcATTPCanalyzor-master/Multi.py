import functools
from Analyzor.DataFactory import DataFactory
from Analyzor import VertexAnalyzor###############
import time
from subprocess import call
from multiprocessing import Pool
import operator
import json
import numpy as np

def runProcess(exe):
    p = subprocess.Popen(exe,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    err = p.stderr.readlines()

    if len(err) != 0:
        for line in err:
            print(line)
        exit(1)

    out = p.stdout.readlines()
    return out

SQLpath ="/afs/crc.nd.edu/group/nsl/activetarget/users/jlai1/SQL/10C/"

def ProcessFile(fname):

    r = fname.split('/')[-1].split('.')[0].split('_')
    runID,fID = map(int,r)

    print(fname)

    dp = DataFactory(fname,SQLpath+'ProtoMap.db')
    dp.InitT3()
    dist = []

    for i in sorted(dp.t3['EventID'].unique()):
        try:
            image = dp.ConstructImage(i)
            image = VertexAnalyzor.FilterBackground(image)
            points,(xc,yc) = VertexAnalyzor.GetEventPositions(image,0) ### was 0 for debug mode
            xscale , yscale = 80/1000.0*2.459,0.1 # 2.549 cm/ns drift velocity
            xc,yc = xc*xscale, yc*yscale
            points = [(float(x*xscale),float(y*yscale)) for x,y in points]+[(xc,yc)]
            dist.append({'runID':runID,'fileID':fID,'eventID':i.item(),'data':points})
            #print(type(runID),type(fID),type(i.item()),type(points),type(points[0][0]))

        except:
            #print("eventID type is:",type(i),"\tproblema")
            pass

    return dist

if __name__ == "__main__":

    pool = Pool(processes=4) #was 4

    param = []
    runs = [85,]

    for run in runs:
        for i in range(1): # was range(1)
            path = SQLpath+'{:04d}_{:04d}.db'.format(run,i)
            param.append(path)

    start_time = time.time()

    res = pool.map(ProcessFile,param)

    end_time = time.time()

    print("total process cost "+str(end_time - start_time))
    data = functools.reduce(operator.add,res)

    with open('dataV2.dat','w') as f:
        json.dump(data,f,indent=4) # the form of data was not serializable, cast to string
