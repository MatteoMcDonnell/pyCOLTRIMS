import numpy as np 
import time
from convertGUI import GUI
import ast
from itertools import combinations

####### Setup / Parameters #######


fileIn,fileOut,channelBounds = GUI()

#convert the input into an actual list of lists



starttime=time.time()

#use if you want to break the loop after a certain number of events.
#if not, set to zero
breakNum = 0


def importData(fileLoc):

    f = open(fileLoc)
    data =  np.fromfile(f,dtype = 'int32', count = -1)
    f.close()
    
    return data

D = importData(fileIn)


####### Data Shaping #######


def makeRow(itt):
    #creates row (or set array of rows) for each event
    
    global D
    count = D[0]
    delay = D[1]

    if count == 0:

        row_ = np.array([delay,np.nan,np.nan,np.nan,itt])
        D = D[2:] 

    else: 
        
        row_ = np.empty((count,5))

        for i in range(count):
           
            #replaces first three columns in row with pos/time data for that hit
            row_[0] = delay 
            row_[i,1:4] = D[2+(3*i):2+(3*(i+1))] 
            row_[i,4] = itt          
       
        #deletes the parsed data from the data series
        D = D[(2+(3*count)):]

    return row_


def dataToArray(data):

    #creates the array of data using the above helper function 
    arrayList = []

    n = 1 

    while len(D) > 1:

        row = makeRow(n)

        arrayList.append(row)

        n += 1

        if  n%500000==0:
            print(n,'events processed')

        #if you want to break early:
        if (breakNum != 0)and(n == breakNum):
            break

    return np.vstack(arrayList)


a = dataToArray(D)



####### POST PROCESSING #######

#if the units of x,y,or tof are not in mm or ns, modify this. 
a[:,:3] *= (1/1000) 


time1 = time.time()

def toCSV(array): #for outputting everything straight to csv files with no cutting at all

    print('Writing to CSV')

    np.savetxt(fileOut[:-3]+'csv',a,delimiter=",") 


def slicing(array):
    #for outputting a set of files based off of the set of parameters 
    #in the channelBounds list

    import pandas as pd

    channelBounds = [ast.literal_eval(lis) for lis in channelBounds]
    
    
    pd.set_option('mode.chained_assignment',None)

    print('Converting to DataFrame')
    frame = pd.DataFrame(array,columns = ['delay','x','y','tof','id'])
    
    #only list coincidences 
    value_counts = frame['id'].value_counts() 
    coinc = frame[frame['id'].isin(value_counts[value_counts > 1].index)]
    coinc.reset_index(inplace=True,drop=True)

    print('Organizing pairwise coincidences, this will take a while...')
    print('(like a really, really long time)')        

    coinc['index0'] = coinc.index
    pairwise = pd.DataFrame([[k, c0, c1] for k, index in coinc.groupby('id').index0
                                         for c0, c1 in combinations(index, 2)
                                       ], columns=['id', 'index1', 'index2'])


    def sumdifGater(data,pairwise,gates):
        
        temp = pairwise
                
        temp['sum'] = data.iloc[temp['index1'].values].tof.values + data.iloc[temp['index2'].values].tof.values
        temp['dif'] = data.iloc[temp['index2'].values].tof.values - data.iloc[temp['index1'].values].tof.values
                            
        temp['checkSum'] = (gates[2] < temp['sum']) & (temp['sum'] < gates[3])
        temp['checkDif'] = (gates[0] < temp['dif']) & (temp['dif'] < gates[1])
                                        
        output = pairwise.loc[temp['checkSum'] & temp['checkDif']][['id','index1','index2']].reset_index(drop=True)
                                                
        return output
   
    
    #I/O is the bit that will take a while, everything else is vectorized
    for bounds in channelBounds:

        current = sumdifGater(coinc,pairwise,bounds)
           
        compList = pd.DataFrame(np.sort(np.append(current['index1'].values,current['index2'].values)),columns=['index0']).drop_duplicates().reset_index(drop=True)  
        
        semi = coinc.loc[compList['index0'].values] 
        id_s = semi.id.drop_duplicates()
       
        
        res = coinc[coinc.id.isin(id_s)]
        res.drop(columns = ['index0'],inplace=True)

        boundString =  '_'.join(str(e) for e in bounds)
        
        outputName = fileOut[:-4] + boundString + '.csv'

        res.to_csv(outputName, index = False) 


if channelBounds == None:
    toCSV(a)

else:
    slicing(a)


print(a.shape[0],"detections, at",a.shape[0]/(time1-starttime),'Hz')
print(time.time()-starttime,'seconds total,',(time.time()-time1),'for post-processing')
