import sys
import re
#import matplotlib.pyplot as plt
#from matplotlib import cm
import math
import numpy as np
import os
import time
import glob




#------------------------------------------------------------------------------------
def check_in_use(name, file_folder):

    if os.path.exists('./data/in_work/'+name+'.work') or os.path.exists(file_folder+'/'+name+'_feature_rips_charge_1D.npy'):
        existing='true'
    else:
        if not os.path.exists("./data/in_work"):
            os.makedirs("./data/in_work")
        file_t =open('./data/in_work/'+name+'.work', 'w')
        file_t.close()
        existing='false'
    return(existing)


#------------------------------------------------------------------------------------
def main():
    
    import random
    #Random time for it to wait so two processes don't start at the same time
    sleep_time=random.uniform(1, 50)

    count_nr=0    
    #python DLIQA.py -input_file
    #filename_pdb = '../../../CnM-dataset/MOAL_Benchmark/D1A2K/D1A2K-a0a-merged.pdb'  #sys.argv[1:]
    for x in range(1,6):
        cross=open('./cross_val_sets/cross_val_set_'+str(x))
        lines = cross.read().splitlines()
        print x
        print lines
    for f in lines:
        pdb_files=[]
        for file_p in glob.glob('../../../CnM-dataset/'+f+'/*.pdb'):
            pdb_files.append(file_p)
        for filename_pdb in pdb_files:

            
            name_split=filename_pdb.split('/')
            file_folder='./data/'+name_split[-2]
            
            name_split2=name_split[-1].split('-')
            
            #This is just for naming just my dataset fine, else remove.      
            if len(name_split2)==3:
                name=name_split[-1][:-11]
            else:
                name=name_split[-1]
            #print(name)     
            
            existing=check_in_use(name,file_folder)
            if existing=='true':
                continue

            else:
                count_nr=count_nr+1
       
    print(count_nr)
if __name__ == '__main__':
    main()





