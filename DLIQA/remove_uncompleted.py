# -*- coding: utf-8 -*-
import sys
import os
import numpy as np

import glob



def main():

    file_nr=0
    removeables=[]
    base_name=[]
    for file in os.listdir('./data/in_work'):
        if file.endswith(".work"):

            removeables.append(file[:-5])
            name_split=file.split('-')

            base_name.append(name_split[0])

    print(len(removeables))
		
    for count,remove_name in enumerate(removeables):
        print(remove_name)
        for dir,_,_ in os.walk('./data/'+base_name[count]):
            for file_r in glob.glob(dir+'/'+remove_name+'*'):
                file_nr=file_nr+1
                os.remove(file_r)
                print file_r
        for file_w in glob.glob('./data/in_work/'+remove_name+'.work'):
			#print file
            file_nr=file_nr+1
            os.remove(file_w)
		#print('-----------------------------------------------------')
	#/data/protein_name/

    print file_nr

if __name__ == '__main__':
    main()

