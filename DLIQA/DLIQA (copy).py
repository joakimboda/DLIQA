import Bio
from Bio.PDB import *
import sys
import re
#import matplotlib.pyplot as plt
#from matplotlib import cm
import math
import numpy as np
import os

import pandas
import multiprocessing

import shutil
import pickle


#This is for removing Bio.python warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

#------------------------------------------------------------------------------------

def Feature_rips_1D(name, file_folder):
    print('Generating Rips Interaction Feature...')
    ProEleList = ['C','N','O','S']
    LigEleList = ['C','N','O','S','P','F','Cl','Br','I']

    Bins = []
    for i in range(201):
        Bins.append(float(i)*0.25)

    def BinID(x, B):
        for iii in range(0,len(B)-1):
            if x >= B[iii] and x <= B[iii+1]:
                y = iii
        return y
    cut = 12.0

    Feature = np.zeros([len(ProEleList), len(LigEleList), len(Bins)-1], float)

    for j in range(len(ProEleList)):
        ep = ProEleList[j]
        for k in range(len(LigEleList)):
            el = LigEleList[k]

            

            if not os.path.exists(file_folder+'/bar_files/'+name+'_'+ep+'-'+el+'.bar'):
                continue
            InFile = open(file_folder+'/bar_files/'+name+'_'+ep+'-'+el+'.bar')
            lines = InFile.read().splitlines()
            for line in lines:
                dim,b,d = line.split()
                if float(d) >= cut:
                    continue
                did = BinID(float(d),Bins)
                Feature[j,k,did] += 1.0

    FeatureFlat = np.zeros([len(Bins)-1, len(ProEleList)*len(LigEleList)], float)
    for i in range(len(ProEleList)):
        for j in range(len(LigEleList)):
            fid = i*len(LigEleList)+j
            for k in range(len(Bins)-1):
                FeatureFlat[k,fid] = Feature[i,j,k]

    OutFile = open(file_folder+'/'+name+'_feature_rips_1D.npy', 'w')
    np.save(OutFile, FeatureFlat)
    OutFile.close()

#------------------------------------------------------------------------------------

def Feature_rips_charge_1D(name, file_folder):
    print('Generating Rips Charge Feature...')
    ProEleList = ['C','N','O','S','H']
    LigEleList = ['C','N','O','S','P','F','Cl','Br','I','H']

    Bins = []
    for i in range(101):
        Bins.append(float(i)*0.01)

    def BinID(x, B):
        for iii in range(0,len(B)-1):
            if x >= B[iii] and x <= B[iii+1]:
                y = iii
        return y
    cut = 1.0

    Feature = np.zeros([len(ProEleList), len(LigEleList), len(Bins)-1], float)

    for j in range(len(ProEleList)):
        ep = ProEleList[j]
        for k in range(len(LigEleList)):
            el = LigEleList[k]
            if not os.path.exists(file_folder+'/bar_files/'+name+'_'+ep+'-'+el+'_charge.bar'):
                continue
            InFile = open(file_folder+'/bar_files/'+name+'_'+ep+'-'+el+'_charge.bar')
            lines = InFile.read().splitlines()
            for line in lines:
                dim,b,d = line.split()
                if float(d) >= cut:
                    continue
                did = BinID(float(d),Bins)
                Feature[j,k,did] += 1.0

    FeatureFlat = np.zeros([len(Bins)-1, len(ProEleList)*len(LigEleList)], float)

    for i in range(len(ProEleList)):
        for j in range(len(LigEleList)):
            fid = i*len(LigEleList)+j
            for k in range(len(Bins)-1):
                FeatureFlat[k,fid] = Feature[i,j,k]

    OutFile = open(file_folder+'/'+name+'_feature_rips_charge_1D.npy', 'w')
    np.save(OutFile, FeatureFlat)
    OutFile.close()


#------------------------------------------------------------------------------------

def Feature_alpha_2D(name, file_folder):
    print('Generating Alpha Feature...')
    
    thr = 12.0
    rs = 0.1
    lth = int(thr/rs)
    small = 0.01

    OrderedName = np.load('PerformanceOrderAlphaHand.npy')
    X = np.zeros([8, lth, 128], float) #X = np.zeros([16, lth, 128], float)

    InFile = open(file_folder+'/'+name+'_alpha.pkl')
    BarCollection = pickle.load(InFile)
    for j in range(len(OrderedName)):
        plname = OrderedName[j]; pname, lname = plname.split('_');
        #f_p_i_j = np.zeros([8,lth], float)
        f_pl_i_j = np.zeros([8,lth], float)
        '''if 'pro_'+pname in BarCollection.keys():
            Bars = BarCollection['pro_'+pname]
            for Bar in Bars:
                if Bar[1] >= thr: continue
                if Bar[2] >= thr: Bar[2] = thr-0.00001
                bid = min(int(np.floor(Bar[1]/rs)), lth-1)
                did = min(int(np.floor(Bar[2]/rs)), lth-1)
                if Bar[0] == 0:
                    f_p_i_j[0,did] += 1.0
                    f_p_i_j[1,bid:did+1] += 1.0
                elif Bar[0] == 1:
                    f_p_i_j[2,bid] += 1.0
                    f_p_i_j[3,did] += 1.0
                    f_p_i_j[4,bid:did+1] += 1.0
                elif Bar[0] == 2:
                    f_p_i_j[5,bid] += 1.0
                    f_p_i_j[6,did] += 1.0
                    f_p_i_j[7,bid:did+1] += 1.0'''
        if 'com_'+plname in BarCollection.keys():
            Bars = BarCollection['com_'+plname]
            for Bar in Bars:
                if Bar[1] >= thr: continue
                if Bar[2] >= thr: Bar[2] = thr-0.00001
                bid = min(int(np.floor(Bar[1]/rs)), lth-1)
                did = min(int(np.floor(Bar[2]/rs)), lth-1)
                if Bar[0] == 0:
                    f_pl_i_j[0,did] += 1.0
                    f_pl_i_j[1,bid:did+1] += 1.0
                elif Bar[0] == 1:
                    f_pl_i_j[2,bid] += 1.0
                    f_pl_i_j[3,did] += 1.0
                    f_pl_i_j[4,bid:did+1] += 1.0
                elif Bar[0] == 2:
                    f_pl_i_j[5,bid] += 1.0
                    f_pl_i_j[6,did] += 1.0
                    f_pl_i_j[7,bid:did+1] += 1.0
        #f_df_i_j = f_pl_i_j[:,:] - f_p_i_j[:,:]
        X[0:8,:,j] = f_pl_i_j[:,:]; #X[8:16,:,j] = f_df_i_j[:,:]

    OutFile = open(file_folder+'/'+name+'_feature_alpha_2D.npy', 'w')
    np.save(OutFile, X)
    OutFile.close()


#------------------------------------------------------------------------------------

def alpha(filename, file_folder, alphaprot):

    print("Computing alpha complex based persistent homology for protein-ligand complex...")

    dt = np.dtype([('dim', int), ('birth', float), ('death', float)])
    ProEleCollection = [['C'],['N'],['O'],['S'],['C','N'],['C','O'],['N','O'],['C','N','O'],['C','N','O','S']]
    LigEleCollection = [['C'],['N'],['O'],['S'],['C','N'],['C','O'],['C','S'],['N','O'],['N','S'],['O','S'],['C','N','O'],['C','N','S'],['C','O','S'],['N','O','S'],['C','N','O','S'],['C','N','O','S','P','F','Cl','Br','I']]

    small = 0.01

    
    PRO = alphaprot['PRO']; LIG = alphaprot['LIG'];
    BarCollection = {}
    # alpha computation for complex
    
    for ep in ProEleCollection:
        for el in LigEleCollection:
            propts = []
            for a in range(len(PRO)):
                if PRO[a][0] in ep:
                    propts.append([ PRO[a][1], PRO[a][2], PRO[a][3]])
            ligpts = []
            for a in range(len(LIG)):
                if LIG[a][0] in el:
                    ligpts.append([ LIG[a][1], LIG[a][2], LIG[a][3] ])
            if len(propts) + len(ligpts) > 3:
                pname = ''
                for eep in ep:
                    pname = pname + eep
                lname = ''
                for eel in el:
                    lname = lname + eel
                name = 'com_'+pname+'_'+lname
                pt = propts + ligpts
                
                
                Bars = []
                tmpoutfile = open(file_folder+'/temp/'+filename+'pt.csv', 'w')
                tmpoutfile.write('x1,x2,x3\n')
                for pp in pt:
                    tmpoutfile.write(str(pp[0])+','+str(pp[1])+','+str(pp[2])+'\n')
                tmpoutfile.close()
                os.system('Rscript alpha.R '+file_folder+'/temp/'+filename+'pt.csv '+file_folder+'/temp/'+filename+'tmp.out')
                
                tmpinfile = open(file_folder+'/temp/'+filename+'tmp.out')
                lines = tmpinfile.read().splitlines()
                for line in lines[1:]:
                    a,b,c,d = line.split()
                    if d!='Inf':
                        if float(d)-float(c) >= small:
                            Bars.append([int(b), float(c), float(d)])
                BarCollection[name] = Bars
                os.system('rm '+file_folder+'/temp/'+filename+'pt.csv '+file_folder+'/temp/'+filename+'tmp.out')

    print('Saving Alpha file...')
    OutFile = open(file_folder+'/'+filename+'_alpha.pkl', 'w')
    #np.save(cwd+'/data/'+filename+'_alpha.npy', BarCollection)
    pickle.dump(BarCollection, OutFile, protocol=pickle.HIGHEST_PROTOCOL)
    OutFile.close()

#------------------------------------------------------------------------------------

def calc_residue_dist(residue1, residue2, outfile, outfile_charge, atom1_used, atom2_used, alphaprot) :

    pro=['C','N','O','S','H']
    lig=['C','N','O','S','P','F','Cl','Br','I','H']
    for atom1 in residue1:
        for atom2 in residue2:
            
            distance=atom1-atom2    

            if distance>20:#If the residues is to far from eachother there is no need to calculate the distance for the rest of the atoms
                return(alphaprot)
            if distance >12:#Go to the next atom
                continue
                
            a1_serial=str(atom1.serial_number)
            a2_serial=str(atom1.serial_number)
            
            if 'H' not in atom1.get_id() and 'H' not in atom2.get_id(): #Check if one of the connected atoms is an H then they are only used in Rips Charge
                  
                x1=str(atom1.get_coord()[0])
                y1=str(atom1.get_coord()[1])
                z1=str(atom1.get_coord()[2])

                x2=str(atom2.get_coord()[0])
                y2=str(atom2.get_coord()[1])
                z2=str(atom2.get_coord()[2])
                

                if a1_serial not in atom1_used:
                    atom1_used.append(a1_serial)
                    try:
                        list_index_p=pro.index(str(atom1.get_name()[0]))
                        ele=str(list_index_p+1)
                        outfile.write('0'+' '+ele+' '+x1+' '+y1+' '+z1+'\n')
                        outfile_charge.write('0'+' '+ele+' '+x1+' '+y1+' '+z1+' '+str(atom1.get_occupancy())+'\n')
                        alphaprot.setdefault('PRO', []).append([atom1.get_name()[0],x1,y1,z1])
                    except:
                    	'odd'
                if a2_serial not in atom2_used:
                    atom2_used.append(a2_serial)
                    try:
                        list_index_l=lig.index(atom2.get_name()[0])
                        ele=str(list_index_l+1)
                        outfile.write('1'+' '+ele+' '+x2+' '+y2+' '+z2+'\n')
                        outfile_charge.write('1'+' '+ele+' '+x2+' '+y2+' '+z2+' '+str(atom2.get_occupancy())+'\n')
                        alphaprot.setdefault('LIG', []).append([atom2.get_name()[0],x2,y2,z2])
                    except:
                            'odd'
            else:#if one of the connecting atoms is an H it only used in Rips charge
                if a1_serial not in atom1_used:
                    atom1_used.append(a1_serial)
                    try:
                        list_index_p=pro.index(atom1.get_name()[0])
                        ele=str(list_index_p+1)
                        outfile_charge.write('0'+' '+ele+' '+x1+' '+y1+' '+z1+' '+str(atom1.get_occupancy())+'\n')
                    except:
                         'odd'
                if a2_serial not in atom2_used:
                    atom2_used.append(a2_serial)
                    try:
                        list_index_l=lig.index(atom2.get_name()[0])
                        ele=str(list_index_l+1)
                        outfile_charge.write('1'+' '+ele+' '+x2+' '+y2+' '+z2+' '+str(atom2.get_occupancy())+'\n')
                    except:
                            'odd'
                
    return(alphaprot)




def calc_dist_matrix(chain_one, chain_two,outfile, outfile_charge, alphaprot) :
    atom1_used=[]
    atom2_used=[]
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            alphaprot=calc_residue_dist(residue_one, residue_two, outfile, outfile_charge, atom1_used, atom2_used, alphaprot)
                                    
    return(alphaprot)


#------------------------------------------------------------------------------------

def main():
    
    
    
    #python DLIQA.py -input_file
    #filename_pdb = '../../../CnM-dataset/MOAL_Benchmark/D1A2K/D1A2K-a0a-merged.pdb'  #sys.argv[1:]
    work_length=len(sys.argv[1:])
    for filename_pdb in sys.argv[1:]:
    
    #filename_pdb = sys.argv[1]  

#if sys_args[0]:
    #    pdb_dir=sys_args[0]
    #    sys_args = sys.argv[1:]
    #if sys_args[0]:
    #    valuefile=sys_args[0] #################################
    #    sys_args = sys.argv[1:]
    #if sys_args[0]:
    #    valuefile=sys_args[0]
    #    sys_args = sys.argv[1:]
    
   
    
    #cwd = os.getcwd()
    
    name_split=filename_pdb.split('/')
    file_folder='./data/'+name_split[-2]
    name=name_split[-1][:-11]
    

    if not os.path.exists(file_folder):
            os.makedirs(file_folder)
            os.makedirs(file_folder + "/pts_dir")
            os.makedirs(file_folder + "/pqr_files")
            os.makedirs(file_folder + "/temp")
    if not os.path.exists(file_folder + "/pts_dir"):
           os.makedirs(file_folder + "/pts_dir")
    if not os.path.exists(file_folder + "/pqr_files"):
           os.makedirs(file_folder + "/pqr_files")
    if not os.path.exists(file_folder + "/temp"):
           os.makedirs(file_folder + "/temp")
    if not os.path.exists(file_folder + "/bar_files"):
           os.makedirs(file_folder + "/bar_files")

    #working_file=open(cwd + '/data/pts_dir/working_file.txt', 'a')
    #working_file.write(name+'\n') #Writes what pdb have been worked on in working_file.txt
    #working_file.close()
           
    print(name)
    
    #filename_pdb=os.path.join(pdb_list[counter][0],pdb_list[counter][1])
    
    PDBobj = PDBParser()
    structure = PDBobj.get_structure(filename_pdb, filename_pdb)
    model = structure[0]
       
    try: 
        PDBobj = PDBParser()
        structure = PDBobj.get_structure(filename_pdb, filename_pdb)
        model = structure[0]
    
    except IOError: 
        print('IO Error', filename_pdb)      
        while 'true':
            input1=raw_input("Error parsing PDB file! Continue(y/n)...")
            if input1.lower()=='y':
                continue
            elif input1.lower()=='n':
                sys.exit()
    
    #Seperate the 2 chains for calculation in PQR
    io = PDBIO()
    structure_pqr=[]
    for chain in model.get_chains():
        io.set_structure(chain)
        io.save(file_folder+'/temp/'+name + "_" + chain.get_id() + ".pdb")
        #If PQR=1 it will try to run pdb2pqr
        PQR=1
        if PQR==1:
            infile=file_folder +'/temp/'+name + "_" + chain.get_id() + ".pdb"
            outfile=file_folder + '/pqr_files/'+name + "_" + chain.get_id() +'.pqr'
            print outfile
            line = '/software/apps/python/2.7.13/anaconda-5.0.0.1/bin/python /proj/wallner/apps/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py --ff=amber '+ infile +' '+ outfile
            os.system(line)
            structure_pqr.append(PDBobj.get_structure(outfile, outfile))
            os.remove(file_folder +'/temp/'+name + "_" + chain.get_id() + ".pdb")
    

    outfile=open(file_folder+ '/pts_dir/'+name + '.pts', 'w')
    outfile_charge=open(file_folder + '/pts_dir/'+name + '_charge.pts', 'w')
    alphaprot={}
    model1=structure_pqr[0][0]
    model2=structure_pqr[1][0]

    for chain1 in model1:
        for chain2 in model2:
            alphaprot=calc_dist_matrix(chain1, chain2, outfile,outfile_charge,alphaprot) #Makes coord files for matlab and returns a dict for alpha complex                      
    outfile.close()
    outfile_charge.close()


    alpha(name, file_folder, alphaprot)
    Feature_alpha_2D(name, file_folder)

    #rips ./data/D1A2K D1A2K-a0a
    line = "matlab -nodisplay -nodesktop -nosplash -r 'rips "+file_folder +" "+ name +"'"
    print('Starting Matlab...')
    os.system(line)
    
    Feature_rips_1D(name, file_folder)
    Feature_rips_charge_1D(name, file_folder)



    

#matlab -nodisplay -nodesktop -nosplash -r PH_complex_interaction.m hej test

if __name__ == '__main__':
    main()





