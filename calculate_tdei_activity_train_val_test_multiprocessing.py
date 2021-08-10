# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:27:36 2021

@author: Gil

Calculate the number of interaction between electrons

      target atom EC
         _ _ _ _
        |       |
source  |       |
atom    |       |
EC      |_______|

1) 0D: atom itself
2) 1D: atom-atoms in 1D topological distance
3) 2D: atom-atoms in 2D 
4) 3D: atom-atoms in 3D
"""

from rdkit import Chem
from multiprocessing import Pool

import pandas as pd
import numpy as np
import h5py
import os

#Step1: load electron configuration (EC) vector file into 'excel' variable. (Figure 1)
excel = os.environ['CSV']
print(excel)
ec = pd.read_csv(excel)
col_selec = ec.columns[1:]

dim = int(os.environ['DIM'])
activity = os.environ['ACTIVITY']
smi = os.environ['SMI']

topological_distance = list(range(dim))

feature_dim = len(col_selec)
distance_dim = len(topological_distance)

#Step 2: electron interaction (Ei) matrix calculation. (Figure 2)
def calculate_elec_int_map (vector1, vector2, atom_tensor):
    for ir in range(len(vector1)):
        for ic in range(len(vector2)):
            atom_tensor[ir, ic] = vector1[ir]*vector2[ic]
    return

#Step 3: TDEi tensor was compiled by collecting Ei matrices
#        in different topological distances.(Figure 3 & 4)
def calculate_topological_tensor(each_dict, atom_tensor, mol_tensor):
    #First element of the list is 0D atom's electron configuration vector
    source_ec = each_dict[0][0] 
    for key, value in each_dict.items():
        for v in value:
            #print (key, source_ec, type(source_ec),v, type(v))
            calculate_elec_int_map(source_ec, v, atom_tensor)
            mol_tensor[:,:,key] = np.add(mol_tensor[:,:,key], atom_tensor)
    return

#For loop to calculate TDEi tensor on multiple molecules
def calculate_tdim_tensor (df):
    tensor_dict = {'indice':[],
                   activity:[],
                   'tec':[]}
    
    atom_tensor = np.zeros((feature_dim, feature_dim))

    for index, row in df.iterrows():
        #print ('Started:',index)
        mol = Chem.MolFromSmiles(row[smi])
        
        #Hydrogen must be added to consider intrinsic hydrgeon in TDEi tensor calculation.
        molH = Chem.AddHs(mol)

        TopDistMat = Chem.rdmolops.GetDistanceMatrix(molH)
        
        mol_tensor = np.zeros((feature_dim, feature_dim, distance_dim))
        #check_flag = 0
        for vector in TopDistMat:
            top_dict = {}
            for t in topological_distance:
                top_dict[t]=[]
            for t in topological_distance:
                atom_index = np.where(vector == t)[0]
                for a in atom_index:
                    atom_symbol = molH.GetAtomWithIdx(int(a)).GetSymbol()

                    #Collect EC vector from atoms at topological distance 't'.
                    ec_atom = ec.loc[ec['atom']==atom_symbol,col_selec].values[0]
                    top_dict[t].append(ec_atom)
                    
            calculate_topological_tensor(top_dict, atom_tensor, mol_tensor)
            #check_flag += top_dict[0][0][0]**2
        
        for t in topological_distance[1:]:
            mol_tensor[:,:,t] = mol_tensor[:,:,t]/2
        
        tensor_dict['indice'].append(index)
        tensor_dict[activity].append(row[activity])
        tensor_dict['tec'].append(mol_tensor)
        #print ('Done:',index)
    
    return tensor_dict

#Multiprocessing was essential to calculate over more than 200,000 molecules.
def parallelize_tec_calculation(df, func, num_cores=4):
    df_split = np.array_split(df, num_cores)
    pool = Pool(num_cores)
    tensor_dict = pool.map(func, df_split)
    pool.close()
    pool.join()
    return tensor_dict

cpu_core = int(os.environ['CPU'])

input_dir = os.environ['DIR']
csv_files = os.listdir(input_dir)
#csv_files = ['mp_test.csv', 'mp_val.csv', 'mp_train.csv']

#For loop to calculate TDEi tensor over multiple files.
for csv in csv_files:
    if csv.endswith('csv'):
        outfile = excel.replace('configuration','interaction_'+str(dim)+'d_'+csv.split('.')[0]).replace('csv','h5')
        if os.path.isfile(os.path.join(input_dir, outfile)):
            print ('File exists:',outfile)
        else:
            print (csv)
            #df = pd.read_csv(csv, dtype={'smi':str,'mp':float})
            df = pd.read_csv(os.path.join(input_dir,csv))
            
            #Topological electron configuration tensor calculation
            
            tec_dict = parallelize_tec_calculation(df, calculate_tdim_tensor, cpu_core)
            
            print ('Combine TEC tensors from each core.')
            final_dict = {'indice':[],
                          activity:[],
                          'tec':[]}
            for t in tec_dict:
                for key, value in t.items():
                    if key=='indice':
                        final_dict['indice']+=value
                    elif key==activity:
                        final_dict[activity]+=value
                    elif key=='tec':
                        final_dict['tec']+=value
                    else:
                        print (t)
            #print (final_dict)        
                        
            #save file
            print ('Input save in h5 format.')
            
            hf = h5py.File(os.path.join(input_dir, outfile),'w')
            hf.create_dataset('indice',data=final_dict['indice'])
            hf.create_dataset(activity.upper(), data=final_dict[activity])
            hf.create_dataset('TEC', data=final_dict['tec'])
            hf.close()
            print ('*'*10)
