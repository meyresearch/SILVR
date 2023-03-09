#!/usr/bin/env python
# coding: utf-8

# # Run below for all molecule 

# In[1]:


"""
#-------old code--------------

import os
import glob

exp_list = glob.glob("../experiments/exp_*")


if not os.path.exists(f'chemgauss_docking/'):
    os.makedirs(f'chemgauss_docking/')
        
        
for exp in exp_list:
    exp_id = exp.split("/")[-1]
    
    if not os.path.exists(f'chemgauss_docking/{exp_id}/'):
        os.makedirs(f'chemgauss_docking/{exp_id}/')
        
        
    
    ligands = glob.glob(f"../experiments/{exp_id}/*.txt")

    for ligand in ligands:
        mol_id = ligand.split("/")[-1].split(".")[0]

        with open(ligand, "r") as readfile:
            file = readfile.read()

        with open("tmp.xyz","w") as writefile:
            writefile.write(file)

        #os.system(f'python DockMolecules.py -in tmp.xyz -out chemgauss_docking/{exp_id}/{mol_id}.sdf -receptor mpro2.oedu')
        
"""


# In[2]:


import os
import glob

#------New code-------
#This could probably be set to run in parallel
#Only molecules that are not fragmented get docked

exp_list = glob.glob("../experiments/exp_*")

if not os.path.exists(f'chemgauss_docking_cleaned/'):
    os.makedirs(f'chemgauss_docking_cleaned/')
        
for exp in exp_list:
    exp_id = exp.split("/")[-1]
    
    if not os.path.exists(f'chemgauss_docking_cleaned/{exp_id}/'):
        os.makedirs(f'chemgauss_docking_cleaned/{exp_id}/')
        
        
    
    ligands = glob.glob(f"../experiments/{exp_id}/*.txt")

    for ligand in ligands:
        mol_id = ligand.split("/")[-1].split(".")[0]

        with open(ligand, "r") as readfile:
            file = readfile.read()

        with open("tmp.xyz","w") as writefile:
            writefile.write(file)

        os.system(f'python CleanThenDockMolecules.py -in tmp.xyz -out chemgauss_docking_cleaned/{exp_id}/{mol_id}.sdf -receptor mpro2.oedu')


# In[136]:





# # For visualisation only

# In[3]:


"""
import py3Dmol

def show_mol_file(mol_file):
    view = py3Dmol.view(query=mol_file)  
    view.setStyle({'stick': {'color':'spectrum'}})
    
    
    with open(mol_file) as readfile:
        try:
            print("score: ", readfile.readlines()[-3])
        except:
            print("Molecule was fragmented")
    
    
    return view.show()
"""


# In[4]:


#all_docks = glob.glob("chemgauss_docking_cleaned/exp_5/*.sdf")
#original_docking = ["chemgauss_docking/exp_5/"+x.split("/")[-1].split(".")[0]+".sdf" for x in all_docks]
#new_docking = ["chemgauss_docking_cleaned/exp_5/"+x.split("/")[-1].split(".")[0]+".sdf" for x in all_docks]
#original_xyz = ["../experiments/exp_5/"+x.split("/")[-1].split(".")[0]+".txt" for x in all_docks]


# In[5]:


#idx = 9


# In[6]:


#Original xyz
#show_mol_file(original_xyz[idx])


# In[7]:


#new docking
#show_mol_file(new_docking[idx])


# In[8]:


#Original docking
#show_mol_file(original_docking[idx])


# In[ ]:





# In[ ]:




