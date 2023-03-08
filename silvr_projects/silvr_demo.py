"""
The maintainability of this code is poor and shoud be improved.
This file was not created with the entent to be read by anyone else
"""


import sys
sys.path
_ = (sys.path.append("/usr/local/lib/python3.8/site-packages"))
_ = (sys.path.append("e3_diffusion_for_molecules"))

import eval_sample

# Rdkit import should be first, do not move it
try:
    from rdkit import Chem
except ModuleNotFoundError:
    pass

import utils
import argparse
from configs.datasets_config import qm9_with_h, qm9_without_h
from qm9 import dataset
from qm9.models import get_model

from equivariant_diffusion.utils import assert_correctly_masked
import torch
import pickle
import qm9.visualizer as vis
from qm9.analyze import check_stability
from os.path import join
from qm9.sampling import sample_chain, sample
from configs.datasets_config import get_dataset_info


################


import torch
from torch.distributions.categorical import Categorical

import numpy as np
from egnn.models import EGNN_dynamics_QM9

from equivariant_diffusion.en_diffusion import EnVariationalDiffusion



####################


import numpy as np
import torch
import torch.nn.functional as F
from equivariant_diffusion.utils import assert_mean_zero_with_mask, remove_mean_with_mask,\
    assert_correctly_masked
from qm9.analyze import check_stability


####################
from qm9.analyze import check_stability
from argparse import Namespace
import datetime
import os


####################
from multiprocessing import Process
import time
import math
import datetime




def read_sxyz(file):
    """
    Read reference.sxyz file and parse into correct format for silvr
    
    ---explanation of parsing---
    Line[0] contains total atoms in the reference. This is the same as for an xyz file
    line[1] contains experimental set up params
        dummy: Number of additional atoms to add that are not contained within the reference
        sample: how many molecules should be generated. Molecules are generated serially, so samples
                defines iterations in a for loop
                nb: future work should look at parallelisation
        comment: This does not affect model running
        
    lines[2:] contain all xyz and silvr rate information: element X Y Z rate
        element is converted to float of atomic number
        X Y Z are retained
        rate is appended to the SILVR vector 
        
    returns
        ref_coords (2D array): atomic number and reference coordinates
        total_atoms (int): total number of atoms to be in each sampled molecule (ie 40 = 40 atoms in generated molecule)
        n_samples: Number of samples to take
        silvr_vector (1D array): containing rates to be used in SILVR protocol
            where values are 0 this means no weighting applied to atom
            where values are 1 this mean total weighting for atom, essentially total replacement of generated 
                atoms with reference atom
            generally rate=0.01 works well
    """
    n_ref = 0
    n_dummy = 0
    n_samples = 0
    ref_coords = []
    silvr_vector = []
    
    txt = file.split("\n")
    n_ref = int(txt[0])
    n_dummy = int(txt[1].split()[0].split(":")[1])
    n_samples = int(txt[1].split()[1].split(":")[1])

    #atom mapping 
    atomic_number_map = {"H":1,"B":5,"C":6,"N":7,"O":8,"F":9,"P":15,"S":16,"Cl":17,"br":35}

    for line in txt[2:]:
        line = line.split()
        xyz = [float(x) for x in line[1:4]]
        atomic_number = atomic_number_map[line[0]]            
        ref_coords.append([atomic_number]+xyz)
        silvr_vector.append(float(line[4]))
            
    total_atoms = n_ref + n_dummy
            
    return ref_coords, total_atoms, n_samples, silvr_vector




#--------------SILVR-------------------
class silvr():
    def __init__(self, 
                args=None,
                device=None,
                flow=None,
                dataset_info=None,
                total_atoms=None,
                silvr_rate=0.01,
                ref_coords=None,
                n_samples=1,
                out_dir_path="outputs/",
                gpu=None):
    
        if total_atoms < len(ref_coords):
            #Total atoms must be at least the number as in the reference
            total_atoms = len(ref_coords)

        self.args = args
        self.device = device
        self.flow = flow
        self.dataset_info = dataset_info
        self.max_n_nodes = self.dataset_info['max_n_nodes']
        self.total_atoms = total_atoms
        self.node_mask, self.edge_mask = self.make_node_mask()
        self.ref_coords = ref_coords
        self.ref_node_mask = self.set_ref_node_mask()   
        self.n_samples = n_samples
        self.gpu = gpu
        
        #if path doesn't exist, make path
        self.out_dir_path = out_dir_path
        if not os.path.exists(self.out_dir_path):
            os.makedirs(self.out_dir_path)
            
        #Convert silvr_rate list to tensor of correct shape
        if isinstance(silvr_rate, list):
            silvr_tensor = torch.zeros(self.max_n_nodes).to(device=self.device)
            silvr_tensor[:len(silvr_rate)] = torch.Tensor(silvr_rate)
            silvr_tensor = torch.unsqueeze(silvr_tensor, 0)
            silvr_tensor = torch.unsqueeze(silvr_tensor, 2)
        self.silvr_rate = silvr_tensor
            

    def make_node_mask(self, number=None):
        if number:
            nodesxsample = torch.tensor([atomic_numbers])
        else:
            nodesxsample = torch.tensor([self.total_atoms])
            max_n_nodes = self.dataset_info['max_n_nodes']  # this is the maximum node_size in QM9

            assert int(torch.max(nodesxsample)) <= max_n_nodes
            batch_size = len(nodesxsample)

            node_mask = torch.zeros(batch_size, max_n_nodes)
            for i in range(batch_size):
                node_mask[i, 0:nodesxsample[i]] = 1

        # Compute edge_mask
        edge_mask = node_mask.unsqueeze(1) * node_mask.unsqueeze(2)
        diag_mask = ~torch.eye(edge_mask.size(1), dtype=torch.bool).unsqueeze(0)
        edge_mask *= diag_mask
        edge_mask = edge_mask.view(batch_size * max_n_nodes * max_n_nodes, 1).to(self.device)
        node_mask = node_mask.unsqueeze(2).to(self.device)

        self.edge_mask = edge_mask
        self.node_mask = node_mask

        return node_mask, edge_mask
    
    
    def set_ref_node_mask(self, number=None):
        if number:
            nodesxsample2=torch.tensor([number])
        elif self.ref_coords:
            nodesxsample2=torch.tensor([len(self.ref_coords)])
        else:
            nodesxsample2=torch.tensor([self.total_atoms])

        batch_size2 = len(nodesxsample2)
        node_mask2 = torch.zeros(batch_size2, self.max_n_nodes)
        for i in range(batch_size2):
            node_mask2[i, 0:nodesxsample2[i]] = 1
        self.ref_node_mask = node_mask2.unsqueeze(2).to(self.device)

        return self.ref_node_mask
    
    

    def sample(self):
        """
        Samples diffusion model with SILVR
        Returns one_hot encoding, charges, coordinates, and atom mask
        These can be converted into an XYZ file
        """
        context = None
        x, h = self.flow.sample(self.n_samples, 
                                   self.max_n_nodes, 
                                   self.node_mask, 
                                   self.edge_mask, 
                                   context, 
                                   fix_noise=False, 
                                   silvr_rate=self.silvr_rate,
                                   ref_coords=self.ref_coords,
                                   ref_node_mask=self.ref_node_mask,
                                    shift_centre=True,
                                  dataset_info=self.dataset_info
                            )

        assert_correctly_masked(x, self.node_mask)
        one_hot = h['categorical']
        charges = h['integer']
        assert_correctly_masked(one_hot.float(), self.node_mask)
        if self.args.include_charges:
            assert_correctly_masked(charges.float(), self.node_mask)

        return one_hot, charges, x, self.node_mask



    def sample_xyz(self):
        """
        Samples the SILVR model
        Outputs XYZ file with stability information
        Returns id of file
            Note files are saves as mol_{mol_id}_000.txt
        """
        #run sampling
        one_hot, charges, x, node_mask = self.sample()
        now = datetime.datetime.now()
        mol_id = now.strftime("%Y_%m_%d_%H%M%S")
        
        #Convert sample to xyz file
        vis.save_xyz_file(
          self.out_dir_path, one_hot, charges, x,
          id_from=0, name=f'demo_mol', dataset_info=self.dataset_info,
          node_mask=node_mask)
        
        #Add stability information to xyz file
        i = 0
        num_atoms = int(node_mask[i:i+1].sum().item())
        atom_type = one_hot[i:i+1, :num_atoms].argmax(2).squeeze(0).cpu().detach().numpy()
        x_squeeze = x[i:i+1, :num_atoms].squeeze(0).cpu().detach().numpy()
        mol_stable = check_stability(x_squeeze, atom_type, self.dataset_info)
        stability_string = f"stable:{mol_stable[0]} satoms:{mol_stable[1]} tatoms:{mol_stable[2]} sratio:{mol_stable[1]/mol_stable[2]}\n"

        with open(f'{self.out_dir_path}demo_mol_000.txt', "r") as readfile:
            data = readfile.readlines()
        readfile.close()
        data[1] = stability_string
        with open(f'{self.out_dir_path}demo_mol_000.xyz', "w") as writefile:
            writefile.writelines(data)
        writefile.close()

        #return name of molecule file
        return(mol_id)


    
    
    
def sample_silvr(reference_file):
    
    ref_coords, total_atoms, n_samples, silvr_vector = read_sxyz(reference_file)
    
    
    #1 Model imports
    #--------------Model path--------------------
    model_path = "e3_diffusion_for_molecules/outputs/edm_geom_drugs/"

    #------------------Get model arguments-------------------
    assert model_path is not None
    with open(join(model_path, 'args.pickle'), 'rb') as f:
        args = pickle.load(f)

    #------------Don't know what this is---------------------
    # CAREFUL with this -->
    if not hasattr(args, 'normalization_factor'):
        args.normalization_factor = 1
    if not hasattr(args, 'aggregation_method'):
        args.aggregation_method = 'sum'

    #-------------Setting Cuda-----------------------
    #Limited CUDA?
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f"cuda" if args.cuda else "cpu")
    args.device = device
    dtype = torch.float32
    utils.create_folders(args)

    #-------------Loading arguments needed for model----------------------
    from configs.datasets_config import geom_with_h, geom_no_h
    dataset_info = get_dataset_info(args.dataset, args.remove_h)

    #---------Getting Model-------------
    flow, nodes_dist, prop_dist = get_model(args, device, dataset_info, None)
    flow.to(device)

    #-----------getting additional params?--------------
    fn = 'generative_model_ema.npy' if args.ema_decay > 0 else 'generative_model.npy'
    flow_state_dict = torch.load(join(model_path, fn),map_location=device)
    flow.load_state_dict(flow_state_dict)
    
    
    #2 Create SILVR object
    silvr_model = silvr(args=args,device=device,flow=flow,
                   dataset_info=dataset_info,total_atoms=total_atoms,
                   silvr_rate=silvr_vector,ref_coords=ref_coords,
                   n_samples=1, out_dir_path="demo/", gpu="0")#Output sample in demo/
    
    t1 = time.time()
    silvr_model.sample_xyz()
    t2 = time.time()
    
    print("Total sampling time: ", t2-t1, " seconds")
        
    return



from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from openbabel import pybel
from rdkit import Chem


def rdkit_fix_radicals(mol,add_h=False,flatten=False,uncharge=True):
    """
    Atoms with unfilled valance get radicals assigned.
    Openbabel will have assigned bond orders based on bond length.
    Here I assume all radical electrons should instead be hydrogen atoms
    """
    for atom in mol.GetAtoms():
        radicals = atom.GetNumRadicalElectrons()
        atom.SetNumRadicalElectrons(0)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + radicals)

    if flatten:
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    if add_h:
        mol = Chem.AddHs(mol,addCoords=True)

    if uncharge:
        un = rdMolStandardize.Uncharger()
        mol = un.uncharge(mol)

    return mol


def xyz_to_mol_clean(xyz, add_h=True, flatten=False):
    """
    add_h - add RDKit hydrogens
    flatten - run Chem.MolFromSmiles(Chem.MolToSmiles(x)) such that geometry infromation is lost
    
    Sometimes these imports fail
    In these cases this function returns False
    """
    try:
        mol_pybel = pybel.readstring("xyz", xyz)
        mol_mol2 = mol_pybel.write("mol2")

        #RDKit - clean radicals
        mol_rdkit = Chem.MolFromMol2Block(mol_mol2)
        mol_final = rdkit_fix_radicals(mol_rdkit, add_h=add_h, flatten=flatten)

        return mol_final
    
    except:
        return False
    
def get_mol():
    with open("demo/demo_mol_000.xyz", "r") as readfile:
        xyz = readfile.read()
        mol = xyz_to_mol_clean(xyz)
        
    if not mol:
        raise Exception("SILVR failed. Please run another sample.")
    return mol