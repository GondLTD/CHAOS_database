###################################################################
######################## Standard Libaries ########################
###################################################################
import json
import subprocess
###################################################################
###################### Third-Party Libaries #######################
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import CalcMolDescriptors
import numpy as np

###################################################################
######################## Custom Libaries ##########################
###################################################################
from atomic_numbers import atomic_number_dict
from paths import paths
from params import (RDKIT_CPU_NUM, RDKIT_RUNTIME, PARTITION_NAME, ENVIRONMENT_NAME)
from aux_functions import (get_comp_time_string, mark_trc_as_failed, move_all_to_failed, change_status, start_job)
from ana_functions import (load_result_json, save_result_json)

###################################################################
########################### FUNCTIONS #############################
###################################################################

def execute_rdk_job(ddb_no, smiles):
    """Execute RDKit job for a given DDB number and SMILES code."""
    if not isinstance(smiles, str):
        mark_trc_as_failed(ddb_no, ['rdk'], 20)
        move_all_to_failed(ddb_no)
    elif is_only_atom(smiles):
        mark_trc_as_failed(ddb_no, ['rdk'], 0)
        move_all_to_failed(ddb_no)
    elif not can_add_Hs(smiles):
        mark_trc_as_failed(ddb_no, ['rdk'], 1)
        move_all_to_failed(ddb_no)
    elif is_ionic(smiles):
        mark_trc_as_failed(ddb_no, ['rdk'], 19)
        move_all_to_failed(ddb_no)
    else:
        job_type = 'rdk'
        cpu_num = RDKIT_CPU_NUM
        runtime = RDKIT_RUNTIME
        dir_path = paths['calc'] / f'{ddb_no}' / 'rdk'
        job_path = dir_path / 'job'
        change_status(ddb_no,['rdk'], ['p'])
        dir_path.mkdir(parents=True, exist_ok=True)
        get_py_rdk(ddb_no, smiles)
        get_job_rdk(ddb_no, runtime, cpu_num)
        start_job(job_path, job_type, ddb_no)

def get_py_rdk(ddb_no, smiles):
    """Create RDKit python script for a given DDB number and SMILES code."""
    xyz_name = f'{ddb_no}_rdk'
    rdk_path = paths['calc'] / f'{ddb_no}' / 'rdk'
    xyz_path = rdk_path / f'{ddb_no}_rdk.xyz'
    py_path = rdk_path / f'{ddb_no}_rdk.py'
    py_str = f"""#!/usr/bin/python
import sys
import pickle
from pathlib import Path
sys.path.append('{paths['src']}')
from rdk_functions import (get_embedded_conformer, check_for_rdk_struct, write_conf_to_xyz_file, extract_rdk_desc)
from aux_functions import (mark_trc_as_failed, move_all_to_failed)
try:
    mol, energy, id = get_embedded_conformer('{smiles}')

    write_conf_to_xyz_file(mol, '{xyz_name}', id, Path('{xyz_path}'))
    check_for_rdk_struct({ddb_no})

except:
    mark_trc_as_failed({ddb_no}, ['rdk'], 4)
    move_all_to_failed({ddb_no})

try:
    extract_rdk_desc(mol, id, {ddb_no})
except:
    mark_trc_as_failed({ddb_no}, ['rdk'], 14)
    move_all_to_failed({ddb_no})
"""
    with py_path.open('w') as file:
        file.write(py_str)

def get_job_rdk(ddb_no, runtime, cpu_num):
    py_job_file_path = paths['calc'] / f'{ddb_no}' / 'rdk' / f'{ddb_no}_rdk.py'
    job_path = paths['calc'] / f'{ddb_no}' / 'rdk' / 'job'
    job_name = f'{ddb_no}_rdk'
    runtime_str = get_comp_time_string(runtime)
    job_str = f"""#!/usr/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={PARTITION_NAME}
#SBATCH --partition=skylake-96
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpu_num}
#SBATCH --time={runtime_str}

module purge
module load anaconda3/latest
conda activate {ENVIRONMENT_NAME}

python3 {py_job_file_path}

conda deactivate
"""
    with job_path.open('w', encoding='utf-8') as file:
        file.write(job_str)

def get_embedded_conformer(smiles):
    def has_overlapping_atoms(mol, conf_id, threshold=1e-3):
        """ Check whether two atoms are too close to each other in a molecule.
        
        Args:
            mol: RDKit Molecule with embedded conformer.
            conf_id: ID of the conformer to check.
            threshold: Tolerance value to consider coordinates as equal (in Angstrom).
            
        Returns:
            bool: True, if two atoms are within the threshold distance, otherwise False.
        """
        conformer = mol.GetConformer(conf_id)
        positions = []
        for atom_idx in range(mol.GetNumAtoms()):
            pos = conformer.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])
        
        positions = np.array(positions)
        
        # Calculate pairwise distances
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist < threshold:
                    return True  # Two atoms are too close to each other
        
        return False
    
    try:
        struct_from_smiles, energy_from_smiles, struct_id = load_best_smiles(smiles)
    except:
        struct_from_smiles, energy_from_smiles, struct_id = None, None, None
    
    # Check for overlapping atoms in the conformer
    if struct_from_smiles and has_overlapping_atoms(struct_from_smiles, struct_id):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol) == 0:  # Success embedding the molecule
            struct_from_smiles = mol
            struct_id = 0
            energy_from_smiles = None  # Energy might need recalculation or updating if needed
    
    return struct_from_smiles, energy_from_smiles, struct_id

def get_embedded_conformer_old(smiles):
    try:
        struct_from_smiles, energy_from_smiles, struct_id = load_best_smiles(smiles)
    except:
        struct_from_smiles, energy_from_smiles, struct_id = None, None, None

    return struct_from_smiles, energy_from_smiles, struct_id

def load_best_smiles(smiles):
    mol = convert_smiles_to_conformers(smiles)
    mol,valid_conf_ids = remove_similiar_conformers(mol)
    mol = optimize_conformers(mol, valid_conf_ids=valid_conf_ids)
    mol, valid_conf_ids = remove_similiar_conformers_after_opt(mol, valid_conf_ids)
    struct, struct_energy, struct_id = lowest_conf_to_mol(mol,valid_conf_ids=valid_conf_ids)

    return struct, struct_energy, struct_id

def convert_smiles_to_conformers(smiles, numConfs = 300):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    try:
        mol = Chem.AddHs(mol)
    except:
        return None
    
    try:
        AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, randomSeed=1050)
    except ValueError as e:
        return None
    
    return mol

def remove_similiar_conformers(mol, threshold=1.0):
    energies = []
    valid_conf_ids  = []
    # Calculate the energy of all conformers
    for conf_id in range(mol.GetNumConformers()):
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if ff is not None:
            energy = ff.CalcEnergy()
            energies.append((conf_id, energy))
    
    # Find energetical unique conformers
    energies.sort(key=lambda x: x[1])
    unique_energies = [energies[0]]
    valid_conf_ids.append(energies[0][0])
    for idx in range(1, len(energies)):
        if energies[idx][1] - unique_energies[-1][1] > threshold:
            unique_energies.append(energies[idx])
            valid_conf_ids.append(energies[idx][0])

    # Remove non unique conformers
    for conf_id, _ in energies:
        if conf_id not in [ue[0] for ue in unique_energies]:
            mol.RemoveConformer(conf_id)
    return mol, valid_conf_ids

def optimize_conformers(mol, valid_conf_ids):
    if mol is None:
        return None

    for conf_id in valid_conf_ids:  # This ensures you're only accessing valid conformers
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if ff is not None:
            ff.Minimize()
    return mol

def remove_similiar_conformers_after_opt(mol, valid_conf_ids, threshold=1.0):
    energies = []
    valid_conf_ids_new  = []
    # Calculate the energy of all conformers
    for conf_id in valid_conf_ids:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if ff is not None:
            energy = ff.CalcEnergy()
            energies.append((conf_id, energy))
    
    # Find energetical unique conformers
    energies.sort(key=lambda x: x[1])
    unique_energies = [energies[0]]
    valid_conf_ids_new.append(energies[0][0])
    for idx in range(1, len(energies)):
        if energies[idx][1] - unique_energies[-1][1] > threshold:
            unique_energies.append(energies[idx])
            valid_conf_ids_new.append(energies[idx][0])

    # Remove non unique conformers
    for conf_id, _ in energies:
        if conf_id not in [ue[0] for ue in unique_energies]:
            mol.RemoveConformer(conf_id)
    return mol, valid_conf_ids_new

def lowest_conf_to_mol(mol, valid_conf_ids):
    lowest_energy = float('inf')
    lowest_conf_id = None

    for conf_id in valid_conf_ids:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if ff is not None:
            energy = ff.CalcEnergy()
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_conf_id = conf_id
    if lowest_conf_id is None:
        return None
    
    for conf_id in list(mol.GetConformers()):
        if conf_id.GetId() != lowest_conf_id:
            mol.RemoveConformer(conf_id.GetId())

    return mol, lowest_energy, lowest_conf_id


def check_for_rdk_struct(ddb_no):
    xyz_file_path = paths['calc'] / f'{ddb_no}' / 'rdk' / f'{ddb_no}_rdk.xyz'

    if xyz_file_path.is_file():
        change_status(ddb_no, ['rdk'], ['s'])
    else:
        mark_trc_as_failed(ddb_no, ['rdk'], 5)
        move_all_to_failed(ddb_no)

def get_number_of_electrons(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Invalid SMILES.')
     
    num_elec = 0
    for atom in mol.GetAtoms():
        num_elec += atomic_number_dict[atom.GetSymbol()]

    return num_elec

def is_only_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    try:
        mol = Chem.AddHs(mol)
    except:
        pass
    return mol.GetNumAtoms() == 1

def can_add_Hs(smiles):
    mol = Chem.MolFromSmiles(smiles)
    try:
        mol = Chem.AddHs(mol)
        return True
    except:
        return False
    
def write_conf_to_xyz_file(mol, name, conf_id, path):
    xyz = f"{mol.GetNumAtoms()}\n{name}\n"
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        symbol = atom.GetSymbol()

        if isotope > 0: 
            symbol = f"{isotope}{symbol}"

        pos = mol.GetConformer(conf_id).GetAtomPosition(atom.GetIdx())
        xyz += f"{symbol:<4} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n"
    
    with open(path, 'w') as file:
        file.write(xyz)

def write_conf_to_xyz_file_old(mol, name, conf_id, path):
    xyz = f"{mol.GetNumAtoms()}\n{name}\n"
    for atom in mol.GetAtoms():
        pos = mol.GetConformer(conf_id).GetAtomPosition(atom.GetIdx())
        xyz += f"{atom.GetSymbol():<2} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n"
    
    with open (path, 'w') as file:
        file.write(xyz)

def extract_rdk_desc(mol, id, ddb_no):
    mol_conf = Chem.Mol(mol)
    mol_conf.RemoveAllConformers()
    mol_conf.AddConformer(mol.GetConformer(id), assignId=True)
    json_dict = load_result_json(ddb_no)
    json_dict['rdk_desc'] = CalcMolDescriptors(mol_conf)
    save_result_json(json_dict, ddb_no)


def is_ionic(smiles):
    return '.' in smiles


def monitor_rdk_status(ddb_no):
    trc_file_path = paths['trc_pending'] / f'{ddb_no}.trc'

    with trc_file_path.open('r', encoding='utf-8') as trc_file:
        trc_dict = json.load(trc_file)
        if trc_dict.get('rdk') == 'r':
            job_name = f"{ddb_no}_rdk"
            result = subprocess.run(
            ['squeue', '--name', job_name, '--noheader'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
            )
            if result.returncode != 0 or not result.stdout.strip():
                change_status(ddb_no, ['rdk'], ['f'])
