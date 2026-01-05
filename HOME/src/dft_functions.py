###################################################################
######################## Standard Libaries ########################
###################################################################
import subprocess
import json
import re

###################################################################
###################### Third-Party Libaries #######################
###################################################################
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

###################################################################
######################## Custom Libaries ##########################
###################################################################
from params import DFT_CPU_NUM, DFT_MEM_HIGH, DFT_MEM_PER_CPU_HIGH, DFT_MEM_LOW, DFT_MEM_PER_CPU_LOW, DFT_RUNTIME_LIMIT, DFT_FUNCTIONAL, PARTITION_NAME
from paths import paths
from atomic_numbers import atomic_symbol_dict, def2_tzvp_basis_functions
from aux_functions import get_comp_time_string, change_status, start_job, get_atom_count_from_smiles, mark_trc_as_failed, move_all_to_failed
from xtb_functions import isotopic_anomaly
from rdk_functions import get_number_of_electrons

###################################################################
########################### FUNCTIONS #############################
###################################################################

def execute_dft_job(ddb_no, smiles,is_retry=False):
    """Execute DFT job for a given DDB number. Has to be run after XTB job is finished."""
    job_type = 'dft'
    cpu_num = DFT_CPU_NUM
    runtime = calc_comp_time_dft(smiles)
    if runtime > DFT_RUNTIME_LIMIT:
        mark_trc_as_failed(ddb_no, ['dft'], 23)
        move_all_to_failed(ddb_no)
    else:
        dft_dir_path = paths['calc'] / f'{ddb_no}' / 'dft'
        dft_job_path = dft_dir_path / 'job'
        change_status(ddb_no,['dft'], ['p'])
        dft_dir_path.mkdir(parents=True, exist_ok=True)
        if is_retry:
            elongate_negative_frequency(ddb_no, Chem.AddHs(Chem.MolFromSmiles(smiles)).GetNumAtoms())
            (paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log').rename(paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_old.log') # rename first try log file
        get_job_dft(ddb_no, smiles, runtime, cpu_num)
        get_com_dft(ddb_no, smiles, cpu_num, is_retry)
        get_py_analyze_dft(ddb_no)
        start_job(dft_job_path, job_type, ddb_no)
        if is_retry:
            change_status(ddb_no, ['dft'], ['r_rt'])

def get_job_dft(ddb_no, smiles, runtime, cpu_num):
    """Create SLURM job script for DFT job."""
    runtime_str = get_comp_time_string(runtime)
    com_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.com'
    job_path = paths['calc'] / f'{ddb_no}' / 'dft' / 'job'
    py_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_analyze_dft.py'
    mem = get_MEM(get_num_basis_functions(smiles))[1]
    name = f'{ddb_no}_dft'
    job_str = f"""#!/usr/bin/bash
#SBATCH --job-name={name}
#SBATCH --account={PARTITION_NAME}
#SBATCH --partition=skylake-96
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpu_num}
#SBATCH --mem-per-cpu={mem}
#SBATCH --time={runtime_str}

module purge
module load gaussian/latest
g16 {com_path}
python {py_path}
"""
    with job_path.open('w') as file:
        file.write(job_str)

def get_com_dft(ddb_no, smiles, cpu_num, is_retry=False):
    """Create Gaussian .com file for DFT job."""
    crest_path = paths['calc'] / f'{ddb_no}' / 'xtb' / 'crest_best.xyz'
    if isotopic_anomaly(smiles) and not is_retry:
        xyz_string = read_xyz_smiles_file_to_string(crest_path, smiles)
    elif not isotopic_anomaly(smiles) and not is_retry:
        xyz_string = read_xyz_file_to_string(crest_path)
    elif isotopic_anomaly(smiles) and is_retry:
        xyz_retry_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_retry.xyz'
        xyz_string = read_xyz_smiles_file_to_string(xyz_retry_path, smiles)
    elif not isotopic_anomaly(smiles) and is_retry:
        xyz_retry_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_retry.xyz'
        xyz_string = read_xyz_file_to_string(xyz_retry_path)
    com_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.com'
    multiplicity = calculate_multiplicity(smiles)
    charge = calculate_charge(smiles)
    islinear = check_linearity(xyz_string)
    mem = get_MEM(get_num_basis_functions(smiles))[0]
    name = f'{ddb_no}_dft'
    if islinear:
        opt_line = f'# opt=cartesian def2TZVP {DFT_FUNCTIONAL} IOp(1/6=100) symm=loose IOp(2/16=3) IOp(2/17=3)'
    else:
        opt_line = f'# opt def2TZVP {DFT_FUNCTIONAL} IOp(1/6=100) symm=loose IOp(2/16=3) IOp(2/17=3)'
    com_str=f"""%Mem={mem}GB
%NProcShared={cpu_num}
%chk={paths['chk'] /name}.chk
{opt_line}

{name}

{charge} {multiplicity}
{xyz_string}

--link1--
%Mem={mem}GB
%NProcShared={cpu_num}
%chk={paths['chk'] /name}.chk
# freq polar geom=Allcheck def2TZVP {DFT_FUNCTIONAL} IOp(1/6=100) symm=loose IOp(2/16=3) IOp(2/17=3)

--link1--
%Mem={mem}GB
%NProcShared={cpu_num}
%chk={paths['chk'] /name}.chk
# nmr=(giao,PrintEigenvectors,Susceptibility) geom=Allcheck def2TZVP {DFT_FUNCTIONAL} IOp(1/6=100) symm=loose IOp(2/16=3) IOp(2/17=3)

--link1--
%Mem={mem}GB
%NProcShared={cpu_num}
%chk={paths['chk'] /name}.chk
#p scrf=COSMORS geom=Allcheck def2TZVP {DFT_FUNCTIONAL}

{name}.cosmo



"""
    
    with com_path.open('w') as file:
        file.write(com_str)

def get_py_analyze_dft(ddb_no):
    """Create Python script to analyze DFT job results for convergence."""
    py_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_analyze_dft.py'
    xyz_retry_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_retry.xyz'
    py_str = f"""#!/usr/bin/python
from pathlib import Path
import sys
sys.path.append('{paths['src']}')
from aux_functions import (change_status, mark_trc_as_failed, move_all_to_failed)
from dft_functions import (check_convergence)

if check_convergence({ddb_no}) and not Path('{xyz_retry_path}').is_file():
    change_status({ddb_no}, ['dft'], ['c'])
elif check_convergence({ddb_no}) and Path('{xyz_retry_path}').is_file():
    change_status({ddb_no}, ['dft'], ['c_rt'])
else:
    mark_trc_as_failed({ddb_no}, ['dft'], 16)
    move_all_to_failed({ddb_no})

"""
    with py_path.open('w') as file:
        file.write(py_str)



def read_xyz_file_to_string(path):
    skip_lines=2
    with path.open('r') as file:
        lines = file.readlines()[skip_lines:]

    xyz = ''.join(lines)

    return xyz

def read_xyz_smiles_file_to_string(path, smiles):
    skip_lines=2
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    isotopes = []
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        symbol = atom.GetSymbol()
        if isotope > 0: 
            isotopes.append(f"{symbol}(iso={isotope})")
        else:
            isotopes.append(symbol)
    
    with path.open('r') as file:
        lines = file.readlines()[skip_lines:]
    
    formatted_xyz = []
    for line, isotope in zip(lines, isotopes):
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f"Ungültiges Format in der XYZ-Datei: {line.strip()}")
        x, y, z = parts[1:4]
        formatted_xyz.append(f"{isotope:<10} {float(x):>12.6f} {float(y):>12.6f} {float(z):>12.6f}")
    
    # Rückgabe als String
    return "\n".join(formatted_xyz)

def calc_comp_time_dft(smiles):
    """Empirical formula to calculate required computation time for DFT job based on number of electrons. Dependend on system."""
    num_elec = get_number_of_electrons(smiles)  
    seconds = 2*int(0.25*num_elec**2.5+600)
    if seconds < 1.5 * 3600:
        return int(1.5 * 3600)
    elif seconds < 48 * 3600:
        return int(48 * 3600)
    elif seconds < 72 * 3600:
        return int(72 * 3600)
    elif seconds > int(10*24*3600):
        return int(10*24*3600)
    else:
        return seconds


def calculate_multiplicity(smiles):
    """Estimation of multiplicity based on SMILES string. Handles O2 as special case. Otherwise assumes low spin states."""
    if smiles == 'O=O':
        multiplicity = 3
        return multiplicity
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        except:
            multiplicity = 1
            return multiplicity
        
        total_valence_electrons = sum(get_valence_electrons(atom) for atom in mol.GetAtoms())

        formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        total_valence_electrons -= formal_charge
        
        unpaired_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        
        if unpaired_electrons > 0:
            multiplicity = unpaired_electrons + 1
        else:
            multiplicity = 1 if total_valence_electrons % 2 == 0 else 2
        
        return multiplicity
    
def calculate_charge(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    
    charge = Chem.GetFormalCharge(mol)
    return charge

def check_linearity(xyz_str):
    atoms, coordinates = read_xyz_for_lin(xyz_str)
    return is_linear_molecule(atoms, coordinates)

def read_xyz_for_lin(xyz_str):
    lines = xyz_str.strip().split('\n')
    atoms = []
    coordinates = []
    for line in lines:
        parts = line.split()
        atoms.append(parts[0])
        coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return atoms, np.array(coordinates)

def is_linear_molecule(atoms, coordinates):
    if len(atoms) < 3:
        return True

    for i in range(1, len(atoms) - 1):
        vec1 = coordinates[i - 1] - coordinates[i]
        vec2 = coordinates[i + 1] - coordinates[i]
        angle = calculate_angle(vec1, vec2)
        if not np.isclose(angle, 180, atol=5):  # Allowing a tolerance of 5 degrees
            return False
    
    return True

def calculate_angle(v1, v2):
    unit_v1 = v1 / np.linalg.norm(v1)
    unit_v2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(unit_v1, unit_v2)
    angle = np.arccos(dot_product)
    return np.degrees(angle)

def get_valence_electrons(atom):
    valence_electrons = {
        1: 1,  2: 2,  3: 1,  4: 2,  5: 3,  6: 4,  # H, He, Li, Be, B, C
        7: 5,  8: 6,  9: 7,  10: 8, 11: 1, 12: 2,  # N, O, F, Ne, Na, Mg
        13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8,  # Al, Si, P, S, Cl, Ar
        19: 1, 20: 2, 21: 3, 22: 4, 23: 5, 24: 6,  # K, Ca, Sc, Ti, V, Cr
        25: 7, 26: 8, 27: 9, 28: 10, 29: 11, 30: 12,  # Mn, Fe, Co, Ni, Cu, Zn
        31: 3, 32: 4, 33: 5, 34: 6, 35: 7, 36: 8,  # Ga, Ge, As, Se, Br, Kr
        37: 1, 38: 2, 39: 3, 40: 4, 41: 5, 42: 6,  # Rb, Sr, Y, Zr, Nb, Mo
        43: 7, 44: 8, 45: 9, 46: 10, 47: 11, 48: 12,  # Tc, Ru, Rh, Pd, Ag, Cd
        49: 3, 50: 4, 51: 5, 52: 6, 53: 7, 54: 8,  # In, Sn, Sb, Te, I, Xe
        55: 1, 56: 2, 57: 3, 58: 4, 59: 5, 60: 6,  # Cs, Ba, La, Ce, Pr, Nd
        61: 7, 62: 8, 63: 9, 64: 10, 65: 11, 66: 12,  # Pm, Sm, Eu, Gd, Tb, Dy
        67: 13, 68: 14, 69: 15, 70: 16, 71: 17, 72: 18,  # Ho, Er, Tm, Yb, Lu, Hf
        73: 19, 74: 20, 75: 21, 76: 22, 77: 23, 78: 24,  # Ta, W, Re, Os, Ir, Pt
        79: 25, 80: 26, 81: 27, 82: 28, 83: 29, 84: 30,  # Au, Hg, Tl, Pb, Bi, Po
        85: 31, 86: 32, 87: 1, 88: 2, 89: 3, 90: 4,  # At, Rn, Fr, Ra, Ac, Th
        91: 5, 92: 6, 93: 7, 94: 8, 95: 9, 96: 10,  # Pa, U, Np, Pu, Am, Cm
        97: 11, 98: 12, 99: 13, 100: 14, 101: 15, 102: 16,  # Bk, Cf, Es, Fm, Md, No
        103: 17, 104: 18, 105: 19, 106: 20, 107: 21, 108: 22,  # Lr, Rf, Db, Sg, Bh, Hs
        109: 23, 110: 24, 111: 25, 112: 26, 113: 27, 114: 28,  # Mt, Ds, Rg, Cn, Nh, Fl
        115: 29, 116: 30, 117: 31, 118: 32  # Mc, Lv, Ts, Og
    }
    return valence_electrons.get(atom.GetAtomicNum(), 0)

def monitor_dft_status(ddb_no):
    trc_file_path = paths['trc_pending'] / f'{ddb_no}.trc'

    with trc_file_path.open('r', encoding='utf-8') as trc_file:
        trc_dict = json.load(trc_file)
        if trc_dict.get('dft') == 'r':
            job_name = f"{ddb_no}_dft"
            result = subprocess.run(
            ['squeue', '--name', job_name, '--noheader'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
            )
            if result.returncode != 0 or not result.stdout.strip():
                change_status(ddb_no, ['dft'], ['f'])


def check_convergence(ddb_no):
    path_calc_dir = paths['calc'] / f'{ddb_no}' / 'dft'

    err_files = [f for f in path_calc_dir.glob("*.err") if f.stat().st_size == 0]

    forbidden_extensions = {".d2e", ".int", ".rwf", ".skr"}
    forbidden_files = [f for f in path_calc_dir.iterdir() if f.suffix in forbidden_extensions]
    
    return bool(err_files) and not forbidden_files

def elongate_negative_frequency(ddb_no, atom_num):
    xyz_file_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_retry.xyz'
    log_file_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'

    atom_symbols, coordinates = read_log_file_to_array(log_file_path)

    first_mode = read_first_mode(log_file_path, atom_num)
    displacement_factor = 1.5

    updated_coordinates = coordinates + (displacement_factor * first_mode)

    comment_line = "Modified geometry for 2"
    xyz_string = f"{len(atom_symbols)}\n{comment_line}\n" + "\n".join(
        f"{symbol:<6} {x:12.6f} {y:12.6f} {z:12.6f}"
        for symbol, (x, y, z) in zip(atom_symbols, updated_coordinates)
    )

    with xyz_file_path.open('w') as f:
        f.write(xyz_string)

def read_com_file_to_string(path, atom_number):
    with path.open('r') as f:
        lines=f.readlines()

    coord_start = None
    for idx, line in enumerate(lines):
        if "0 1" in line:
            coord_start = idx + 1
            break

    coord_end = coord_start + atom_number

    return lines[coord_start:coord_end]

def read_frequencies(path):
    with path.open('r') as f:
        content = f.read()

    # search for all frequencies
    frequencies = re.findall(r"Frequencies --\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", content)

    return frequencies


def read_first_mode(path, n_atoms):
    with path.open('r') as f:
        content = f.read()

    # Find the normal modes section
    mode_match = re.search(r"Atom\s+AN\s+X\s+Y\s+Z([\s\S]+?)(?=\n\n|\Z)", content)
    if not mode_match:
        raise ValueError("No normal modes found.")

    mode_lines = mode_match.group(1).strip().splitlines()


    if len(mode_lines) < n_atoms:
        raise ValueError(f"The number of mode data lines is less than the number of atoms ({n_atoms}).")

    # Extract the first mode
    first_mode = []
    for i in range(len(mode_lines)):
        parts = mode_lines[i].split()
        # Skip invalid lines (e.g., headers like 'X Y Z')
        if len(parts) < 5 or not parts[0].isdigit():
            continue
        try:
            # Extract only the first mode (first three columns X, Y, Z)
            mode_values = list(map(float, parts[2:5]))  # X, Y, Z of the first mode
            first_mode.append(mode_values)
        except ValueError:
            raise ValueError(f"Invalid line while parsing: {mode_lines[i]}")

        # Stop after `n_atoms` lines
        if len(first_mode) == n_atoms:
            break

    first_mode = np.array(first_mode)

    return first_mode


def read_log_file_to_array(log_file_path):
    with log_file_path.open('r') as file:
        lines = file.readlines()

    stationary_found = False
    start_collecting = False
    atom_symbols = []
    coordinates = []

    for line in lines:
        if "Stationary point found." in line and not stationary_found:
            stationary_found = True
            continue

        if stationary_found and "Standard orientation:" in line:
            start_collecting = True
            continue

        if start_collecting:
            if line.strip().startswith("-----"):
                if coordinates:  # Table ended
                    break
                continue

            parts = line.split()
            if len(parts) == 6 and parts[0].isdigit():  # Valid line
                atom_symbols.append(atomic_symbol_dict[int(parts[1])])  # Atom number or symbol
                coordinates.append([float(parts[3]), float(parts[4]), float(parts[5])])

    return np.array(atom_symbols), np.array(coordinates)

def get_num_basis_functions(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    num_basis_functions = 0
    for atom in mol.GetAtoms():
        element = atom.GetSymbol()
        if element not in def2_tzvp_basis_functions:
            return -1
        num_basis_functions += def2_tzvp_basis_functions[element]
    
    return num_basis_functions

def get_MEM(num_basis_functions):
    if num_basis_functions <= 250:
        return DFT_MEM_LOW, DFT_MEM_PER_CPU_LOW
    else:
        return DFT_MEM_HIGH, DFT_MEM_PER_CPU_HIGH
