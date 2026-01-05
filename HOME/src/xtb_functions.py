###################################################################
######################## Standard Libaries ########################
###################################################################
import shutil
import re
import subprocess
import json

###################################################################
###################### Third-Party Libaries #######################
###################################################################
import numpy as np
from rdkit import Chem

###################################################################
######################## Custom Libaries ##########################
###################################################################
from aux_functions import (get_comp_time_string, mark_trc_as_failed, move_all_to_failed, change_status, start_job)
from rdk_functions import get_number_of_electrons
from paths import paths
from atomic_numbers import atomic_number_dict
from params import (XTB_CPU_NUM, EWIN, XTB_MEM, XTB_RUNTIME_LIMIT, PARTITION_NAME, ENVIRONMENT_NAME)

###################################################################
########################### FUNCTIONS #############################
###################################################################

def execute_xtb_job(ddb_no, smiles, is_retry=False):
    """Execute xtb job for a given DDB number and SMILES code."""
    job_type = 'xtb'
    has_iso = False
    cpu_num = XTB_CPU_NUM
    runtime = calc_comp_time_xtb(smiles)
    if runtime > XTB_RUNTIME_LIMIT:
        mark_trc_as_failed(ddb_no, ['xtb'], 22)
        move_all_to_failed(ddb_no)
    else:
        method  = choose_xtb_method(runtime)
        xyz_path = paths['calc'] / f'{ddb_no}' / 'rdk' /  f'{ddb_no}_rdk.xyz'
        xtb_dir_path = paths['calc'] / f'{ddb_no}' / 'xtb'
        xtb_coord_path = paths['calc'] / f'{ddb_no}' / 'xtb' / f'{ddb_no}_rdk.xyz'
        xtb_job_path = xtb_dir_path / 'job'
        change_status(ddb_no,['xtb'], ['p'])

        xtb_dir_path.mkdir(parents=True, exist_ok=True)
        shutil.copy2(xyz_path, xtb_coord_path)
        try:
            if isotopic_anomaly(smiles):
                has_iso = True
                get_isotopic_details(ddb_no, smiles)
        except:
            mark_trc_as_failed(ddb_no, ['xtb'], 21)
            move_all_to_failed(ddb_no)
        get_py_xtb(ddb_no, cpu_num, method, has_iso, is_retry)
        get_job_xtb(ddb_no, runtime, cpu_num)
        start_job(xtb_job_path, job_type, ddb_no)

def get_py_xtb(ddb_no, cpu_num, method, has_iso=False, is_retry=False):
    """Create python script for xtb job."""
    xtb_dir_path = paths['calc'] / f'{ddb_no}' / 'xtb'
    py_path = xtb_dir_path / f'{ddb_no}_xtb.py'
    out_path = xtb_dir_path / f'{ddb_no}_xtb.out'
    err_path = xtb_dir_path / f'{ddb_no}_xtb.err'
    xyz_path = xtb_dir_path / f'{ddb_no}_rdk.xyz'
    ewin = EWIN

    if has_iso and not is_retry:
        cmd_line = f"""result = subprocess.run(["crest", '{xyz_path}', "--T", "{cpu_num}", "--{method}", "--ewin", "{ewin}", "details.inp"], check=True, stdout=out, stderr=err)"""
    elif not has_iso and not is_retry:
        cmd_line = f"""result = subprocess.run(["crest", '{xyz_path}', "--T", "{cpu_num}", "--{method}", "--ewin", "{ewin}"], check=True, stdout=out, stderr=err)"""
    elif not has_iso and is_retry:
        cmd_line = f"""result = subprocess.run(["crest", '{xyz_path}', "--T", "{cpu_num}", "--{method}", "--ewin", "{ewin}", "--noreftopo"], check=True, stdout=out, stderr=err)"""
    elif has_iso and is_retry:
        cmd_line = f"""result = subprocess.run(["crest", '{xyz_path}', "--T", "{cpu_num}", "--{method}", "--ewin", "{ewin}", "details.inp", "--noreftopo"], check=True, stdout=out, stderr=err)"""

    
    commandStr = f"""#!/usr/bin/python
import subprocess
import sys
sys.path.append('{paths['src']}')
from pathlib import Path
from xtb_functions import (check_for_xtb_struct)

try:
    with open(Path('{out_path}'), 'w') as out, open(Path('{err_path}'), 'w') as err:
        {cmd_line}
except subprocess.CalledProcessError as e:
    print(f"Error in xtb calculation: {{e}}")

check_for_xtb_struct({ddb_no})"""

    with py_path.open('w') as file:
        file.write(commandStr)

def get_job_xtb(ddb_no, runtime, cpu_num):
    """Create SLURM job script for xtb job."""
    job_name = f'{ddb_no}_xtb'
    runtime_str = get_comp_time_string(runtime)
    mem = XTB_MEM
    job_path = paths['calc'] / f'{ddb_no}' / 'xtb' / 'job'
    py_path = paths['calc'] / f'{ddb_no}' / 'xtb' / f'{ddb_no}_xtb.py'
    job_str = f"""#!/usr/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={PARTITION_NAME}
#SBATCH --partition=skylake-96
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpu_num}
#SBATCH --mem-per-cpu={mem}
#SBATCH --time={runtime_str}

module purge
module load anaconda3/latest
conda activate {ENVIRONMENT_NAME}

python3 {py_path}

conda deactivate
"""
    with job_path.open('w') as file:
        file.write(job_str)


def calc_comp_time_xtb(smiles):
    """Empiric estimation of the computation time for xtb job based on number of electrons. Dependend on system."""
    num_elec = get_number_of_electrons(smiles)  
    seconds = int(1700*np.exp(0.0301*num_elec)+600)
    if seconds < 1.5 * 3600:
        return int(1.5 * 3600)  # 1.5 hours in seconds
    elif seconds < 48 * 3600:
        return int(48 * 3600)  # 48 hours in seconds
    elif seconds < 72 * 3600:
        return int(72 * 3600)  # 72 hours in seconds
    elif seconds > int(10 * 24 * 3600):
        return int(10 * 24 * 3600)  # 10 days in seconds
    else:
        return seconds 

def get_number_of_electrons(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Invalid SMILES.')
     
    num_elec = 0
    for atom in mol.GetAtoms():
        num_elec += atomic_number_dict[atom.GetSymbol()]

    return num_elec

def choose_xtb_method(runtime):
    "Deprecated, currently only gfn2 is used."
    # if runtime < 36000:
    #     return 'gfn2'
    # elif runtime < 86400:
    #     return 'gfn2//gfnff'
    # else:
    #     return 'gfnff'
    return 'gfn2'
    
def check_for_xtb_struct(ddb_no):
    xyz_path = paths['calc'] / f'{ddb_no}' / 'xtb' / 'crest_best.xyz'
    xtb_out_path = paths['calc'] / f'{ddb_no}' / 'xtb' / f'{ddb_no}_xtb.out'

    with xtb_out_path.open('r') as file:
        xtb_out_lst = file.readlines()
    
    if xyz_path.is_file() and ('CREST terminated normally.' in xtb_out_lst[-1]):
        change_status(ddb_no, ['xtb'], ['s'])
    elif has_topo_error(xtb_out_path):
        change_status(ddb_no, ['xtb'], ['wrt']) # Waiting for ReTry (wrt)    
    else:
        mark_trc_as_failed(ddb_no, ['xtb'], 6)
        move_all_to_failed(ddb_no)

def has_topo_error(xtb_out_path):
    with xtb_out_path.open('r') as file:
        for line in file:
            if "WARNING! Change in topology detected!" in line:
                return True
    return False


def isotopic_anomaly(smiles):
    isotope_pattern = re.compile(r"\[\d+[A-Za-z]+\]")
    return bool(isotope_pattern.search(smiles))

def get_isotopic_details(ddb_no, smiles):
    detail_file_path = paths['calc'] / f'{ddb_no}' / 'xtb' / 'details.inp'
    mol = Chem.MolFromSmiles(smiles)    
    isotopes = {}
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope() 
        if isotope > 0:  # Only isotopes with different mass number
            symbol = atom.GetSymbol()
            isotope_symbol = f"{isotope}{symbol}"  # 13C, 2H etc.
            isotopes[isotope_symbol] = isotope * Chem.GetPeriodicTable().GetAtomicWeight(symbol)

    with open(detail_file_path, "w") as f:
        f.write("$hess\n")
        f.write("  element mass:\n")
        for symbol, mass in isotopes.items():
            f.write(f"    {symbol}, {mass:.5f}\n")
        f.write("$end\n")

def monitor_xtb_status(ddb_no):
    trc_file_path = paths['trc_pending'] / f'{ddb_no}.trc'

    with trc_file_path.open('r', encoding='utf-8') as trc_file:
        trc_dict = json.load(trc_file)
        if trc_dict.get('xtb') == 'r':
            job_name = f"{ddb_no}_xtb"
            result = subprocess.run(
            ['squeue', '--name', job_name, '--noheader'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
            )
            if result.returncode != 0 or not result.stdout.strip():
                change_status(ddb_no, ['xtb'], ['f'])