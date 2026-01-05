###################################################################
######################## Standard Libaries ########################
###################################################################
import pickle
import json
import re
###################################################################
###################### Third-Party Libaries #######################

###################################################################
######################## Custom Libaries ##########################
###################################################################
from params import ANA_CPU_NUM, ANA_CPU_TIME, ANA_MEM, PARTITION_NAME, ENVIRONMENT_NAME
from paths import paths
from atomic_numbers import atomic_symbol_dict
from aux_functions import change_status
from aux_functions import (mark_trc_as_failed, move_all_to_failed, move_trc_to_finished, get_comp_time_string, start_job, get_trc_dict,
                           move_log, move_csm, move_com, move_xtb_out, move_rdk_xyz, move_xtb_xyz, delete_calc_dir, delete_chk_file)

###################################################################
########################### FUNCTIONS #############################
###################################################################

def execute_ana_job(ddb_no):
    """Execute analysis job for a given DDB number. Has to be run after DFT job is finished."""
    job_type = 'ana'
    ana_dir_path = paths['calc'] / f'{ddb_no}' / 'ana'
    ana_job_path = paths['calc'] / f'{ddb_no}' / 'ana' / 'job'
    change_status(ddb_no, ['ana'],['p'])

    ana_dir_path.mkdir(parents=True, exist_ok=True)

    get_py_ana(ddb_no)
    get_job_ana(ddb_no)
    start_job(ana_job_path, job_type, ddb_no)

def get_job_ana(ddb_no):
    """Create SLURM job script for analysis job."""
    job_name = f'{ddb_no}_ana'
    runtime_str = get_comp_time_string(ANA_CPU_TIME)
    mem = ANA_MEM
    job_path = paths['calc'] / f'{ddb_no}' / 'ana' / 'job'
    py_path = paths['calc'] / f'{ddb_no}' / 'ana' / f'{ddb_no}_ana.py'
    job_str = f"""#!/usr/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={PARTITION_NAME}
#SBATCH --partition=skylake-96
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={ANA_CPU_NUM}
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

def get_py_ana(ddb_no):
    """Create Python script for analysis job."""
    ana_dir_path = paths['calc'] / f'{ddb_no}' / 'ana'
    py_path = ana_dir_path / f'{ddb_no}_ana.py'

    commandStr = f"""#!/usr/bin/python
import subprocess
import sys
sys.path.append('{paths['src']}')
from pathlib import Path
from ana_functions import (extract_results)
from aux_functions import (mark_trc_as_failed, move_all_to_failed)

try:
    extract_results({ddb_no})
except:
    mark_trc_as_failed({ddb_no}, ['ana'], 17)
    move_all_to_failed({ddb_no})
"""

    with py_path.open('w') as file:
        file.write(commandStr)

def extract_results(ddb_no):
    "Main function of ana_functions.py to extract results from finished DFT calculation. Also, moves files and changes statuses."
    save_rdk_coord(ddb_no)
    save_xtb_coord(ddb_no)
    if is_log_success(ddb_no):
        save_dft_coord(ddb_no)
        save_dft_vib(ddb_no)
        if is_equi_struc(ddb_no):
            change_status(ddb_no, ['dft'], ['s'])
        elif get_trc_dict(ddb_no)['dft'] != 'c_rt' and not is_equi_struc(ddb_no):
            change_status(ddb_no, ['dft'], ['wrt'])
            change_status(ddb_no, ['ana'], ['w'])
            return
        elif get_trc_dict(ddb_no)['dft'] == 'c_rt' and not is_equi_struc(ddb_no):
            change_status(ddb_no, ['dft'], ['s_t'])
        move_log(ddb_no)
        move_csm(ddb_no)
        move_com(ddb_no)
        move_xtb_out(ddb_no)
        move_xtb_xyz(ddb_no)
        move_rdk_xyz(ddb_no)
        delete_calc_dir(ddb_no)
        change_status(ddb_no, ['ana'], ['s'])
        move_trc_to_finished(ddb_no)
        delete_chk_file(ddb_no)
    else:
        mark_trc_as_failed(ddb_no, ['dft'], 16)
        move_all_to_failed(ddb_no)


def load_mol_pkl(ddb_no):
    mol_path = paths['mol_pkl'] / f'{ddb_no}_mol.pkl'
    with mol_path.open('wb') as f:
        mol = pickle.load(f)

    return mol

def load_result_json(ddb_no):
    json_path = paths['json'] / f'{ddb_no}.json'
    with json_path.open('r') as f:
        json_dict = json.load(f)
    return json_dict

def save_result_json(result_dict, ddb_no):
    json_path = paths['json'] / f'{ddb_no}.json'
    with json_path.open('w') as f:
        json.dump(result_dict, f, indent=4)

def save_rdk_coord(ddb_no):
    rdk_xyz_path = paths['calc'] / f'{ddb_no}' / 'rdk' / f'{ddb_no}_rdk.xyz'
    json_dict = load_result_json(ddb_no)

    coord_list = []
    with rdk_xyz_path.open('r') as f:
        cart_vec = f.readlines()[2:]
    for line in cart_vec:
        parts = line.split()
        atom_type = parts[0]
        x, y, z = map(float, parts[1:])
        coord_list.append([atom_type, x, y, z])

    json_dict['rdk_coord'] = coord_list

    save_result_json(json_dict, ddb_no)

def save_xtb_coord(ddb_no):
    xtb_xyz_path = paths['calc'] / f'{ddb_no}' / 'xtb' / 'crest_best.xyz'
    json_dict = load_result_json(ddb_no)

    coord_list = []
    with xtb_xyz_path.open('r') as f:
        cart_vec = f.readlines()[2:]
    for line in cart_vec:
        parts = line.split()
        atom_type = parts[0]
        x, y, z = map(float, parts[1:])
        coord_list.append([atom_type, x, y, z])

    json_dict['xtb_coord'] = coord_list

    save_result_json(json_dict, ddb_no)

def is_log_success(ddb_no):
    log_file_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'
    with log_file_path.open('r') as f:
        content = f.read()
    return content.count('Normal termination of Gaussian 16') == 4

def save_dft_coord(ddb_no):
    json_dict = load_result_json(ddb_no)
    log_file_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'
    with log_file_path.open('r') as f:
        content = f.read()
    stat_point_found = content.find('-- Stationary point found.')
    if stat_point_found == -1:
        return []
    opt_block = re.search(
                r"Standard orientation:\s+-+\s+Center\s+Atomic\s+Atomic\s+Coordinates \(Angstroms\)\s+-+\s+((?:.+\n)+?) +-+\s",
                content[stat_point_found:]
    )

    if not opt_block:
        return []
    
    coord_list = []
    for line in opt_block.group(1).strip().splitlines():
        parts = line.split()
        if len(parts) >= 6:
            atom_symbol = atomic_symbol_dict[int(parts[1])]
            x, y, z = map(float, parts[3:6])
            coord_list.append([atom_symbol, x, y, z])

    json_dict['dft_coord'] = coord_list

    save_result_json(json_dict, ddb_no)

def save_dft_vib(ddb_no):
    json_dict = load_result_json(ddb_no)
    log_file_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'

    with log_file_path.open('r') as f:
        content = f.read()

    # Search for all frequencies and intensities
    frequencies = re.findall(r"Frequencies --\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", content)
    intensities = re.findall(r"IR Inten\s+--\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", content)

    ir_spectrum = []
    for freq_set, inten_set in zip(frequencies, intensities):
        for freq, inten in zip(freq_set, inten_set):
            ir_spectrum.append([float(freq), float(inten)])

    json_dict['dft_vib'] = ir_spectrum

    save_result_json(json_dict, ddb_no)

def is_equi_struc(ddb_no):
    json_dict = load_result_json(ddb_no)
    log_file_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'

    with log_file_path.open('r') as f:
        content = f.read()

    # Search for negative frequencies
    frequencies = re.findall(r"Frequencies --\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)", content)
    has_negative_frequency = any(float(freq) < 0 for freq_set in frequencies for freq in freq_set)

    # Set 'is_equi' and save the result
    json_dict['is_equi'] = not has_negative_frequency
    save_result_json(json_dict, ddb_no)

    return not has_negative_frequency