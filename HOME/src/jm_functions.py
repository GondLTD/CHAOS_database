###################################################################
######################## Standard Libaries ########################
###################################################################
import pickle
import subprocess
import json

###################################################################
###################### Third-Party Libaries #######################
###################################################################
from rdkit import Chem

###################################################################
######################## Custom Libaries ##########################
###################################################################
from paths import paths
from aux_functions import (mark_trc_as_failed,move_all_to_failed,write_error_to_err, init_json_dict)
from rdk_functions import (execute_rdk_job, monitor_rdk_status)
from xtb_functions import (execute_xtb_job, monitor_xtb_status)
from dft_functions import (execute_dft_job, monitor_dft_status)
from ana_functions import (execute_ana_job)

###################################################################
########################### FUNCTIONS #############################
###################################################################


def read_num_jobs(user_name):
    """
    Returns the number of the running SLURM jobs for the given user 'user_name'.
    Args:
    user_name - user name within SLURM system, must be string
    """
    try:
        result = subprocess.run(['squeue', '-u', user_name], stdout=subprocess.PIPE, text=True)
        output = result.stdout
        job_count = len(output.strip().split('\n')) - 1

        return job_count
    
    except Exception as e:
        return None
    
def process_jobs(jobs_to_run):
    """
    This function is the central unit of Abinitiomaker. It delegates the different types of jobs to be dealt with.
    It iterates through trc_pending an searches for .trc file, which are pending for calculation.
    Some filters are applied: 1. SMILES code is available, 2. a mol object can be made from the given smiles code.
    Then, a json file for all results of a molecule is created (if missing). Then, the job is specified via reading the
    trc file (get_job_type). Finally, the specified job is executed.
    """
    with open(paths['key_ddb_to_smiles'], 'rb') as f:
        smiles_dict = pickle.load(f)

    count = 0  
    for file_path in paths['trc_pending'].iterdir():
        if count >= jobs_to_run:
            break
        
        ddb_no = get_ddb_from_path(file_path)
        try:
            smiles = smiles_dict[ddb_no]
        except KeyError:
            mark_trc_as_failed(ddb_no, ['rdk'], 12)
            move_all_to_failed(ddb_no)
            continue
        try:
            test_mol = Chem.MolFromSmiles(smiles)
        except:
            mark_trc_as_failed(ddb_no, ['rdk'], 18)
            move_all_to_failed(ddb_no)
            continue
        try:
            if not (paths['json'] / f'{ddb_no}.json').is_file():
                init_json_dict(ddb_no)
        except:
            mark_trc_as_failed(ddb_no, ['rdk'], 15)
            move_all_to_failed(ddb_no)
            continue

        try:
            job_type = get_job_type(file_path)
        except:
            write_error_to_err(ddb_no, 13)
            move_all_to_failed(ddb_no)
            continue

        if job_type == 'idle':
            monitor_rdk_status(ddb_no)
            monitor_xtb_status(ddb_no)
            monitor_dft_status(ddb_no)

        if job_type not in ['idle', 'fail']:
            execute_job(ddb_no, smiles, job_type)
            count += 1
        else:
            if job_type == 'fail':
                mark_trc_as_failed(ddb_no, ['rdk', 'xtb', 'dft', 'ana'], 2)
                move_all_to_failed(ddb_no)

def get_ddb_from_path(file_path):
    """
    Returns the DDB number of a given trc path'.
    Args:
    file_path - from the trc_pending directory
    """
    return int(file_path.stem)

def get_job_type(file_path):
    """
    Specifies the job type. Returns one of the following:
        rdk    - A rdkit job will be executed
        xtb    - A crest conformer search will be executed
        xtb_rt - A retry of the conformer search will be executed
        dft    - A quantum mechanical density functional theory calculation will be performed, including optimization, vibrational analysis, NMR calculation and cosmo calculation
        dft_rt - A retry of dft will be executed, due to a imaginary frequency
        ana    - Analysis will be performed on the results of the previously performed dft jobs
        idle   - Everything what can be done is currently done regarding the specific component
        fail   - this component failed somewhere, now a clean up job is performed.
    Args:
        file_path - from the trc_pending directory
    """
    with file_path.open('r', encoding='utf-8') as trc_file:
                trc_dict = json.load(trc_file)

    rdk_status = trc_dict.get('rdk')
    xtb_status = trc_dict.get('xtb')
    dft_status = trc_dict.get('dft')
    ana_status = trc_dict.get('ana')

    status_list = [rdk_status, xtb_status, dft_status, ana_status]

    if all(status == 'w' for status in status_list):
        return 'rdk'
    if rdk_status == 's' and all(status == 'w' for status in status_list[1:]):
        return 'xtb'
    if rdk_status == 's' and xtb_status == 'wrt' and all(status == 'w' for status in status_list[2:]):
        return 'xtb_rt'
    if rdk_status == 's' and xtb_status == 's' and all(status == 'w' for status in status_list[2:]):
        return 'dft'
    if rdk_status == 's' and xtb_status == 's' and dft_status == 'wrt' and ana_status == 'w':
        return 'dft_rt'
    if rdk_status == 's' and xtb_status == 's' and dft_status == 'c' and ana_status == 'w':
        return 'ana'
    if rdk_status == 's' and xtb_status == 's' and dft_status == 'c_rt' and ana_status == 'w':
        return 'ana'
    if any(status in ['p', 'r', 'r_rt'] for status in status_list) and 'f' not in status_list:
        return 'idle'
    if 'f' in status_list:
        return 'fail'
    
    return 'fail'

def execute_job(ddb_no, smiles, job_type):
    """
    Executes the job for the given job_type. The jobs are defined in seperate file (rdk_functions, xtb_functions, dft_functions, ana_functions)
    """
    if job_type == 'rdk':
        execute_rdk_job(ddb_no, smiles)
    elif job_type == 'xtb':
        execute_xtb_job(ddb_no, smiles)
    elif job_type == 'xtb_rt':
        execute_xtb_job(ddb_no, smiles, is_retry=True)
    elif job_type == 'dft':
        execute_dft_job(ddb_no, smiles)
    elif job_type =='dft_rt':
        execute_dft_job(ddb_no, smiles, is_retry=True)
    elif job_type == 'ana':
        execute_ana_job(ddb_no)
