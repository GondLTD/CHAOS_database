###################################################################
######################## Standard Libaries ########################
###################################################################
import shutil
import json
import subprocess
from datetime import datetime

###################################################################
###################### Third-Party Libaries #######################
from rdkit import Chem

###################################################################
######################## Custom Libaries ##########################
###################################################################
from paths import paths
from error_dict import error_dict

###################################################################
########################### FUNCTIONS #############################
###################################################################

def mark_trc_as_failed(ddb_no, failed_status_list, error_code):
    """
    Manipulates the trc file to 'f'.
    Args:
        ddb_no: DDB number
        failed_status_list: a list of the states which failed, which are to be chosen from ['rdk', 'xtb', 'dft', 'ana']
        error_code: The number of the error, can be found in error_dict
        """
    trc_file_path = paths['trc_pending'] / f'{ddb_no}.trc'

    with trc_file_path.open('r', encoding='utf-8') as trc_file:
        trc_dict = json.load(trc_file)
    for err_idx, failed_status in enumerate(failed_status_list):
        trc_dict[failed_status] = 'f'
    write_error_to_err(ddb_no, error_code)
    with trc_file_path.open('w', encoding='utf-8') as trc_file:
        json.dump(trc_dict, trc_file, indent=4)

def move_all_to_failed(ddb_no):
    """Moves the failed trc fail to trc_failed. Therefore, it is no longer considered during process_job().
    Moreover, all results found so far are transferred to the failed directory for further insights on errors."""
    src_path = paths['trc_pending'] / f'{ddb_no}.trc'
    dst_path = paths['trc_failed'] / f'{ddb_no}.trc'
    shutil.move(str(src_path), str(dst_path))
    ### possible paths to transfer
    rdk_path        = paths['calc'] / f'{ddb_no}' / 'rdk' / f'{ddb_no}_rdk.xyz'
    xtb_out_path    = paths['calc'] / f'{ddb_no}' / 'xtb' / f'{ddb_no}_xtb.out'
    xtb_xyz_path    = paths['calc'] / f'{ddb_no}' / 'xtb' / 'crest_best.xyz'
    com_path        = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.com'
    csm_path        = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_dft.cosmo'
    log_path        = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'
    chk_path        = paths['chk']  / f'{ddb_no}_dft.chk'

    if rdk_path.exists() and rdk_path.is_file():
        move_rdk_xyz(ddb_no, fail_mode=True)
    if xtb_out_path.exists() and xtb_out_path.is_file():
        move_xtb_out(ddb_no, fail_mode=True)
    if xtb_xyz_path.exists() and xtb_xyz_path.is_file():
        move_xtb_xyz(ddb_no, fail_mode=True)
    if com_path.exists() and com_path.is_file():
        move_com(ddb_no, fail_mode=True)
    if csm_path.exists() and csm_path.is_file():
        move_csm(ddb_no, fail_mode=True)
    if log_path.exists() and log_path.is_file():
        move_log(ddb_no, fail_mode=True)
    if chk_path.exists() and chk_path.is_file():
        move_chk_file(ddb_no)

    add_all_errors(ddb_no)
    delete_calc_dir(ddb_no)

def add_all_errors(ddb_no):
    """Collects all .err files from the calculation directory and appends them to the central .err file."""
    ddb_path = paths['calc'] / f'{ddb_no}'
    err_file_path = paths['err'] / f'{ddb_no}.err'
    
    with err_file_path.open('a') as err_file:
        for subdir_path in ddb_path.iterdir():
            if subdir_path.is_dir(): 
                error_files = subdir_path.glob('slurm-*.err')
                
                err_file.write(f'####{subdir_path.name.upper()}_ERROR####\n')
                
                for error_file in error_files:
                    with error_file.open('r') as file:
                        content = file.read()
                        err_file.write(content + '\n')


def write_error_to_err(ddb_no, error_code):
    """Writes the error from error_dict to a given .err file with a timestamp"""
    err_file_path = paths['err'] / f'{ddb_no}.err'
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # YYYY-MM-DD HH:MM:SS
    error_message = f"{current_time} - {error_dict[error_code]}"
    with err_file_path.open('a', encoding='utf-8') as file:
        file.write(error_message)


def get_comp_time_string(comp_time):
        seconds = comp_time
        days, seconds = divmod(seconds, 86400)  # 86400 seconds per day
        hours, seconds = divmod(seconds, 3600)  # 3600 seconds per hour
        minutes, seconds = divmod(seconds, 60)  # 60 seconds per minute

        # Formatierung der Zeit als D-hh:mm:ss, wobei D den ersten Tag als 0 zählt
        return f"{days}-{hours:02}:{minutes:02}:{seconds:02}"

def start_job(job_path, job_type, ddb_no):
    subprocess.run(['sbatch', str(job_path)], cwd=job_path.parents[0])
    change_status(ddb_no, [job_type], ['r'])

def change_status(ddb_no, status_list, state_list):
    trc_file_path = paths['trc_pending'] / f'{ddb_no}.trc'
    with trc_file_path.open('r', encoding='utf-8') as trc_file:
        trc_dict = json.load(trc_file)
    for status_idx, status in enumerate(status_list):
        trc_dict[status] = state_list[status_idx]
    with trc_file_path.open('w', encoding='utf-8') as trc_file:
        json.dump(trc_dict, trc_file, indent=4)

def init_json_dict(ddb_no):
    json_path = paths['json'] / f'{ddb_no}.json'
    init_dict = dict()
    with json_path.open('w') as f:
        json.dump(init_dict, f, indent=4)

def move_trc_to_finished(ddb_no):
    src_path = paths['trc_pending'] / f'{ddb_no}.trc'
    dst_path = paths['trc_finished'] / f'{ddb_no}.trc'
    shutil.move(str(src_path), str(dst_path))

def get_atom_count_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms()

def get_trc_dict(ddb_no):
    trc_file_path = paths['trc_pending'] / f'{ddb_no}.trc'
    with trc_file_path.open('r', encoding='utf-8') as trc_file:
        trc_dict = json.load(trc_file)
    
    return trc_dict

def move_xtb_out(ddb_no, fail_mode=False):
    src_path = paths['calc'] / f'{ddb_no}' / 'xtb' / f'{ddb_no}_xtb.out'
    if fail_mode:
        dst_path = paths['failed_xtb_out'] / f'{ddb_no}_xtb.out'
    else:
        dst_path = paths['xtb_out'] / f'{ddb_no}_xtb.out'
    shutil.move(str(src_path), str(dst_path))

def move_rdk_xyz(ddb_no, fail_mode=False):
    src_path = paths['calc'] / f'{ddb_no}' / 'rdk' / f'{ddb_no}_rdk.xyz'
    if fail_mode:
        dst_path = paths['failed_rdk_xyz'] / f'{ddb_no}_rdk.xyz'
    else:
        dst_path = paths['rdk_xyz'] / f'{ddb_no}_rdk.xyz'
    shutil.move(str(src_path), str(dst_path))

def move_xtb_xyz(ddb_no, fail_mode=False):
    src_path = paths['calc'] / f'{ddb_no}' / 'xtb' / 'crest_best.xyz'
    if fail_mode:
        dst_path = paths['failed_xtb_xyz'] / f'{ddb_no}_xtb.xyz'
    else:
        dst_path = paths['xtb_xyz'] / f'{ddb_no}_xtb.xyz'
    shutil.move(str(src_path), str(dst_path))

def move_log(ddb_no, fail_mode=False):
    src_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.log'
    if fail_mode:
        dst_path = paths['failed_log'] / f'{ddb_no}.log'
    else:
        dst_path = paths['log'] / f'{ddb_no}.log'
    shutil.move(str(src_path), str(dst_path))

def move_csm(ddb_no, fail_mode=False):
    src_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}_dft.cosmo'
    if fail_mode:
        dst_path = paths['failed_csm'] / f'{ddb_no}_dft.cosmo'
    else:
        dst_path = paths['csm'] / f'{ddb_no}_dft.cosmo'
    shutil.move(str(src_path), str(dst_path))

def move_com(ddb_no, fail_mode=False):
    src_path = paths['calc'] / f'{ddb_no}' / 'dft' / f'{ddb_no}.com'
    if fail_mode:
        dst_path = paths['failed_com'] / f'{ddb_no}.com'
    else:
        dst_path = paths['com'] / f'{ddb_no}.com'
    shutil.move(str(src_path), str(dst_path))

def move_chk_file(ddb_no):
    src_path = paths['chk'] / f'{ddb_no}_dft.chk'
    dst_path = paths['failed_chk'] / f'{ddb_no}_dft.chk'
    shutil.move(str(src_path), str(dst_path))

def delete_calc_dir(ddb_no):
    dir_path = paths['calc'] / f'{ddb_no}'
    if dir_path.exists() and dir_path.is_dir():
        shutil.rmtree(dir_path)

def delete_chk_file(ddb_no):
    chk_file_path = paths['chk'] / f'{ddb_no}_dft.chk'
    if chk_file_path.exists() and chk_file_path.is_file():
        chk_file_path.unlink()