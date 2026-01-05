from pathlib import Path

home_path    = XXX
scratch_path = XXX

paths = {'home'             : home_path,
         'src'              : home_path / 'src',
         'scratch'          : scratch_path,
         'calc'             : scratch_path / 'calc',
         'json'             : home_path    / 'json',
         'err'              : home_path / 'err',
         'trc_pending'      : home_path / 'trc' / 'trc_pending',
         'trc_failed'       : home_path / 'trc' / 'trc_failed',
         'trc_finished'     : home_path / 'trc' / 'trc_finished',
         'key_ddb_to_smiles': home_path / 'data'/ 'KEY_ddb_to_smiles.pkl',
         'log'              : home_path / 'log',
         'csm'              : home_path / 'csm',
         'xtb_out'          : home_path / 'xtb_out',
         'xtb_xyz'          : home_path / 'xtb_xyz',
         'rdk_xyz'          : home_path / 'rdk_xyz',
         'com'              : home_path / 'com',
         'failed_log'       : home_path / 'failed' / 'log',
         'failed_csm'       : home_path / 'failed' / 'csm',
         'failed_xtb_out'   : home_path / 'failed' / 'xtb_out',
         'failed_xtb_xyz'   : home_path / 'failed' / 'xtb_xyz',
         'failed_rdk_xyz'   : home_path / 'failed' / 'rdk_xyz',
         'failed_com'       : home_path / 'failed' / 'com',
         'failed_chk'       : home_path / 'failed' / 'chk',
         'chk'              : scratch_path.parent / 'chk'}