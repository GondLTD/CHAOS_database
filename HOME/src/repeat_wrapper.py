import time
import subprocess
from pathlib import Path
from params import ENVIRONMENT_NAME
from paths import paths

commands = f'module purge; module load anaconda3/latest; conda activate {ENVIRONMENT_NAME}; python3 main.py'

while True:
    # repeat the main job processing script every 60 seconds, arbitrary time interval
    subprocess.run(commands, shell=True, cwd=paths['src'])
    time.sleep(60)