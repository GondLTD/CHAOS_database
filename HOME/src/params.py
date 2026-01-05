########## GENERAL ##########
JOB_LIMIT = 201
PARTITION_NAME = 'XXX'
USER_NAME =  'XXX'
############ RDK ############
RDKIT_CPU_NUM = 1
RDKIT_RUNTIME = 2 * 24 * 3600       # 2 days in seconds
RDKIT_MEM = 10
############ XTB ############
XTB_CPU_NUM = 6
XTB_MEM = 2000
EWIN = 1
XTB_RUNTIME_LIMIT = 3 * 24 * 3600     # 3 days in seconds
############ DFT ############
DFT_CPU_NUM = 6
DFT_MEM_HIGH = 20 # GB
DFT_MEM_LOW  = 10 # GB
DFT_MEM_PER_CPU_HIGH = 4000 #MB
DFT_MEM_PER_CPU_LOW = 2000 #MB
DFT_RUNTIME_LIMIT = 3 * 24 * 3600     # 3 days in seconds
# DFT_FUNCTIONAL = 'bv5lyp IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)' # korrektes b3lyp
DFT_FUNCTIONAL = 'wB97XD'
############ ANA ############
ANA_CPU_NUM = 1
ANA_CPU_TIME = 600
ANA_MEM = 1000
############ MISC ###########
ENVIRONMENT_NAME = 'XXX'