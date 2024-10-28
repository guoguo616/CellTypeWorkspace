from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

DATABASES = {
    'ENGINE': 'django.db.backends.postgresql',
    'NAME': 'cell-type-workspace',
    'USER': 'cell-type-workspace',
    'PASSWORD': 'cell-type-workspace',
    'HOST': 'localhost',
    'PORT': '5432',
}

APIURL = 'https://scdbapi.deepomics.org/'
FILEAPI = 'https://scdbapi.deepomics.org/tasks/getoutputfile/'

USER_TASK_PATH = BASE_DIR / 'workspace/user_data/'
WORKSPACE_MODULE = BASE_DIR / 'workspace/module/'
# SCQUERY_SCRIPT = BASE_DIR / 'workspace/module/sc_query_old/run.sh'

# ABSUSERTASKPATH = '/home/platform/phage_db/phage_api/workspace/user_task'
# PHAGEFASTA = '/home/platform/phage_db/phage_data/data/phage_sequence/phage_fasta/'
# PHAGEGBK = '/home/platform/phage_db/phage_data/data/phage_sequence/phage_gbk/individual_data/'
# PHAGEGFF = '/home/platform/phage_db/phage_data/data/phage_sequence/phage_gff3/individual_data/'
# PROTEINSEQUENCE = '/home/platform/phage_db/phage_data/data/phage_sequence/proteins/'
# TEMPPATH = '/home/platform/phage_db/phage_data/data/tmp/'


# CLUSTERTREEPATH = '/home/platform/phage_db/phage_data/data/phage_sequence/cluster_tree_v2/'
# CLUSTERALIGNMENTPATH='/home/platform/phage_db/phage_data/data/analysis_data/alignment/result'
# CLUSTERSEQUENCEPATH = '/home/platform/phage_db/phage_data/data/phage_sequence/group/'
# METADATA = '/home/platform/phage_db/phage_api/workspace/csv_data/'
# FASTAPATH = '/home/platform/phage_db/phage_data/data/'
# # Analysis script path
# ANALYSIS = '/home/platform/phage_db/phage_api/workspace/analysis_script/'
# TASKLOG = '/home/platform/phage_db/phage_api/workspace/task_log/'
# DEMOFILE = '/home/platform/phage_db/phage_api/demo_file/'
