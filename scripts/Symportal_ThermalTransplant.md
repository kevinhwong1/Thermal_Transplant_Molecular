# ITS2 Analysis with Symportal

Project: Thermal Transplant Molecular

References:
* [SymPortal Framework](https://github.com/didillysquat/SymPortal_framework)

This analysis was performed on the URI HPC Andromeda cluster.

# Copying SymPortal into my own directory for permission access

```
cd /data/putnamlab/kevin_wong1/
rsync -rl /opt/software/SymPortal/0.3.21-foss-2019b/ ./SymPortal/
```

# Modifying settings and config file

```
module load SymPortal
mv settings_blank.py ./settings.py
mv sp_config_blank.py ./sp_config.py
base64 /dev/urandom | head -c50
```

`nano settings.py`

pasted secret key in `settings.py` file in `SECRET_KEY = 'YFW6bSD/qlLZ6rZ5rEDwsyPEqFOfC2fbRSfIAy6GcoRrnRN+N7'` and saved it.

`nano sp_config.py`

* user_name = "undefined" --> user_name = "kevinhwong1"

* user_email = "undefined" --> user_email = "kevinhwong1@gmail.com"

# Creating the SymPortal conda environment

```bash
interactive
module load Miniconda3/4.9.2
module load SymPortal

conda env create -f $EBROOTSYMPORTAL/symportal_env.yml #try doing this as a script next time?

exit
```

#### Creating the database and testing

```bash

interactive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal

python3 manage.py migrate

##### Populating reference_sequences

python3 populate_db_ref_seqs.py

module unload SymPortal

exit
```

`nano symportal_setup.sh`

```bash
#!/bin/bash
#SBATCH --job-name="SP_setup"
#SBATCH -t 500:00:00
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/SymPortal
#SBATCH --exclusive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal/0.3.21-foss-2020b
module unload

export PYTHONPATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/lib/python3.7/site-packages:$PYTHONPATH

export PATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/bin:$PATH

##### Checking installation
python3 -m tests.tests

echo "Mission Complete!" $(date)

```

```
ANALYSIS COMPLETE: DataAnalysis:
       name: testing
       UID: 2

DataSet analysis_complete_time_stamp: 20211210T153514



Cleaning up after previous data analysis test: 2
Deleting /glfs/brick01/gv0/putnamlab/kevin_wong1/SymPortal/outputs/analyses/2
Cleaning up after previous data loading test: 7
Deleting /glfs/brick01/gv0/putnamlab/kevin_wong1/SymPortal/outputs/loaded_data_sets/7
Mission Complete! Fri Dec 10 16:35:52 EST 2021
```

#### Loading Data

```bash
scp Wong_Kevin_ThermalTrans_ITS2_Meta.csv kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/20210609_Thermal_Transplant_ITS2/Wong_Kevin_ThermalTrans_ITS2_Meta.csv
```

`nano symportal_load.sh`

```bash
#!/bin/bash
#SBATCH --job-name="SP_load"
#SBATCH -t 500:00:00
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/SymPortal
#SBATCH --exclusive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal/0.3.21-foss-2020b
module unload SymPortal/0.3.21-foss-2020b

export PYTHONPATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/lib/python3.7/site-packages:$PYTHONPATH

export PATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/bin:$PATH

main.py --load /data/putnamlab/kevin_wong1/20210609_Thermal_Transplant_ITS2/raw_data \
--name Thermal_transplant_1 \
--num_proc $SLURM_CPUS_ON_NODE \
--data_sheet /data/putnamlab/kevin_wong1/20210609_Thermal_Transplant_ITS2/Wong_Kevin_ThermalTrans_ITS2_Meta.csv

```

#### Running analysis

```bash
#!/bin/bash
#SBATCH --job-name="SP_analysis"
#SBATCH -t 500:00:00
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/SymPortal
#SBATCH --exclusive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal/0.3.21-foss-2020b
module unload SymPortal/0.3.21-foss-2020b

export PYTHONPATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/lib/python3.7/site-packages:$PYTHONPATH

export PATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/bin:$PATH

# Checking dataset number
./main.py --display_data_sets

# Running analysis
./main.py --analyse 8 --name ThermalTrans_analysis --num_proc $SLURM_CPUS_ON_NODE

# Checking data analysis instances
./main.py --display_analyses

```

```bash
scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/SymPortal/outputs/analyses/3/20211216T114913/its2_type_profiles /Users/kevinwong/MyProjects/Thermal_Transplant_Molecular/output/ITS2/.
```
