MT Data Processing (L0) 
=====

# 1. Raw- and meta-data 
     
Virtual links are used here for accessing raw- and meta-data.
 
```
00_data_virtual_link/01_data_location -> /g/data/my80/proc_mus/Musgraves/

00_data_virtual_link/02_data_location -> /g/data/my80/proc_mus/Musgraves/

00_data_virtual_link/03_data_location_SA -> /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/SA/Level_0_Concatinated_Time_Series_ASCII

00_data_virtual_link/03_data_location_WA -> /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/WA/Level_0_Concatinated_Time_Series_ASCII

00_data_virtual_link/03_metadata_dir -> /g/data/uc0/nre900/MT/MT_jupyter_notebook/WA_and_SA/
```

# 2. How to run
Here are step by step runs:

```bash
###
#  This is in the sh terminal, for example, a terminal in VDI at NCI.
###

###
#  Load modules
#
#  The key modules is python3, openmpi and mpi4py. 
#  Be careful for the verion.
###
module load pbs   setuptools/23.2.1-py3.5   openmpi/1.8.4
module load python3/3.5.2     mpi4py/2.0.0-ompi-1.8.4-py3.5 
module load gcc/5.2.0
module load netcdf4-python/1.2.4-ncdf-4.3.3.1-py3.5

###
#  Make working space
###
mkdir 01_workspace 02_workspace 03_workspace

###
#  Parallel Runing
#
#  7 cores are used below. More cores (e.g., 16 or 32) can be used for further accelerations.
#  However, too many cores (e.g., 64 or much more) will not introduce substantial accelerations.
###
### 1. packing
mpirun -np 7 python3 ./01_pack.py # output and logs are in 01_workspace
### 2. for time series
mpirun -np 7 python3 ./02_time_series.py # output and logs are in 02_workspace
### 3. for netcdf
mpirun -np 7 python3 ./03_netcdf.py # output and logs are in 03_workspace

```




