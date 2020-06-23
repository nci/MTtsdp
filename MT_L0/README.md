## Packed Raw Data and Level 0 MT time series data processing levels 

We present the following Python3 codes:

   1. **01_pack.py** - Produces a zip file of the raw instrument data for each site in the survey (e.g. station1.zip, station2.zip, ...).
     
   2. **02_time_series.py** - Generates Level 0 concatenated EX, EY, BX, BY, BZ ASCII files at a per station per day granularity. These ASCII files have associated metadata available. 
   
   3. **03_netcdf.py** - Creates a single Level 0 concatenated netCDF file per station per day for variables EX, EY, BX, BY, BZ. The MT time series metadata is also available in the netCDF file.

## Raw data and metadata 
     
Virtual links are used here for accessing the __Earth Data Logger__ raw MT time series data and associated metadata.
 
```
00_data_virtual_link/01_data_location -> /g/data/my80/proc_mus/Musgraves/

00_data_virtual_link/02_data_location -> /g/data/my80/proc_mus/Musgraves/

00_data_virtual_link/03_data_location_SA -> /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/SA/Level_0_Concatinated_Time_Series_ASCII

00_data_virtual_link/03_data_location_WA -> /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/WA/Level_0_Concatinated_Time_Series_ASCII

00_data_virtual_link/03_metadata_dir -> /g/data/uc0/nre900/MT/MT_jupyter_notebook/WA_and_SA/
```

## Example of how to run the MT_L0 codes at the NCI

To run these codes, you will need to open up a terminal, for example, using NCI's Virtual Desktop Infrastructure (VDI). You will then need to load in the appropriate system modules - the key modules for these codes are [python3](https://www.python.org/), [openmpi](https://www.open-mpi.org/) and [mpi4py](https://mpi4py.readthedocs.io/en/stable/).  

```bash
module load pbs   setuptools/23.2.1-py3.5   openmpi/1.8.4
module load python3/3.5.2     mpi4py/2.0.0-ompi-1.8.4-py3.5 
module load gcc/5.2.0
module load netcdf4-python/1.2.4-ncdf-4.3.3.1-py3.5

```
Next you will need to create a work space folder for each of the three time series processing codes:

```bash
mkdir 01_workspace 02_workspace 03_workspace
```
###  Running the codes in parallel

For this example we use 7 cores, but more cores could be requested for further acceleration.

#### 1. Packed raw time series data
```bash
mpirun -np 7 python3 ./01_pack.py # output and logs are in 01_workspace
```
#### 2. Level 0 ASCII time series
```bash
mpirun -np 7 python3 ./02_time_series.py # output and logs are in 02_workspace
```
#### 3. Level 0 NetCDF time series with associated metadata
```bash
mpirun -np 7 python3 ./03_netcdf.py # output and logs are in 03_workspace
```




