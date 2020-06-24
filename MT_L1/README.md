Level 1 (L1) MT time series data processing
===========================================
We present the following Python3 codes which were used for generating the [AusLAMP Musgraves Level 1 MT time series products](http://dx.doi.org/10.25914/5eaa30de53244):

  1. `01_check.py` 
  2. `02_merge.py`
  3. `03_ascii_2_bin.py`
  4. `04_rotate.py`
  5. `05_make_netCDF.py`

Note that these scripts were designed to work on outputs from **Earth Data Logger** instruments only.

### Level 0 time series ASCII data and metadata
Virtual links (`00_raw_ASCII_new/`) are used to access the Level 0 ASCII time series data and the associated metadata text files.
 
```bash 
00_raw_ASCII_new/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018 -> /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/?A/Level_0_Concatinated_Time_Series_ASCII/

00_raw_ASCII_new/WA_and_SA_metadata -> /g/data/uc0/my80_dev/sheng_data/WA_and_SA_metadata/
```

## L1 processing Pipeline

### Concatenate Level 0 daily ASCII time series files to ASCII time series files over the whole recording period 
The daily Level 0 **Earth Data Logger (EDL)** ASCII time series files for each station were concatenated over the whole recording period. This step involved:

+ Checking each folder and associated files for any potential issues including incorrect number of samples per day, missing files, evidence of instrument drift, etc.
   
+ For each station, the daily concatinations were merged into a single file encompassing the whole recording interval (minus the first and final day) for each electromagnetic component (EX, EY, BX, BY, BZ). 

+ Rotation and downsampling (from 10Hz to 1Hz) routines were performed based on information provided in the station metadata text files (e.g. `00*/WA_metadata`). A Matlab code was provided which described how to run the rotation and downsampling routines.  

### Generate NetCDF files
For each station, there are individual files for the variables `.B*`, `.E*`, `.ambientTemperature` as well as an associated metadata text file. All variables and metadata are incorporated into a single NetCDF file. That is, one `NetCDF` file per station over the whole recording period.

### Example of how to run the MT_L1 codes at the NCI
To run these codes, you will need to open up a terminal, for example, using NCI's Virtual Desktop Infrastructure. You will then need to load in the appropriate system modules - the key modules for these codes are [python3](https://www.python.org/), [openmpi](https://www.open-mpi.org/) and [mpi4py](https://mpi4py.readthedocs.io/en/stable/).

```bash
module load pbs   setuptools/23.2.1-py3.5   openmpi/1.8.4
module load python3/3.5.2     mpi4py/2.0.0-ompi-1.8.4-py3.5 
module load gcc/5.2.0
module load netcdf4-python/1.2.4-ncdf-4.3.3.1-py3.5
```
Work space folders will need to be created for each of the L1 time series processing steps:

```bash
mkdir 01_workspace 02_workspace 03_workspace 04_workspace 05_workspace
```
The following commands will generate an information file that will be used to check for any potential issues with the Level 0 ASCII time series files (e.g. incorrect number of samples per day, missing files, evidence of instrument drift):

```bash
wc -l /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/WA/Level_0_Concatinated_Time_Series_ASCII/*/*/* | sed  '$d' > 01_workspace/WA.log

wc -l /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/SA/Level_0_Concatinated_Time_Series_ASCII/*/*/* | sed  '$d' > 01_workspace/SA.log
```
Now we can run the `01_check.py` script:

```bash
python3 ./01_check.py > 01_workspace/01.log
```
Next, check the output file `01_workspace/01.log` to see if any stations have issues with incorrect number of samples or missing files. If there are any issues encountered, the identified stations will need to be rectified (or excluded) before proceeding.

Now we can run `02_merge.py` which will merge the daily L0 ASCII files into a single ASCII file over the whole recording period for each station. As a default, the first and final days of recording are excluded when running `02_merge.py`. However, you can change what days are excluded by editing the `merge` function in `02_merge.py` (i.e. changing the values for `cut_head` and `cut_end` - only whole numbers are accepted). Let's now run `02_merge.py` using mpi: 

```bash
mpirun -np 8 python3 ./02_merge.py # check the merged data and log files in 02_workspace/
```
The next step is to run `03_ascii_2_bin.py`, which makes use of `ascii2bin.cc` to converte the merged ASCII files into intermediate binary files. This step is fundamental for accelerating the subsequent I/O operations. 

```bash
### ASCII to binary

cp ascii2bin.cc 03_workspace
cd 03_workspace/ 
g++ -O2 ascii2bin.cc -o ascii2bin  # compile and generate the c language program
cd ..
wc -l 02_workspace/merged_data_WA/* | sed '$d' > ./03_workspace/WA.tmp
wc -l 02_workspace/merged_data_SA/* | sed '$d' > ./03_workspace/SA.tmp
cat 03_workspace/WA.tmp  | awk '{s = $2; gsub("02", "03", s); gsub("data","data_bin", s); printf("./03_workspace/ascii2bin %s %s %s.bin \n", $2, $1, s)}' > ./03_workspace/WA.cmd
cat 03_workspace/SA.tmp  | awk '{s = $2; gsub("02", "03", s); gsub("data","data_bin", s); printf("./03_workspace/ascii2bin %s %s %s.bin \n", $2, $1, s)}' > ./03_workspace/SA.cmd
```

To run 03_ascii_2_bin.py:

```bash
mpirun -np 6 ./03_ascii_2_bin.py
```
To check if all transformation were performed as expected, look inside the log files - there should be `None` and nothing else.

```bash
cat 03_workspace/*log  | awk '{print $NF}' | sort | uniq 
```

Following this, rotations and downsampling are performed on the binary data and the outputs and logs are stored in `04_workspace/`:

To run 04_rotate.py:

```bash
mpirun -np 6 python3 ./04_rotate.py # the output are processed time series stored as binary files in 04_workspace.
```

The final step is to generate the L1 NetCDF files with the associated metadata. The outputs and logs will be stored inside 05_workspace/ (one .nc file per station): 

```bash
mpirun -np 6 python3 ./05_make_netCDF.py 
```
