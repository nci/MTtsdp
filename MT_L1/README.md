MT Data Processing (L1)
=====

# 1. Raw- and meta-data 
 is fold is for new ASCII data (time series) and meta data (attributes).
     
Virtual links inside `00_raw_ASCII_new/` are used here to access raw- and meta-data.

1.1 
``` 
00_raw_ASCII_new/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018 -> /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/?A/Level_0_Concatinated_Time_Series_ASCII/

00_raw_ASCII_new/WA_and_SA_metadata -> /g/data/uc0/my80_dev/sheng_data/WA_and_SA_metadata//
```

# 2. Processing Pipeline

## 2.1 Raw ASCII data to level0 ASCII data
To merge raw data (time series and attributes), there are several steps to transfer them into level0 ASCII data.

+ 0) Check each fold and files for possible wrongness that includes wrong number of samples per day, missing of some files, etc.
   
+ 1) At each station, merge daily recordings into a single recording for each component, and cut out
     the data in the leading and the final day.

+ 2) Rotate, downsampling. The rotation need a few parameters that are 
     stored in `00*/WA_metadata` for each station. Those parameters are used
     to rotate, scale, etc. The downsampling is to change from 10Hz to 1Hz. 
     All these three operations are explained by a tutorial matlab script. 

*NOTE1: some stations at WA are using different recorders and hence operations above do not apply to them. These four stations should be abandoned in the first place. Those four stations are declared in Nigel's jupyter notebooks.*

*NOTE2: a few stations at WA and SA have wrong recordings. They are recognized in the 0th step.*

## 2.2 Generate netCDF files
Specifically, at each station, there are individual files for channels of `.B?`, `.E?`, `.ambientTemperature`, and metadata parameters,
here we need to include all those data into a single netCDF file. That is one `netcdf` file for one station.

# 3. How to run
Here are step by step run:

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
mkdir 01_workspace 02_workspace 03_workspace 04_workspace

###
#  Step 0 check
#
#  Check possible wrongness at each station
###
# generate information file that is used for check
wc -l /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/WA/Level_0_Concatinated_Time_Series_ASCII/*/*/* | sed  '$d' > 01_workspace/WA.log
wc -l /g/data/my80/proc_mus/States_and_Territories/SA/Long_period/SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018/SA/Level_0_Concatinated_Time_Series_ASCII/*/*/* | sed  '$d' > 01_workspace/SA.log
# run the check script
python3 ./01_check.py > 01_workspace/01.log 
# !!! Now, check the output file 01_workspace/01.log to find out wrong stations.
#     Please fix the wrong stations. Do not enter the following steps until there
#     is not any problems reported by 01_workspace/01.log.

###
#  Step 1 merge
#
#  This merge many ASCII files into a single ASCII file at each station.
#  A few stations are excluded after checking (done in step 0).
#  Please edit the days that should be dropped, for example, 
#    the first day and the last day. That can be done by editing
#    the parameter at the head of `02_merge.py`.
#    The default setting is to drop the first and the last day.
#  mpi running is used. That can be mpirun in VDI, or through
#    qTorch system in gadi. VDI comes with 8 cores, and it is 
#    powerful enough to finish the running in a few minutes that can
#    be less than the queue time in gadi.
#  NOTE: double check the version of python3, openmpi, and mpi4py.
###
mpirun -np 8 python3 ./02_merge.py # check the merged data and log files in 02_workspace/

###
#  Step 2 Rotation and downsampling
#  
#  First, we transfer merged ASCII file into binary to accelerate 
#    following IO operations. The output and logs are in 03_workspace/
#  Second, we read from binary data to do rotation and downsample.
#    The output and logs are in 04_workspace/
###
# ASCII to binary
cp ascii2bin.cc 03_workspace
cd 03_workspace/ 
g++ -O2 ascii2bin.cc -o ascii2bin  # compile and generate the c language program
cd ..
wc -l 02_workspace/merged_data_WA/* | sed '$d' > ./03_workspace/WA.tmp
wc -l 02_workspace/merged_data_SA/* | sed '$d' > ./03_workspace/SA.tmp
cat 03_workspace/WA.tmp  | awk '{s = $2; gsub("02", "03", s); gsub("data","data_bin", s); printf("./03_workspace/ascii2bin %s %s %s.bin \n", $2, $1, s)}' > ./03_workspace/WA.cmd
cat 03_workspace/SA.tmp  | awk '{s = $2; gsub("02", "03", s); gsub("data","data_bin", s); printf("./03_workspace/ascii2bin %s %s %s.bin \n", $2, $1, s)}' > ./03_workspace/SA.cmd
# run
mpirun -np 6 ./03_ascii_2_bin.py
# check if all transformation are done and correct. Inside the log files,  there should be `None` and nothing else
cat 03_workspace/*log  | awk '{print $NF}' | sort | uniq 
# rotate
mpirun -np 6 python3 ./04_rotate.py # the output are processed time series stored as binary files in 04_workspace.

###
#  Step 3 get netcdf
#
#  Output and logs are inside 05_workspace/ (one .nc file per station)
#  NOTE: We include all the metadata attributes into netcdf. Those attributes correspond to raw time series.
#        Some of those attributes are meaningless to the processed time series. That is discribed in netcdf files.
###
mpirun -np 6 python3 ./05_make_netCDF.py 

```




