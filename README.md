# TDEi_tensor
Codes for calculation of topological distance based electron interaction (TDEi) tensor.

*Required python libraries (version)
1) RDKit (2018.09.01)
2) h5py (2.9.0)
3) Pandas (0.25.0)
4) Numpy (1.15.4)

*Dependencies
Tensorflow 2 GPU version in Anaconda environment (tensorflow 2.2.0 was used in development).

How to use
1) edit 'sbatch_tdei_activity_calcultion.sh' file to give specifications for calculation.
2) specify each variables
    !!first line is to activate tensorflow-gpu anaconda envirnoment. Please modify this line accordingly!!
    'CSV': electron configuration file in ec_vector directory.
    'DIM': topological distance to be considered in TDEi tensor calculation
    'ACTIVITY': column name of target endpoint in the data file.
    'SMI': column name for smiles of the compounds in the data file.
    'CPU': the number of cpu cores for multiprocessing.
    'DIR': a directory for data (example is csv files in data directory)
3) Run the sh file (expected operating system to run the code is Linux. Code was run on Centos 7 and 8)
