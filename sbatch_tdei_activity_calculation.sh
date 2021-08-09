source activate tf-gpu #anaconda tensorflow-gpu activation
export CSV=name_of_electron_configuration_file
export DIM=topological_distance_information
export ACTIVITY=target_endpoint
export SMI=columne_name_for_smiles_code_in_inputfile
export CPU=the_number_of_available_cpu
export DIR=directory_for_input_files
python calculate_tdei_activity_train_val_test_multiprocessing.py

exit 0
