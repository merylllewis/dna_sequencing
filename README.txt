compare_dna_seq_methods.py is a python script to compare the accuracy of two competing DNA sequencing methods.

The script takes the CSV files with data pertaining to each biochemistry as command line arguments.

The following outputs will be stored in the same directory as the input data :
1) 2 new CSV files containing the new base calls assigned to each spot in both cycles
2) 2 new txt files containing the error in assigning base calls, the contrast (confidence) values in both cycles and the number 
of spots that were missed (no signal)

base_call_assigner.py contains all the core functions that evaluate base calls, calculate error in base call assignment and write to disk

To run, in current form,

$ python compare_dna_seq_methods.py data/sequencing_data_biochem1.csv data/sequencing_data_biochem2.csv

This will produce files:
1) data/sequencing_data_biochem1_new_calls.csv
2) data/sequencing_data_biochem1_analysis_log.txt
3) data/sequencing_data_biochem2_new_calls.csv
4) data/sequencing_data_biochem2_analysis_log.txt