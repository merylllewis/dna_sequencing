#!/usr/bin/env python3

""" compare_dna_seq_methods.py 

Analyzes dye signal intensities of DNA sequenced through a particular biochemistry and assigns a base call to each spot 
of the sequence based on the base that best matches the dye with the maximum intensity at that spot. 

CSV files with data pertaining to each biochemistry must be passed in as command line arguments

"""

__author__ = " Meryl Lewis "
__email__ = " merylllewis@gmail.com "


import basecall_assigner as ba
import numpy as np
import pandas as pd
import sys


for i in range(1, len(sys.argv)):

    filename = sys.argv[i]

    print("\nReading file " + filename + "\n")

    # Read csv data into Pandas dataframe
    data_read = pd.read_csv ( filename, dtype = 'str' )

    # Extract dye intensities for both cycles from data_read
    dye_intensities_1 = ( data_read.ix[ :, 1:5 ] ).astype(float)
    dye_intensities_2 = ( data_read.ix [ :, 7:11 ] ).astype(float)

    # Extract reference sequences for both cycles from data_read
    ref_1 = np.array(data_read.ix[ :, 0 ])
    ref_2 = np.array ( data_read.ix [ :, 6 ] )

    # Find best dye to base map and assign correct base calls to all spots in both cycles. Also find dye contrast
    [ base_calls_1, dye_map_1, error_1 ] = ba.find_min_error_base_calls ( dye_intensities_1, ref_1 )
    [ contrast_1, no_signal_1 ] = ba.compute_dye_contrast( np.array( dye_intensities_1 ) )

    [ base_calls_2, dye_map_2, error_2 ] = ba.find_min_error_base_calls ( dye_intensities_2, ref_2 )
    [ contrast_2, no_signal_2 ] = ba.compute_dye_contrast( np.array( dye_intensities_2 ) )

    # This is where the new csv file and analysis log file will be saved on disk
    new_calls_csv_filename = filename [ 0:-4 ] + "_new_calls.csv"
    analysis_log_filename = filename [ 0:-4 ] + "_analysis_log.txt"

    # Write new base calls to CSV
    ba.write_new_base_calls_to_csv( data_read, base_calls_1, base_calls_2, new_calls_csv_filename )

    # Write analysis log to disk
    errors = [ error_1, error_2 ]
    dye_maps = [ dye_map_1, dye_map_2 ]
    contrasts = [ contrast_1, contrast_2 ]
    no_signals = [ no_signal_1, no_signal_2 ]

    ba.create_analysis_log_file( errors, dye_maps, contrasts, no_signals, analysis_log_filename )