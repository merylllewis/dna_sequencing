#!/usr/bin/env python3

""" dna_sequencing.py 

Analyzes dye signal intensities of DNA sequenced through a particular biochemistry and assigns a basecall to each spot 
of the sequence based on the base that best matches the dye with the maximum intensity at that spot. 

CSV files with data pertaining to each biochemisty must be passed in as command line arguments

"""

__author__ = " Meryl Lewis "
__email__ = " merylllewis@gmail.com "


import basecall_reader as br
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

    # Find best dye to base map and assign correct basecalls to all spots in both cycles
    basecalls_1 = br.find_min_error_basecalls ( dye_intensities_1, ref_1 )
    basecalls_2 = br.find_min_error_basecalls ( dye_intensities_2, ref_2 )

    # Find dye contrast for both cycles
    avg_contrast_1 = br.get_dye_contrast( np.array( dye_intensities_1 ) )
    avg_contrast_2 = br.get_dye_contrast( np.array( dye_intensities_2 ) )

    print( avg_contrast_1 )
    print( avg_contrast_2 )

    # This is where the new csv file will be saved on disk
    new_filename = filename [ 0:-4 ] + "_new_calls.csv"

    # Write new basecalls to CSV
    br.write_to_csv( data_read, basecalls_1, basecalls_2, new_filename )
