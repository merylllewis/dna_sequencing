import numpy as np
import pandas as pd


def assign_basecalls( best_base_dye_map, max_intensity_dyes ):
    """
    For a given dye to base map, assigns basecalls for all spots
    
    :param best_base_dye_map: base to dye map that results in the least error
    :param max_intensity_dyes: dyes that produce the highest intensity signal at each spot in the DNA sequence
    
    :return: assigned basecalls to each spot of the DNA sequence
    
    """

    basecalls = []
                                 
    for i in range( 0, max_intensity_dyes.shape[ 0 ] ):

        # max_intensity_dye = N implies all dye intensities were 0 for this spot
        if max_intensity_dyes[ i ] != 'N' :
            basecalls.append( best_base_dye_map[ max_intensity_dyes[ i ] ] )
        else:
            basecalls.append( 'N' )
            
    return basecalls


def find_best_dye_base_map( max_intensity_dyes, dye_intensities, ref, initial_base_dye_map ):
    """
    Find the basecall to dye map that results in the least error
    
    :param max_intensity_dyes: dyes that produce the highest intensity signal at each spot in the DNA sequence
    :param dye_intensities: array containing all four dye intensities for each spot of the DNA sequence
    :param ref: reference DNA sequence
    :param initial_base_dye_map: initial guess the base to dye map
    
    :return: base to dye map that produces the least error and the corresponding error
    
    """
    # This data structure helps to map integers to bases
    test_base_combination = { 0:'A',  1:'C', 2:'G', 3:'T' }

    best_base_dye_map = initial_base_dye_map.copy( )
    min_error = max_intensity_dyes.shape[ 0 ] + 1 # initialize to max possible value + 1

    dyes = list( initial_base_dye_map.keys( ) )
    num_spots = max_intensity_dyes.shape[ 0 ]

    # Iterate through all possible dye to base maps.
    for i in range( 0, 4 ):
        for j in range( 0, 4 ):
        
            if j == i :
                # Two dyes cannot map to the same base
                continue
        
            for k in range( 0, 4 ):
            
                if k == i or k == j :
                    # Two dyes cannot map to the same base
                    continue
            
                for l in range( 0, 4 ):
                
                    if l == i or l == j or l == k :
                        # Two dyes cannot map to the same base
                        continue

                    # Test one possible dye to base combination ( based on i, j, k, l )
                    current_test_map = {
                        dyes[ 0 ] : test_base_combination[ i ],
                        dyes[ 1 ] : test_base_combination[ j ],
                        dyes[ 2 ] : test_base_combination[ k ],
                        dyes[ 3 ] : test_base_combination[ l ]
                    }
                                
                    error = 0

                    # Iterate through all spots to find the total error for this dye to base combination
                    for spot in range( 0, num_spots ):

                        # Check dye intensities for zeros. If all dye intensities are 0, mark spot with 'N'
                        if np.sum( dye_intensities[ spot, : ] ) == 0 :
                            max_intensity_dyes[ spot ] = 'N'
                            error += 1  # Since no data was recorded for this spot, this is definitely an error
                        elif current_test_map[ max_intensity_dyes[ spot ] ] != str( ref.item( spot ) ) :
                            # With the given dye to base map, the dye with max intensity does not match with reference base
                            error += 1
                    
                    if error < min_error :
                        best_base_dye_map = current_test_map.copy()
                        min_error = error

    min_error = 100 * min_error / max_intensity_dyes.shape[ 0 ]  # Reassigns the percentage error to min_ error
    return best_base_dye_map, min_error


def find_min_error_basecalls( dye_intensities, ref ):
    """
    For a given set of dye intensities and a reference DNA sequence, finds the dye to base mapping that results in 
    basecalls with the least error, as compared to the reference sequence
    
    :param dye_intensities: dataframe holding float values of intensities for four dyes corresponding to A, C, G, T bases
    :param ref: reference DNA sequence
    :return: basecalls for each spot in the DNA sequence based on the dye to base map that minimizes error
    """

    # Find the dye with the highest intensity for each spot
    max_intensity_dyes = dye_intensities.idxmax( axis = 1 )

    # Get column names (dyes) from data
    column_names = list( dye_intensities.columns.values )

    # Assign an initial guess of dye to base map. This will be modified when the correct mapping has been found
    dye_base_map = {column_names[ 0 ]:'A', column_names[ 1 ]:'C', column_names[ 2 ]:'G', column_names[ 3 ]:'T'}

    # Find the best possible dye to base map and find the error in assignment with this map
    [best_dye_base_map, percent_assignment_error] = find_best_dye_base_map( max_intensity_dyes, np.array( dye_intensities ), ref, dye_base_map )

    # Display best dye to base map and corresponding error
    print( "Percentage Error in mapping = " + str ( percent_assignment_error ) + "%"  )
    print( "Reassigned dye to base map = " + str ( best_dye_base_map ) )

    # Assign new basecalls to all spots based on new dye to base map
    basecalls = assign_basecalls(best_dye_base_map, max_intensity_dyes)

    return basecalls, best_dye_base_map, percent_assignment_error


def compute_dye_contrast( dye_intensities ):
    """
    Returns root mean square (RMS) contrast values for all spots. For each spot, contrast is computed as a ratio of the 
    intensity of the dye with maximum intensity to the sum of all four dye intensities for that spot.
    
    :param dye_intensities: intensities of all dyes for each spot in a cycle 
    :return: contrast for dye with maximum intensity for all spots in a cycle
    
    """

    num_spots_no_signal = 0 # Number of spots for which no signal was recorded
    sum_rms_contrast = 0    # Sum of individual RMS contrast values for all spots

    # Square of dye intensities needed for RMS contrast measurement
    sq_dye_intensities = np.square( dye_intensities )

    # Iterate through all spots to find the RMS contrast for each spot
    for spot in range( 0, sq_dye_intensities.shape[0 ] ):

        # This is the dye whose corresponding base was called
        max_dye_intensity = np.max( sq_dye_intensities[ spot, : ], axis = 0 )

        if max_dye_intensity != 0 :
            sum_sq_contrast = 0 # Sum of square contrasts for all dyes at this spot

            for dye in range( 0, sq_dye_intensities.shape[ 1 ] ):
                # Square contrast for a specific dye at a specific spot
                sq_contrast = max_dye_intensity - sq_dye_intensities[ spot, dye ]
                sum_sq_contrast += sq_contrast

            # RMS contrast is the square root of the average square contrast at this spot
            rms_contrast = np.sqrt( sum_sq_contrast / sq_dye_intensities.shape[ 1 ] )
            sum_rms_contrast += rms_contrast
        else:
            # If max_dye_intesity = 0, no signal was received for this spot and this dye must be excluded
            num_spots_no_signal += 1

    avg_rms_contrast = sum_rms_contrast / (sq_dye_intensities.shape[0 ] - num_spots_no_signal)

    print( "Dye contrast = ", round( avg_rms_contrast, 3 ) )
    print( "Number of spots for which no signal was received = ", num_spots_no_signal )

    return avg_rms_contrast, num_spots_no_signal


def write_new_basecalls_to_csv ( data_read, basecalls_1, basecalls_2, filename ):
    """
    Creates a new CSV file with the reassigned basecalls for each spot of each cycle

    :param data_read: data read from input csv file
    :param basecalls_1: new basecalls assigned to first cycle
    :param basecalls_2: new basecalls assigned to second cycle
    :param filename: filename to which new CSV file containing basecalls and dye intensities must be written

    """

    # Create individual dataframes for data from each cycle
    cycle_1 = pd.DataFrame ( data_read.ix [ :, 0:6 ] )
    cycle_2 = pd.DataFrame ( data_read.ix [ :, 6:11 ] )

    # Create dataframes for newly assigned basecalls for each cycle
    assigned_basecalls_1 = pd.DataFrame ( basecalls_1, columns = [ 'new_calls_1' ] )
    assigned_basecalls_2 = pd.DataFrame ( basecalls_2, columns = [ 'new_calls_2' ] )

    # Concatenate dataframes
    data_write = pd.concat ( [ cycle_1, assigned_basecalls_1, cycle_2, assigned_basecalls_2 ], axis = 1 )

    # Write new basecalls data to disk
    print ( "\nWriting new basecalls data to " + filename )
    data_write.to_csv ( filename )


def create_analysis_log_file ( errors, dye_maps, contrasts, no_signals, filename ):
    """
    Writes an analysis log of errors, dye to base maps and contrasts to disk
    
    :param errors: percentage error in basecalls for all cycles
    :param dye_maps: dye to base map which minimizes basecall error for all cycles
    :param contrasts: contrast of dye with maximum intensity for each spot for all cycles
    :param filename: filename of the analysis log
     
    """
    # Writing analysis of error, dye to base maps and intensity contrasts to disk
    print( "\nWriting analysis log to " + filename + "\n\n" )

    analysis_log = open( filename, 'w' )

    for i in range( 0, len( errors ) ):
        analysis_log.write( "For cycle " + str( i + 1 ) + ":\n")
        analysis_log.write( "Dye to base map that minimizes error : " + str( dye_maps[ i ] ) + "\n" )
        analysis_log.write( "Number of spots for which no signal was received: " + str( no_signals [ i ] ) + "\n" )
        analysis_log.write( "Error in basecalls (including spots where no signal was received): " + str( errors[ i ] ) + "%\n" )
        analysis_log.write( "Dye signal contrast with the given biochemistry: " + str( round( contrasts[ i ], 3) ) + "\n\n" )

    analysis_log.close()


