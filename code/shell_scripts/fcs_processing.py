#! /usr/bin/env python2.7
"""
This script parses and cleans up a provided Flow Cytometry Standard (fcs) file
and saves it as a Comma Separated Value (csv).
"""
import os
import re
import numpy as np
import pandas as pd
import optparse
import FlowCytometryTools as flow

# #########################################################################
def main():
    # Initialize the option parser
    parser = optparse.OptionParser()
    
    #Add options.
    parser.add_option('-i', '--input_file', dest='filename', help='name of single\
            file to be processed.', metavar="filename")
    parser.add_option('-d', '--directory', dest='inputdir', help='name of\
            input directory to be processed')
    parser.add_option('-p', '--pattern', dest='pattern', help='filename\
            pattern to parse files.')
    parser.add_option('-o', '--output', dest='out',
            help='name of output directory')
    parser.add_option('-c', '--channel', action='append', dest='channels',
            help=' individual channels to extract. Each channel must have its\
            own -c flag.') 
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose',\
            help='print progress to stdout', default=False)
    parser.add_option('-f', '--force', action='store_true', dest='force',
            help='force saving of files to output directory if needed.',
            default=False)
    
    # get the ops and args
    ops, args = parser.parse_args()
    
    # List files
    if (ops.inputdir == None) & (ops.filename == None):
        raise ValueError('no input directory/file provided! Please indicate\
                the input directory that contains the fcs files')
    
    # get all the files in the directory
    files = []
    if ops.inputdir != None:
        usr_files = np.array(os.listdir(ops.inputdir))
        # Use the pattern to identify all of the files.
        files_idx = np.array([ops.pattern in f for f in usr_files])
        file_names = usr_files[files_idx]
        
        #Add the input directory ahead of each file.
        for f in file_names:
            files.append('%s/%s' %(ops.inputdir, f)) 
    else:
        files.append(ops.filename)

    # Test that the output directory exists and is empty. 
    if ops.out != None:
        if os.path.isdir(ops.out) == False:
            os.mkdir(ops.out)
            print("Made new ouptut directory %s. I hope that's okay..." %ops.out)
        elif len(os.listdir(ops.out)) != None: 
            if ops.force == True:
                cont = 'y'
            else:
                cont = raw_input('Output directory is not empty! Continue? [y/n]: ')


    # loop through the files
    for i,f in enumerate(files):
        # consider only the fcs files
        if f.endswith('.fcs'):
            # read the file
            fcs_file = flow.FCMeasurement(ID=f, datafile=f)
            # if there are not set channels, get all the channels
            if ops.channels == None:
                fcs_data = fcs_file.data
            # if channels are indicated, collect only those
            else:
                fcs_data = fcs_file.data[ops.channels]
   
            #parse the file name to change the extension
            filename  = re.sub('.fcs', '.csv', f)
            

            #Determine if they should be saved to an output directory or not.
            if ops.out == None:
                fcs_data.to_csv(filename)
                if ops.verbose == True:
                    print(f + ' -> ' + filename)

            else:
                find_split = filename.rsplit('/', 1)
                if len(find_split) != 1:
                    filename = filename.rsplit('/', 1)[1]
                
                # Determine how to save the file.
                if len(os.listdir(ops.out)) != None:
                    if cont.lower() == 'y':
                        fcs_data.to_csv(ops.out + '/' + filename)
                    if ops.verbose == True:
                        print(f + ' -> ' + ops.out + '/' + filename)
                    else:
                        raise ValueError('output directory is not empty.')
    

if __name__ == '__main__':
    main()
    print('thank you -- come again')
