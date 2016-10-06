#! /usr/bin/env python
"""
This script renames a directory of Flow Cytometry Standard (fcs) files in
accordance to a user-provided template of file names.
"""
import os
import optparse
import shutil
import numpy as np
import pandas as pd

# ########################################################################
def main():
    #Initialize the option parser
    parser = optparse.OptionParser()

    #Add the various options.
    parser.add_option('-d', '--dir', dest='targ', help='directory\
        containing files to be renamed', metavar='DIRECTORY')
    parser.add_option('-t', '--temp', dest='temp', help='path to csv\
        template for renamed files', metavar='FILE')
    parser.add_option('-e', '--ext', dest='ext', help='target file\
        extension', default='.fcs', metavar='EXTENSION')
    parser.add_option('-p', '--pattern', dest='pattern', help='target file\
        pattern for parsing', metavar='PATTERN')
    parser.add_option('-o', '--output', dest='out', help='path to output\
    directory', metavar='DIRECTORY')
    parser.add_option('-v', '--verbose', action='store_true',\
            dest='verbose', default=False, help='print progress to stdout')
    parser.add_option('-f', '--force',  action='store_true', dest='force',
            default=False, help='force creation of any necessary output\
            directories')
    #Get the options and arguments.
    ops, args = parser.parse_args()

    #Ensure that the required directory and template files are provided.
    if (ops.targ == None) | (ops.temp == None):
        raise ValueError('directory and template file are required.')

    #Load the renaming template file.
    fname_temp = pd.read_csv(ops.temp, header=None)
    
    #Read the old files that contain the file_extension pattern.
    old_files = np.array(sorted(os.listdir(ops.targ)))

    if ops.pattern != None:
    #Identify all of the target files which contain the proper file extension
    #and the pattern to parse the files
        ext_bool = np.array([(ops.ext in f) and (ops.pattern in f) for f in old_files])
    else:
        ext_bool = np.array([ops.ext in f for f in old_files])

    #Consider only the files which have the correct file extension.
    targets = old_files[np.array(ext_bool)]

    #Flatten the list of new file names in the template file.
    new_names = fname_temp.values.flatten(order='F')

    #Make sure the list of new names is the same length as the targets.
    if len(new_names) != len(targets):
        raise ValueError('length mismatch between template and target files.\
                Did you forget a sample?')
    #Iterate through all files and rename.
    for i, f in enumerate(targets):
        #The the original file name and define the new name.
        orig = ops.targ + f
        renamed = new_names[i] + ops.ext 

        #Determine if it will be renamed or copied.
        if ops.out != None:
            if os.path.isdir(ops.out) == False:
                os.mkdir(ops.out)
                print("Made new output directory %s. Hope that's okay..."
                        %ops.out)
                shutil.copy(orig, ops.out + '/' + renamed)

            elif len(os.listdir(ops.out))!=0:
                if ops.force==True:
                    shutil.copy(orig, ops.out + renamed)
                else:
                    cont = input('Output directory is not empty! Continue? [y/n] :')
                    if cont.lower() == 'y':
                        shutil.copy(orig, ops.out + '/' + renamed)
                    else: 
                        raise ValueError('output directory is not empty.')
        else:
            os.rename(orig, ops.targ + renamed)
        if ops.verbose == True:
            print(f + ' -> ' + renamed)


if __name__ == '__main__':
    main()
    print('thank you -- come again')

