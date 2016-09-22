# shell scripts
Here lie bits of shell and python code which are useful for saving, renaming,
and storing data generated during this project. Below are short descriptions of
each function and how they should be used. **For all python scripts,
information of usage can be found by typing `./script.py --help`.**

## `deposit_data`
**language** : shell

**purpose** : Transfer a single data file or group of files to the appropriate
 	folder on our group data server. Note that you must have already
	set up an rsa key with Griffin. 

**arguments** :

	$1 name of source file.
	$2 username on dataserver.
	$3 destination within `mwc_induction` on data server.

**example** :

	>$ sh deposit_data 20160714_wt_O1_RBS1_0uMIPTG.csv <username> flow/csv/

## `fetch_data`
**language** : shell

**purpose** : Retreive a single datafile or a group of file from your account
	or another users account on the group data server. NOte that you must have
	already set up an rsa key with Griffin.

**arguments** : 

	$1 Your username.
	$2 Username of the owner of desired data files.
	$3 Name of or pattern within desired file(s).
	$4 Destination on your local machine for transferred files.

**example** :

	>$ sh fetch_data <username> <other_username> flow/fcs/201607*_wt	~/Desktop/

## `file_rename.py`
**language** : python

**purpose** : Easily rename a large group of files generated using MACSQuant
	flow cytometer. 

**arguments** :

	-d, --dir : DIRECTORY
		Name of directory conntaining files to be renamed.

	-t, --template : CSV FILE
		Path to csv containing new names of files. The structure of
		this csv file should be the same as the layout of samples in
		the 96-well plates. Please stick to the standardized naming
		convention of `YYYYMMDD_strain_operator_rbs_xuMIPTG``

	-e, --ext : PATTERN
		Extension to append to each file. If not specified, '.fcs' will
		be used.

	-o, --output : DIRECTORY
		Path to output directory. If none is provided, files will be
		renamed in place. If this directory does not exist, it will be
		automatically made. If it exists but is not empty, the user
		will be prompted for confirmation that the renaming should
		proceed.

	-v, --verbose : NONE
		Optional flag to print progress of renaming to screen.

	-f, --force : NONE
		Optional flag to force renaming of files is ouput directory is
		not empty.

**example** :
	
	>$ ./file_rename.py -d ./20160714_flow_data/ -t ./20160714_template.csv -o ./201601714_renamed/ -vf
		
## `fcs_processing.py`
**language** : python 2.7

**purpose** : Read a provided Flow Cytometry Standard (fcs) file or directory
	of fcs files, extract the desired channels, and save them as a Comma Separated
	Value (csv) file in a specified output folder.

**arguments** : 

	-i, --input_file :  FILENAME
		Name of individual file to process.

	-d, --directory : DIRECTORY
		Name of directory containing files to be processed.	

	-p, --pattern : PATTERN
		Pattern within desired filenames to be processed.
	
	-o, --output : DIRECTORY
		Path to output directory. If this directory does not exist, it
		will be made. If this directory exists but is not empty, user
		input will be required to continue with execution.

	-c, --channel : CHANNEL NAME
		Name of channel (exact) to be extracted from the fcs file. If
		multiple channels are desired, multiple calls of -c must be
		called.

	-v, --verbose : NONE
		Optional flag to print progress of processing to screen.

	-f, --force : NONE
		Optional flag to force renaming of files if output directory is
		not empty.

**example** : 
	
	>$ ./fcs_processing.py -d 20160714_flow_data/ -p 20160714_wt* -o ./ -c FITC-A -vf



