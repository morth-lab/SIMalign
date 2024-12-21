############### READ ME!!! #################

###################################
### Download the Git repository ###
###################################
to download the SIMAlign commandline program, you need to clone it from Github.
If you have GitHub Desktop installed, navigate to "Repository" in the toolbar of the window, click "clone a repository", go to the URL banner and paste the following WRL into the first empty bar:
https://github.com/TheGamingEngineer/SIMAlign_commandline.git

in the lower bar, enter the destination path, where you want the repository to be cloned to and then click "clone". 

If you want to download the repository using a command prompt, first ensure that you have git installed. 
then navigate to the destination folder, where you want the folder to be and enter the following code: 

$ git clone https://github.com/TheGamingEngineer/SIMAlign_commandline.git

###########################
### Running the Program ###
###########################
SIMAlign comamndline is designed to work within a virtual environment with the necessary installments for the application to work. 

In order to initially access the virtual environment, open a command prompt and navigate to the folder:
SIMALign_commandline

To activate this environment, use
    $ conda activate simalign

To deactivate an active environment, use
     $ conda deactivate

OBS!: certain functionality is only available with proper internet connection. 
Hence it is impossible to run any other homology method than user-specified in offline mode. 


to run SIMAlign as a commandline write:
	$ python SIMALign_parser_v2 -- QUERY path/to/query_file --TEMPLATE path/to/template_file_1 path/to/template_file_2 [...] --JOB_KEY job_key --HOMOLOGY_SEARCH_METHOD method --MAX_DISTANCE int --MAX_INITIAL_RMSD int --afdb50 True --afdb_swissprot True --afdb_proteome True --pdb100 True --FOLDSEEK_MODE "tmalign" --THRESHOLD float --NUMB_HOMO int --SEQUENCE_IDENTITY float --RESULT_DIR path/to/result_directory


int = an integer
float = a float number

The Following of the variables has default values: 
<variable>               == <default value>
--TEMPLATE               == None 
--JOB_KEY                ==  program-generated job-key 
--HOMOLOGY_SEARCH_METHOD == foldseek 
--MAX_DISTANCE           == 6 
--MAX_INITIAL_RMSD       == 5
--afdb50                 == True 
--afdb_swissprot         == False 
--afdb_proteome          == False 
--pdb100                 == False 
--FOLDSEEK_MODE          == tmalign 
--THRESHOLD              == 0.7 
--NUMB_HOMO              == 0 
--SEQUENCE_IDENTITY      == 0.6 
--RESULT_DIR             == User_data/{JOB_KEY}

the only obligatory option is --QUERY. 
When defining route to file, be sure to include the whole path from root. 

#################################
### Data Structure of Results ###
#################################

Then the program is done, the user should have a result directory of the following structure: 

RESULT_DIR:
|
|- SIMAlign_{JOB_KEY}.zip
|  |- tar.gz
|     |- m8-file
|  |- alignment.aln
|  |- m8-file
|  |- hotspots_mode_1.html
|  |- hotspots_mode_2.html
|  |- pymol_output.pse
|  |- AF-<x>-F1-model-v4.pdb
|  |-  <-||-> *x



