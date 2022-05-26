# WGGC Single Cell

## IMPORTANT!!!
The pipeline and this README.md is designed to work on the HPC of the HHU, Germany.
General instructions for simply running with snakemake exist, but are very short and not completely tested yet.

Update 15.05.2022: There has been a small change in the config.yaml format. For more information, checkout configfiles/config_example.yaml.

Update 26.05.2022: Fix the issue with the DGE rule that came with the change from the 15.05.2022.

## Setup directory tree:
### the scRNAseq_Pipeline folder containing the pipeline
|  
|- path/to/scRNAseq_Pipeline  
&nbsp;&nbsp;&nbsp;&nbsp;|- cluster  
&nbsp;&nbsp;&nbsp;&nbsp;|- configfiles  
&nbsp;&nbsp;&nbsp;&nbsp;|- data  
&nbsp;&nbsp;&nbsp;&nbsp;|- dependencies  
&nbsp;&nbsp;&nbsp;&nbsp;|- envs  
&nbsp;&nbsp;&nbsp;&nbsp;|- errorMessage  
&nbsp;&nbsp;&nbsp;&nbsp;|- scripts  
&nbsp;&nbsp;&nbsp;&nbsp;|- snakemakeScripts  
&nbsp;&nbsp;&nbsp;&nbsp;|- workDirectory  
  
### an example of an output folder
|  
|- path2/two/pipeline_output_of_project  
&nbsp;&nbsp;&nbsp;&nbsp;|- clusterLogs  
&nbsp;&nbsp;&nbsp;&nbsp;|- csv  
&nbsp;&nbsp;&nbsp;&nbsp;|- logs  
&nbsp;&nbsp;&nbsp;&nbsp;|- plot  
&nbsp;&nbsp;&nbsp;&nbsp;|- outputs  
&nbsp;&nbsp;&nbsp;&nbsp;|- shinyApp  
&nbsp;&nbsp;&nbsp;&nbsp;|- workDirectory  

Note: Some of the folders may not exists if you download the repository. The pipeline should automatically create them.  
- cluster: contains files needed for using the Pipeline on the HPC-Cluster
- configfiles: config.yamls of the files used for testing and the config_example.yaml
- data: contains the folders with your CellRanger outputs
- dependencies: folder in which the downloaded ShinyCell-master.zip and DoubletFinder-master.zip are to be stored
- envs: containing the envs.yaml needed to use snakemake with conda-env via "--use-conda" option
- errorMessage: containing the error messages
- scripts: contains the R scripts used in the pipeline run
- snakemakeScripts: contains the scripts used in the Snakefile
- workDirectory: contains a txt file with instructions on how to run the shinyApp which is copied to the shinyApp/ folder containing the outputs

- clusterLogs: directory containing logs created by using the pipeline in cluster mode
- csv: contains the pipelines .csv outputs
- logs: contains logs from snakemake
- plots: contains the plot outputs of the pipeline
- outputs: .rds files that the pipeline creates
- shinyApp: contains a program making it possible to interact with the data over a graphical interface
- workDirectory: where temporary outputs and R residues are saved

## Setup config.yml for data:
To adapt the pipeline to your data you just need to adjust the config.yaml file.
Check "config_example.yaml" in the configfiles/ directory.
It is an empty version of the config.yaml and explains how the config.yaml is to be filled out.
Check config354.yaml for an example that includes the comment descriptions.

## Setup conda environment via Snakemake (only if not working on HHU HPC):
In order to run the pipeline you need to install "Anaconda" or "Miniconda" and create a Snakemake environment on it (Pipeline currently works best with Snakemake/5.10.0).
Download the following two packages as .zip from Github and move them to "scRNAseq_Pipeline/dependencies/":
  - https://github.com/chris-mcginnis-ucsf/DoubletFinder
  - https://github.com/SGDDNB/ShinyCell 
  - https://github.com/genomics-kl/seurathelpeR
Then you need to switch into the "scRNAseq_Pipeline" directory via your terminal.  
Next you need to use the command "snakemake -p --cores 1 --use-conda" so that Snakemake can create the environments with the two downloaded packages itself.  
After the environments are installed you can use "snakemake -p --cores X --use-conda" to run the pipeline properly where X is the number of cores you would like to use.
On how to use the the pipeline see the last points of the HPC version.

## Use on the Hilbert HPC cluster of the HHU
- First download/clone the repository and put it on you HPC folder (a /scratch_gs/your_hpc_username/folder_name/ folder is recommended a single run can need over 50GB ram of space for the outputs).
- Create a config.yaml (see section above) [It can also be named any other name as long as it is a config.yaml and contains all parameters]
- **IMPORTANT**: Make sure the projectDirectoryPath is a path you have access to on the HPC. In the example config.yaml in configfiles/ the projectDirectoryPath is a path only I can access.
- Move your Cellranger output data into the data/ folder.
- Download the following two packages as .zip and move them into the dependencies/ folder:
  - https://github.com/chris-mcginnis-ucsf/DoubletFinder
  - https://github.com/SGDDNB/ShinyCell 
  - https://github.com/genomics-kl/seurathelpeR
  - https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html (ver.1.30.0 as tar.gz)
  - https://bioconductor.org/packages/release/data/annotation/html/GenomeInfoDbData.html (ver.1.2.7 as tar.gz)
  - https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html (ver.1.46.1 as tar.gz)
  - https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html (ver.1.24.0 as tar.gz)
  - https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html (ver.1.16.0 as tar.gz)
- the versions are important for compatibility with the environment for the HPC
- It is possible that you have to change "account" in cluster/cluster.json to a HPC-project name of yours
- Log into the HPC via terminal
- create a screen via "screen -S screen_name" command, this way you won't lose any progress even if you lose internet connection. You can reconnect to your screen via "screen -r screen_name".
- Use the "ssh snakemake-node" command to switch to the snakemake node.
- On the snakemake-node switch to the folder containing the repository.
- Use "bash clusterExecution.sh path/to/config.yaml" to start the pipeline. If this is the first time you start the pipeline, snakemake will need to create the conda environments for the rules. This can take a while (around 20 minutes to a day on the HPC). Then it will install the two packaged in dependencies in the environments. Afterwards you need to repeat the "bash.clusterExecution.sh path/to/config.yaml" command to continue the pipeline.
- The pipeline may stop and throw errors after a while. This is completely normal. If you are missing an essential input in your config.yaml "numberOfCells" or it is time to to look at a plot to fill out parts of the config.yaml you could not fill out before a certain part of the pipeline finished running, the pipeline will stop and throw an error to get your attention. To check if this is a "normal" error check in the clusterLogs/ the .errors file in which the error occurred. If the Error under the RuleException is "MissingInputError" then it is normal and the error message will tell you what you need to do next. If it is another kind of error please contact me or raise an issue since this would mean there is an issue with the pipeline.
- After filling out a new input field in the config.yaml continue running the pipeline with the "bash clusterExecution.sh path/to/config.yaml".
- After the pipeline is completely done running, use the exit command twice to first exit the snakemake-node and then your screen.
  (Use screen -list on the login-node to see the screens you have if you want to check if you forgot a screen).
- You can run the pipeline with two different dataset inputs at the same time. Just create two screens with different names and run on the respective screens "bash clusterExecution.sh path/to/config1.yaml" and "path2/two/config2.yaml".
- **Important**: The pipeline approximates the resources needed to request walltime and RAM from the HPC meaning there are times it won't ask for enough resources. In the snakemakeScripts/ directory is now a Python script called "PipelineOptions.py" where you can add additional RAM or walltime. You could also simply increase "numberOfCells" in the config.yaml, but if you increase it too drastically (like to 400000 cells) then there will be an HPC error because more resources were requested than existing on the HPC. If the screen with says the job is still running even though the HPC says the job is not running then not enough walltime was requested. **For some reason the HPC no longer throws an error when time runs out before the job is finished.** If near the end of the .errors file is a line with "/bin/bash: line 1:  ____ Killed" or something similar with "killed" then not enough RAM was requested. If this happens snakemake might lock the project folder. In this case use the resumePipeline.sh instead of the clusterExecution.sh to continue running the pipeline. You should be able to switch back to clusterExecution.sh after using resumePipeline.sh once unless the same error here happens again. resumePipeline.sh is used the same way as clusterExecution.sh.
- **Important**: Sometimes there is an issue with the .sh files when moving from a windows PC to the HPC. In this case use the command "sed -i 's/\r//g' filename.sh" to fix the 'command \r not found' error with bash.

## To-do list:
- documentation
- DGE (reparation almost done, HTO does not perform DGE, might need to fix that)
- new and improved benchmarking
- Seurat parallelization with future.
