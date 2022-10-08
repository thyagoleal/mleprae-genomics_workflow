# Workflow for bacterial genomics

### About

This workflow takes in unprocessed fastq.gz files from short-read sequencers then pre-process, map, and call variants against the *M. leprae* genome TN reference (pre-configured). The workflow was built upon python2 code from Chloe Loiseau and adapted with python3 code and Snakemake workflow manager by me. The workflow is delivered through a docker image  and run using singularity locally with all dependencies automatically.

The workflow was designed to be user-friendly to those without coding skills. However, a basic understanding is still necessary to run things without problems and to fix any issues that may happen.

If you have any questions, please write-me up. 

#### Pre-requisites

* A computer with any linux distribution or MacOS (not tested in Windows, but should work);
* Snakemake and singularity installed;
* FASTQ.gz files sample-demultiplexed, separated by lanes or not, single-ended (SE) or paired-end (PE). Please, make sure that your files follow the [official naming scheme](https://raw.githubusercontent.com/thyagoleal/mleprae-genomics_workflow/main/docs/illumina.txt?token=GHSAT0AAAAAABVB3FP7ZRMT7ZQOCPQIDX5GYX6RXEA) from Illumina;
* Internet connection for downloading software (at least for the first time).

### Instructions

#### Setting up

1. Download this repository by using ` git clone github.com/thyagoleal/mleprae-genomics_workflow/` or through the Code -> Download ZIP green button.
2. Unzip the contents and navigate inside the directory;

![image-20220715165418518](.imgs/image-20220715165418518.png)

2. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) following the instructions in the software page;

3. Install [Mamba](https://github.com/mamba-org/mamba) with: `conda install mamba -n base -c conda-forge`;

4. Launch a terminal making sure you're inside the *mleprae_genomics_workflow* directory. Then,  create the environment with the specified dependencies using mamba:

   ```shell
   conda activate base
   mamba env create --name mleprae_genomics_ENV --file environment.yml
   ```

5. If everything installed correctly run `conda activate mleprae_genomics_ENV`;
6. Now you can use the workflow and its installed tools. It is likely that  all of these steps have been done for you already;

#### Configuring and running the workflow

Now that you have a conda environment with the installed tools, you can configure your samples.

1. Create a directory named `data` if you don't have one. It should be inside the  `mleprae-genomics_workflow` directory.  :warning: Careful about not overwriting an already populated `data`dir. 

2. Inside of `data` create two another directories named `samples` and `refs`;

3. Make sure you are in the root workflow directory `mleprae-genomics_workflow`. Then, to create the sample files and feed the workflow, execute the code below:

   ``` shell
   scripts/prepare_input.py /path/to/samples/dir/
   ```

   This script will check the names of your fastq files, merge lanes (default), fix suffixes (i.e., standardize `fastq.gz` instead of `fq.gz` etc), and check which samples are paired-end or not. Then, it will build the `samples.tsv` and `units.tsv` files necessary to start the workflow with soft links inside `data/samples` pointing to the original files. 

   > The original files will remain **untouched** in the data/samples dir. The workflow-ready files will be in data/samples/prepared-fastqs/ by default (you can change it, check the prepare_input.py script --help). 

   

