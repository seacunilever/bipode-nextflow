# Bifrost Nextflow workflow

Workflow for Bifrost allowing local and Azure cloud batch execution. The main motivation is to allow flexible specification of Azure batch compute nodes to reduce costs through use of Low Priority instances which cost 20% of Dedicated instances. To do this the workflow must be able to resume any jobs which fail due to Low Priority nodes being pre-empted. The large number of short running independent computations performed by Bifrost lend themselves well to to this approach.

A secondary motivation is workflow portability across different execution environments (eg. local, HPC, cloud providers). Currently two profiles exist, local and Azure cloud.

This workflow takes in pre-formatted json data files for one cell line/chemical (referred to as one dataset).

## Todo
- Look at running from Apache Airflow with DAG triggering an Azure Container Instance (ACI) to run nextflow. This means terminal session does not need to be kept alive.

## Prerequisites

- Linux (required for Nextflow, can be WSL2 https://learn.microsoft.com/en-us/windows/wsl/install)
- Docker (note: Azure ML VMs have docker pre-installed)
- Nextflow. Azure Cloud requires version 21.04.0 or later (https://www.nextflow.io/docs/latest/getstarted.html)

## Setup

### Install Nextflow

See instructions here: https://www.nextflow.io/docs/latest/getstarted.html

For a Azure ML compute instance, this would look like
```bash
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
sudo mv ./nextflow /usr/local/bin
```
This is included in file `nextflow-install.sh`

### Containers

For local runs, build the container from the dockerfile in `./docker` tagged as `bifrost:2023`
```bash
cd ./docker
docker build -t bifrost:2023 .
```

For Azure cloud use, it is likely the container already exists in the Azure Container Repository (ACR), so just check details in `nextflow.config` are correct. If not, you must build the image locally as above and then push the built image to ACR with tag `bifrost:2023`. These instructions will show you how https://learn.microsoft.com/en-us/azure/container-registry/container-registry-get-started-docker-cli?tabs=azure-cli. For more details see https://www.nextflow.io/docs/latest/azure.html#private-container-registry.

Eg.
```bash
#login to azure cli
az login --use-device-code
#login to azure container registry (ACR)
az acr login --name bnlwees01d57280sb01acr.azurecr.io
#create an alias of the image with path to registry
docker tag bifrost:2023 bnlwees01d57280sb01acr.azurecr.io/bifrost:2023
#push image
docker push bnlwees01d57280sb01acr.azurecr.io/bifrost:2023
```

If you use different container names or details, amend `nextflow.config` accordingly

### Azure cloud config

This readme assumes these Azure resources already exist: batch account, blob storage container (associated with batch account), Azure Container Repository (ACR).

Set your Azure account details for storage, container registry and batch in a file `credentials.json` with format
```json
{
    "storageAccountName": "",
    "storageAccountKey": "",
    "batchAccountName": "",
    "batchAccountKey": "",
    "AcrUserName": "",
    "AcrPassword": ""
}
```
These account details can be found on the Azure portal http://portal.azure.com and are not included in this repo for security reasons.

A template is included in this repo as `credentials_example.json`

## Run the workflow

> **Note**
> If running on a compute cluster infrastructure, `nextflow` must be able to communicate
> with the workload manager at all times, otherwise tasks will be cancelled.
> The best way to do this is to run `nextflow` using a `screen` or `tmux`
> terminal.
>
> E.g. Screen
>
> ```bash
> # Open a named screen terminal session
> screen -S my_nextflow_run
> # run nextflow
> nextflow run -c <config> -profile <profile> <nextflow_script>
> # "Detach" screen terminal
> <ctrl + a> <ctrl + d>
> # list screen sessions
> screen -ls
> # "Attach" screen session
> screen -r my_nextflow_run
> ``

`tmux` comes pre-installed on Azure ML compute instances

### Locally

Set the run parameters in the appropriate parameter settings file `params_local.yml`
```bash
nextflow run bifrost.nf -profile local -params-file params_local.yml
```

Note that you can also run locally on an Azure VM, doing so will run nextflow on the node itself and can be useful for debugging since you do not have to wait for nodes to be provisioned. Be sure to check `n_cores` in `params_local.yml` matches the size of node you are using. Note that the default path `~/cloudfiles/code/` is a fileshare, you will want to run nextflow from the local filesystem of the VM, `~/localfiles` for performance as the working directory will be local.

### Azure cloud

Set the run parameters in the appropriate parameter settings file `params_az.yml` 
Make sure `n_cores` matches with the cloud machine set in `nextflow.config` 
Note that Azure cloud files should be prefixed with `az://`. You can also files local to where you are running the workflow, they do not need to be on azure cloud storage.

```bash
nextflow run bifrost.nf -profile az -params-file params_az.yml
```

Nextflow creates logs in the directory you run from as `.nextflow.log`, `.nextflow.log.1` etc. for each run. Filenames starting with `.` are hidden in Linux, you can see them with `ls -a`. To unhide them you can just change the name eg.

`mv .nextflow.log nextflow.log`

#### Resuming

Add `-resume` when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

#### Cleanup

Nextflow keeps all logs and files generated for runs in the work directory unless they are removed, so the workflow can be resumed. A convenience script `nextflow-clean.sh` is included for removing files for all runs. Note that this is for Azure cloud execution only, not local runs.

## Known issues

- Running locally using the Azure `az` profile is a problem, as download of azure plugin, like so many things, gets blocked by Zscalar
```
Downloading plugin nf-azure@0.14.2
Plugin with id nf-azure not found in any repository
```

- Conda installed Nextflow causes `java.lang.UnsupportedOperationException: Not a valid Azure Blob Storage file attribute view: interface`. At least 22.10.1 did anyway.

- Do not use `scratch = true` directive in any nextflow process as this causes issues writing files to cloud blob store as working directory

- If using azure cloud, run from an Azure VM. Using `SEACserver` has incurred timeout errors like
```
Error executing process > 'CONC_RESPONSE_ANALYSIS (2010)'

Caused by:
  timeout
```
`.nextflow.log`
```
java.lang.RuntimeException: java.net.SocketTimeoutException: timeout
	at rx.exceptions.Exceptions.propagate(Exceptions.java:57)
```

## Azure cloud notes

A pool is created on workflow execution and deleted after finishing. Low priority nodes are used to save cost. The compute pool scaling formula will scale up and down the pool size based on the amount of work in the queue up to a provided maximum number of nodes.  To change the maximum, amend `max_nodes` in `params_az.yml`, note that at peak times this number of low priority nodes may not be available.  Processes are set to retry 3 times before failing so if any tasks are kicked by Azure they will restart. If all tasks do not complete re-running with `-resume` will process these as long as to work directory has not been touched or settings changed.

Costs for compute can be found here https://azure.microsoft.com/en-gb/pricing/details/batch/windows-virtual-machines/. Azure machine type `Standard_F4s_v2` is used by default and set in `nextflow.config`.

The batch jobs nextflow submits can be viewed in the Azure portal which also shows compute pool status and load
While jobs are running and the pool nodes are active, the files including logs can be viewed which helps job introspection
The output blob store can be interacted with from the Portal or Azure Storage Explorer (GUI application)

![alt text](blog-nextflow-and-azure-batch-part-1-of-2-1.png "Azure")

## Developer notes

Be sure to check `azure.storage.tokenDuration` is set long enough for the duration of your run

Many execution platforms are inefficient if a workflow tries to execute many short running processes. It can take more time to schedule and request resources for each small instance than bundling the short processes into a larger process task. Nextflow channel operators collate groups together inputs into batches which can run for longer with the short tasks themselves parallelised inside the process script which is what is done in `conc_response_analysis.py`.

Azure machine type `Standard_F4s_v2` is used by default but can be changed to any machine type available to the Azure batch account by amending `vmType` in `nextflow.config`. Make sure to set `n_cores` in `params_az.yml` to use all the cores on the machine.

The tens of thousands of probe pkl files from `split_data.py` were taking almost 2hrs to be transferred to the cloud blob storage work directory due to the large number of files
The solution implemented here compresses the pkl files into a `.tar.gz`, then each of the concentration response modelling processes then extract the files they need from it.
This approach could have data transfer implications if huge numbers of parallel jobs are kicked off in the cloud, as all the `.tar.gz` of all probe files are part of each job, rather than only the ones required to be processed. A similar approach has been applied to the files for the fits to the process that combines them with `compress_output.py` with Nextflow task ids being added to the `tar.gz` files to avoid name clashes.

In order to facilitate inspection of dataset results as early as possible during running, `groupTuple()` is used with a `groupKey()` which is produced from the number of fit tar.gz files or each dataset. Doing this means Nextflow knows when all the of fits are available to run the compress process for the dataset rather than waiting until all the fits for all the datasets are complete. This is covered as a tip in Nextflow docs on `groupTuple()` https://www.nextflow.io/docs/latest/operator.html#grouptuple. The fits channel is sorted to have results available as soon as possible allowing inspection and an early decision to be made to terminate the run if problems are found with the fits rather than waiting until all datasets complete.

The docker command is given the `--cpus 1` argument by Nextflow if `cpus` is not explicitly set in `executor.cpus`

## Useful links

Nextflow docs for Azure https://www.nextflow.io/docs/latest/azure.html

Setup for Azure
https://shaunchuah.github.io/posts/setting-up-azure-with-nextflow
https://seqera.io/blog/nextflow-and-azure-batch-part-1-of-2/
https://techcommunity.microsoft.com/t5/healthcare-and-life-sciences/covid-variant-analysis-on-azure-using-nextflow-part-2/ba-p/3075741

Batch pool scaling https://learn.microsoft.com/en-us/azure/batch/batch-automatic-scaling
Pool scaling formula was taken from here https://github.com/Azure/doAzureParallel/blob/master/R/autoscale.R

Nextflow slack channel for help
https://www.nextflow.io/slack-invite.html