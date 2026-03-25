# README for Unilever internal use
This is an extension to the general docs since the Bipode Nextflow workflow is used both internally in SERS and is also released externally.

This document covers:

__How to run on Unilever infrastructure__
- Running on a local machine
- Running from local machine with Azure Batch compute
- Running using Azure Batch compute orchestrated using nf-runner (https://dev.azure.com/SEAC-Projects/_git/nf-runner)

__Making an external release__
1. SERS Github (https://github.com/seacunilever/bipode-nextflow-httr)
1. PyPI (https://pypi.org/project/bipode-httr)
1. Bioconda (https://anaconda.org/bioconda/bipode-httr)
1. Public container on Seqera containers (https://seqera.io/containers/)
1. SERS Github (https://github.com/seacunilever/bipode-nextflow)

## Running internally from a local machine
For typical use you would use a released version with nf-runner (see later below) but for development work you will need to run the workflow from a local machine.

You can run either entirely locally with profile 'local' for local compute or with Azure Batch with profile 'az'.

Nextflow run commands in the examples in the general docs will pull the Nextflow workflow files from GitHub, so the commands need to be changed to point to the internal workflow files. It should be possible to point them to a DevOps repo https://www.nextflow.io/docs/stable/git.html#azure-repos but I have not got it to work. Often you will be developing with the files locally anyway, so change your run command to point to the local workflow path eg. `nextflow run .` for the current directory.

By default, the container (if docker is used) will be pulled from the public container on Seqera Containers. If changes need to made to it for internal use it needs to come from the Azure Container Registry (ACR) (see "Updating 'bipode-httr' container in the Azure Container Registry (ACR)" later)

#### Local machine
1. Check that `params.workDir` in `nextflow.config` is set somewhere sensible on the blob storage

2. Check that the docker image with container tag in `conf/modules.config` and profile 'local' in `nextflow.config` exists locally (pull if required)

3. Use profile 'local' in Nextflow run command like
```bash
nextflow run . -profile local etc.
```

change nextflow cmd to clone and .

#### Azure Batch compute
1. Set the following env vars which are picked up in `nextflow.config`
##### Service Principal
```bash
CLIENT_ID
CLIENT_SECRET
TENANT_ID
```

##### Storage account name
```bash
STORAGE_ACCOUNT
```

##### Container registry (ACR)
```bash
CONTAINER_REGISTRY_USER_NAME
CONTAINER_REGISTRY_PASSWORD
CONTAINER_REGISTRY_SERVER
```

##### Batch account
```bash
BATCH_ACCOUNT_NAME
BATCH_ACCOUNT_LOCATION
SUBSCRIPTION_ID
VNET
SUBNET
SUBNET_RESOURCE_GROUP_NAME
```

2. Check that `params.workDir` in `nextflow.config` is set somewhere sensible on the blob storage

3. Check that the the docker image with container tag in `conf/modules.config` and profile 'az' in `nextflow.config` exists in the Azure Container Registry (ACR)

4. Use profile 'az' in Nextflow run command
```bash
nextflow run . -profile az etc.
```

### Updating 'bipode-httr' container image in the Azure Container Registry (ACR)
The workflow is set up to pull the container from Seqera which is built from the Bioconda release of 'bipode-httr' (see [./docs/containers.md](./docs/containers.md)) for external use.

__For internal running, the container needs to come from SERS ACR__. This allows for internal development of this code

To update the container in the ACR (change version numbers and ACR registry name as required)
1. Clone the files from branch you want from 'bipode-httr' repo into `./docker` eg.
```bash
git clone --branch <branch-name> --single-branch https://SEAC-Projects@dev.azure.com/SEAC-Projects/BIFROST/_git/bipode-httr
```
2. Build the image with dockerfile in `./docker`, eg. while in `./docker` run
```bash
docker build -t bipode-httr:0.4.2 .
```
Make sure the version tag matches the version in `pyproject.toml` in the Python library

3. Tag the image with the fully qualified path (ie. includes the ACR login server). Make sure the tag is the same as `params.container` for 'az' profile in `nextflow.config`
```bash
docker tag bipode-httr:0.4.2 bnlweu57679acr04.azurecr.io/bipode-httr:0.4.2
```
4. Push image to repo in ACR
You will need the username and password from Azure portal (settings > access keys)
```bash
# login to ACR
docker login bnlweu57679acr04.azurecr.io

# push to repo
docker push bnlweu57679acr04.azurecr.io/bipode-httr:0.4.2
```

5. Check/update the container tag in Nextflow config
    - `conf/modules.config`
    - `nextflow.config` (if required)

## nf-runner (in work)
nf-runner takes Nextflow workflows from the artifact feed in Azure DevOps. In brief, __for each new version__ you will need to:
1. Push the new code to the appropriate arifact feed with the correct version
2. Update the Azure Table 'pipelines' to install the pipeline and have nf-runner pick it up and make it available

See full instructions in https://dev.azure.com/SEAC-Projects/_git/nf-runner

## Making an external release
This section covers the release of both the Python library 'bipode-httr' and the Nextflow workflow 'bipode'.

Note that version of 'bipode' and 'bipode-httr' should match.

### PyPI
Repo is https://pypi.org/project/bipode-httr

1. Increment version in `pyproject.toml` and push to 'main' branch on SERS GitHub repo (https://github.com/seacunilever/bipode-nextflow-httr)
2. For a change to the description on PyPI page update `DISCLAIMER.md`
3. Make a new release in GitHub tagged with the version number you changed in `pyproject.toml`
4. GitHub action will trigger on making a new release and push the PyPI repo. You Check progress of 'pypi' under 'Deployments' in GitHub

### Bioconda
For full instructions, see this guide https://bioconda.github.io/tutorials/2024-updating-bioinformatic-software-to-bioconda.html

1. Fork https://github.com/bioconda/bioconda-recipes in github, use your individual account not 'seacunilever' organisation
2. Make changes to `bioconda-recipes/recipes/bipode-httr/meta.yaml` (usually just version to match PyPI and sha256 sum). Make sure you generate the sha256 hash from the actual .tar.gz on PyPI, not GitHub, as it seems to make changes to it.
3. Commit and push to 'main' on your fork
4. Open a PR on https://github.com/bioconda/bioconda-recipes with 'from' being your fork (ie. choose 'compare across forks' and find it)
5. Wait for the CI to run, you can see progress on the PR
6. PRs get manually reviewed and merged

Note that Bioconda's autobump system can detect updates on PyPI and automatically create pull requests to update the corresponding recipe, but in their words it's "not guaranteed to happen immediately or reliably". Docs are at https://bioconda.github.io/contributor/updating.html. If you are quick, you shouldn't get two PRs created.

### Rebuild public container on Seqera Containers
See [./docs/containers.md](./docs/containers.md)

### Nextflow workflow
'bipode' repo in SERS GitHub (https://github.com/seacunilever/bipode-nextflow) uses 'main' branch as the latest version. To allow internal working changes, __'main' branch there should match 'release' branch in the internal repo in Azure DevOps__ (https://dev.azure.com/SEAC-Projects/BIFROST/_git/bipode-nextflow) 

- Make sure the `bipode-httr` container is also updated on Seqera Containers (see [./docs/containers.md](./docs/containers.md))
- Check/update the container tag in
    - `conf/modules.config`
    - `nextflow.config` (if required)

#### 1. Update `CHANGELOG.md`

#### 2. Push to 'release' branch in DevOps repo
Commits should be squashed, internal files removed and changed committed to 'release' branch like

```bash
# 1. Checkout the release branch
git checkout release

# 2. Merge new changes from main (without squashing yet)
git merge main --no-commit

# 3. Remove unwanted files for this release
git rm -r ./docker/
git rm README_unilever.md
```

Manually remove profile 'az' and 'local' from `nextflow.config`

```bash
# 4. Stage everything
git add .

# 5. Commit as a single release commit, change your message
git commit -m "Release 0.4.2"

# 6. Push the updated release branch
git push origin release
```

#### 3. Push to 'main' branch in SERS GitHub repo
Make sure you have a remote 'github' which points to https://github.com/seacunilever/bipode-nextflow.git

```bash
git remote -v
```

If you need to add it
```bash
git remote add github https://github.com/seacunilever/bipode-nextflow.git
```

Now push `release` branch to GitHub as `main`
```bash
git push github release:main --force
```
This pushes your local release branch to GitHub and **overwrites** the main branch there

#### 4. Make the release in GitHub
Create a new release in GitHub, the tag should match version number for the 'bipode-httr' Python package