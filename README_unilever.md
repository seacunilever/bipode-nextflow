# Notes for Unilever use
This file details:
- How to run this Nextflow workflow on Unilever infrastructure
- What needs to be removed from this repo for external release on Github (https://github.com/seacunilever/bifrost)

## Running internally
Use profile 'az' in Nextflow run command after setting the following env vars:

Service Principal
```bash
CLIENT_ID
CLIENT_SECRET
TENANT_ID
```

Storage account name
```bash
STORAGE_ACCOUNT
```

Container registry (ACR)
```bash
CONTAINER_REGISTRY_USER_NAME
CONTAINER_REGISTRY_PASSWORD
CONTAINER_REGISTRY_SERVER
```

Batch account
```bash
BATCH_ACCOUNT_NAME
BATCH_ACCOUNT_LOCATION
SUBSCRIPTION_ID
VNET
SUBNET
SUBNET_RESOURCE_GROUP_NAME
```

Check that `params.workDir` is set somewhere sensible

## Container from Azure Container Registry (ACR)
The workflow is set up to pull the container from Seqera which is built from the Bioconda release of 'bifrost-httr' (see `./docs/containers.md`). For internal running, the container needs to come from SERS ACR.

To update the container:
- Get files from branch you want from 'bifrost-httr' repo into `./docker` eg.
    - `git clone --branch <branch-name> --single-branch https://SEAC-Projects@dev.azure.com/SEAC-Projects/BIFROST/_git/bifrost-httr`
- Build with dockerfile in `./docker`
- Make sure the tag is the same as `params.container` for 'az' profile in `nextflow.config`
- Push to repo in ACR

## Make a release
- Make sure that any changes to 'bifrost-httr' Python package are published on PyPI and Bioconda
- Make sure the bifrost-httr container is also updated on Seqera Containers (see `./docs/containers.md`)
- Update the container tag in
    - `conf/modules.config`
    - `nextflow.config` (if required)

### Push to release branch
Commits should be squashed and internal files removed

```bash
# 1. Checkout the release branch
git checkout release

# 2. Merge new changes from main (without squashing yet)
git merge main --no-commit

# 3. Remove unwanted files for this release
git rm -r ./docker/
git rm ./docs/containers.md README_unilever.md
```

Manually remove profile 'az' and 'local' from `nextflow.config`

```bash
# 4. Stage everything
git add .

# 5. Commit as a single release commit, change to your message
git commit -m "Release 0.4.2"

# 6. Push the updated release branch
git push origin release
```

### Push to SERS GitHub (seacunilever/bifrost)
Make sure you have a remote 'github' which points to https://github.com/seacunilever/bifrost.git

```bash
git remote -v
```

If you need to add it
```bash
git remote add github https://github.com/seacunilever/bifrost.git
```

Now push `release` branch to GitHub as `main`
```bash
git push github release:main --force
```
This pushes your local release branch to GitHub and **overwrites** the main branch there

### Make the release
Create a new release in GitHub, the tag should match version number for the 'bifrost-httr' Python package