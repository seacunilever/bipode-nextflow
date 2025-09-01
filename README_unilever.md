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

## Container from Azure Contain Registry (ACR)
The workflow is set up to pull the container from Seqera which is built from the Bioconda release of 'bifrost-httr' (see ./docs/containers.md). For internal running, the container needs to come from SERS ACR.

To update the container:
- Get files from branch you want from 'bifrost-httr' repo into `./docker` eg.
    - `git clone --branch <branch-name> --single-branch https://SEAC-Projects@dev.azure.com/SEAC-Projects/BIFROST/_git/bifrost-httr`
- Build with dockerfile in `./docker`
- Make sure the tag is the same as `params.container` for 'az' profile in `nextflow.config`
- Push to repo in ACR

## Checklist for external release
- Remove `./docker/`
- Remove profile 'az' from `nextflow.config`
- Remove this file
