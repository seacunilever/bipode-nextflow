# Notes for Unilever use
This file details:
- How to run this Nextflow workflow on Unilever infrastructure
- What needs to be removed from this repo for external release on Github (https://github.com/seacunilever/bifrost)

## Running internally
Set the following env vars:

Service Principal
```bash
CLIENT_ID
CLIENT_SECRET
TENANT_ID
```

Storage account
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

Use profile 'az' in Nextflow run command

## Container from Azure Contain Registry (ACR)
The workflow is set up to pull the container from Seqera which is built from the Bioconda release of 'bifrost-httr' (see ./docs/containers.md). For internal running, the container needs to come from SERS ACR.

To update the container
- Clone repo https://dev.azure.com/SEAC-Projects/BIFROST/_git/bifrost-httr into `./docker` and checkout the branch you want
- Build with dockerfile in `./docker`
- Make sure the name is the same as `params.container` in 'az' profile in `nextflow.config`

## Checklist for external release
- Remove `./docker/`
- Remove profile 'az' from `nextflow.config`
- Remove this file
