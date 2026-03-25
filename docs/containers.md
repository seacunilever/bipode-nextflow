# seacunilever/bipode-nextflow: Container Management

## Overview

This guide explains how to create and manage containers for the Bipode pipeline using Seqera Containers. Seqera Containers provides an easy way to generate both Docker and Singularity container images from Conda packages, which can then be used in Nextflow modules.

## Prerequisites

- Access to [Seqera Containers](https://seqera.io/containers/)
- Understanding of Nextflow module structure
- Basic knowledge of Docker and Singularity containers

## Creating Containers with Seqera Containers

### Step 1: Access Seqera Containers

1. Navigate to [https://seqera.io/containers/](https://seqera.io/containers/)
2. You'll see a search interface for finding packages

### Step 2: Search for the Package

1. Enter `bipode-httr` into the search box
2. Two entries will appear:
   - One for **PyPI**
   - One for **Conda**
3. **Important**: Always select the **Conda** entry, not the PyPI one

![Selecting the Conda package for bipode-httr](images/bifrost_containers_conda.png)

### Step 3: Generate Container URIs

#### For Docker Containers

1. Select the Conda `bipode-httr` package
2. Choose **Docker** as the container type
3. Click **"Get container"**
4. The Docker URI will be presented immediately (e.g., `community.wave.seqera.io/library/bipode-httr:0.4.0--3e1755e45da93297`)
5. Copy this URI for use in your configuration

#### For Singularity Containers

1. Select the Conda `bipode-httr` package
2. Choose **Singularity** as the container type
3. Click **"Get container"**
4. **Wait for the build to complete** - the container is being built behind the scenes
5. Once ready, you'll see options for the container URI
6. **Important**: Select the **"https"** checkbox to get the HTTPS link instead of the ORAS one
7. Copy the HTTPS URI (e.g., `https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f3fddbc020a7295009575826d86375cd3a78b52a1ed859911022f2b315d723e4/data`)

![Singularity container with HTTPS option selected](images/bifrost_containers_singularity.png)

> **Note**: The first time a container is requested, it needs to be built. Subsequent requests will use the cached version and be available immediately.

## Updating Container Configuration

The container and conda configurations are centrally managed in `conf/modules.config`. This ensures consistent versions across all processes in the pipeline.

### Configuration Structure

The central configuration in `conf/modules.config` looks like this:

```nextflow
process {
    // Set container and conda for all processes
    conda = "bipode-httr=0.4.0"
    container = workflow.containerEngine == 'singularity' ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f3fddbc020a7295009575826d86375cd3a78b52a1ed859911022f2b315d723e4/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.4.0--3e1755e45da93297'

    // ... rest of process configuration
}
```

### Updating Container URIs

To update the container configuration:

1. Generate new Docker and Singularity URIs using Seqera Containers
2. Update both URIs in `conf/modules.config`
3. Update the conda version to match the container version

## Best Practices

1. **Always use the Conda package** when searching in Seqera Containers
2. **Wait for Singularity builds** to complete before copying URIs
3. **Use HTTPS URIs** for Singularity containers, not ORAS URIs
4. **Keep conda and container versions in sync** to ensure consistency
5. **Test the pipeline** after updating container URIs to ensure compatibility
6. **Version control** your changes to track container updates over time

## Troubleshooting

### Container Build Failures

- If a Singularity build fails, wait a few minutes and try again
- Check that you selected the correct Conda package, not PyPI

### Configuration Errors

- Ensure URIs are copied correctly without extra spaces or characters
- Verify that both Docker and Singularity URIs are updated
- Check that the Nextflow syntax is valid after the update
- Ensure conda version matches the container version

### Performance Issues

- New containers may take longer to pull on first use
- Consider pre-pulling containers in production environments

## Additional Resources

- [Seqera Containers Documentation](https://seqera.io/containers/)
- [Nextflow Container Documentation](https://www.nextflow.io/docs/latest/container.html)
- [Docker Documentation](https://docs.docker.com/)
- [Singularity Documentation](https://docs.sylabs.io/)
