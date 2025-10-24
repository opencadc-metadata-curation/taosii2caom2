# taosii2caom2
A Docker image containing the tools to update and archive files, and to generate CAOM2 Observations HDF5 files from the TAOSII telescopes.

# Table of contents
1. [Set Up (One Time Only)](#set_up)
    1. [Build the Image](#build)
    1. [Credentials](#creds)
    1. [File Location](#working_dir)
1. [Initialize Workspace (One Time Only)](#initialize)
1. [Ingest Sky Camera Images](#taosii_run)
1. [Debugging](#debugging)


# How To Use taosii2caom2

These are Linux-centric instructions.

# Set Up (one time only) <a name="set_up"></a>

## Docker Image <a name="build"></a>

Set up a working directory location, on a machine with Docker installed. The machine should have POSIX access to the files to be ingested.

1. To retrieve the container image from hub.docker.com, run this:

   ```
   docker pull opencadc/taosii2caom2
   ```

## Credentials <a name="creds"></a>

1. Create a CADC User account. This user account will have authorization to read and write files and metadata records for the files at CADC.

1. Run the following command to create a proxy certificate file for that CADC user, in the working directory. You will be prompted for the password:

   ```
   docker run --rm -ti --user $(id -u):$(id -g) -e HOME=/usr/src/app -v ${PWD}:/usr/src/app opencadc/taosii2caom2 cadc-get-cert --cert-filename /usr/src/app/cadcproxy.pem --days-valid 10 -u <CADC User Name here>
   ```

1. The `taosii_run*.sh` scripts described later will attempt to copy `$HOME/.ssl/cadcproxy.pem` to the 'working directory', so copy `cadcproxy.pem` to `$HOME/.ssl/`. 

## File Location <a name="working_dir"></a>

`taosii2caom2` can store files from local disk to CADC storage. This behaviour is controlled by configuration
information, located in a file named `config.yml`. Most of the `config.yml` values are already set appropriately, but there are a few values that need to be
set according to the execution environment. For a complete description of the `config.yml` content, see
https://github.com/opencadc/collection2caom2/wiki/config.yml.

1. Copy the file `config.yml` to the working directory. e.g.:

   ```
   wget https://github.com/opencadc-metadata-curation/taosii2caom2/blob/5a5e79a33677c15cbcf15975b0ddb85d5c6892f4/config/config.yml
   ```
1. In `config.yml`, set the value of `tap_id` to `ivo://cadc.nrc.ca/ams/shared`

1. In `config.yml`, tell the application what to do with files on disk after the files have been stored to CADC:

   1. Set `cleanup_files_when_storing` to `True` or `False`. 
       1. If this is set to `False`, `taosii2caom2` will do nothing with the files on disk.
       1. If this is set to `True`, `taosii2caom2` will move stored files to either a success or failure location.

   2. If `cleanup_files_when_storing` is set to `True`, set `cleanup_failure_destination`  and `cleanup_success_destination` to fully-qualified directory names that are visible within the Docker container. A directory is visible within a Docker container if it
      is one of the values on the right-hand-side of the colon in a `-v` `docker run` parameter.

1. Tell `taosii2caom` in `config.yml` whether to re-submit duplicate files.
   1. Set `store_modified_files_only` to `True` or `False`. If this is set to `False`, there is no effect on execution. If this is set to true, `taosii2caom2`
      checks that the local version of the file has a md5 checksum that is different from the file at CADC before transferring the file to CADC storage. This affects only the `store` `task_types`.


# Initialize Execution Location (one time only) <a name="initialize"></a>

1. In the main branch of this repository, find the scripts directory, and copy the files `taosii_run.sh`  and `taosii_run_incremental.sh` to the working directory. e.g.:

   ```
   wget https://github.com/opencadc-metadata-curation/taosii2caom2/blob/5a5e79a33677c15cbcf15975b0ddb85d5c6892f4/scripts/taosii_run.sh
   wget https://github.com/opencadc-metadata-curation/taosii2caom2/blob/5a5e79a33677c15cbcf15975b0ddb85d5c6892f4/scripts/taosii_run_incremental.sh
   ```

1. Ensure the scripts are executable:

   ```
   chmod +x taosii_run.sh
   chmod +x taosii_run_incremental.sh
   ```

1. Edit the scripts to specify the file location:

   1. `taosii_run.sh`:
      1. Find this line: `docker run --rm --name ${COLLECTION}_todo  --user $(id -u):$(id -g) -e HOME=/usr/src/app -v ${PWD}:/usr/src/app/ -v /data:/data ${IMAGE} ${COLLECTION}_run || exit $?`
      2. Replace the `/data/:` portion of the command with the fully-qualified directory name of where the application should find the data. This directory will be called the "data source directory" in these instructions. 

   1. `taosii_run_incremental.sh`:
      1. Find this line: `docker run --rm --name ${COLLECTION}_state  --user $(id -u):$(id -g) -e HOME=/usr/src/app -v ${PWD}:/usr/src/app/ -v /data:/data ${IMAGE} ${COLLECTION}_run_incremental || exit $?`
      2. Replace the `/data/:` portion of the command with the fully-qualified directory name of where the application should find the data.

# How to Store Files and Create CAOM2 Records at CADC <a name="taosii_run"></a>

`taosiicaom2` may be run so that it processes files incrementally, according to their timestamp on disk, or so that is processes all the files it finds.

1. To run the application incrementally:

   ```
   ./taosii_run_incremental.sh
   ```
   By default, incremental mode will start 24 hours prior to the current execution time. This can be changed by modifying the `state.yml` file content that is created on the first run.

1. To run the application on all the files it finds:

    ```
    ./taosii_run.sh
    ```

# Debugging <a name="debugging"></a>

1. To debug the application from inside the container, run the following command. Replace the `<data directory here>` with the fully-qualified path name of the directory where the data to be processed is located.

   ```
   user@dockerhost:<cwd># docker run --rm -ti -v ${PWD}:/usr/src/app -v <data directory here>:/data --user $(id -u):$(id -g) -e HOME=/usr/src/app --name taosii_run taosii2caom2_app /bin/bash
   cadcops@53bef30d8af3:/usr/src/app# taosii_run
   ```

1. For some instructions that might be helpful on using containers, see:
   https://github.com/opencadc/collection2caom2/wiki/Docker-and-Collections

1. For some insight into what's happening, see: https://github.com/opencadc/collection2caom2

1. For Docker information, see: https://www.docker.com
