[![license](https://img.shields.io/badge/license-BSD%203-green?logo=Open-Source-Initiative)](https://github.com/AgPipeline/transformer-soilmask/blob/master/LICENSE)

[![Enforcing testing](https://github.com/AgPipeline/transformer-soilmask/workflows/Enforcing%20testing/badge.svg)](https://github.com/AgPipeline/transformer-soilmask/actions?query=workflow%3A%22Enforcing+testing%22)
[![testing](https://github.com/AgPipeline/transformer-soilmask/workflows/Testing%20Docker%20image/badge.svg)](https://github.com/AgPipeline/transformer-soilmask/actions?query=workflow%3A%22Testing+Docker+image%22)

# Transformer Plot Clip

Clip GeoTIFF or LAS files according to plots without merging.

## Authors

* Christophe Schnaufer, University of Arizona, Tucson, AZ
* Max Burnette, National Supercomputing Applications, Urbana, Il

## Use 

### Sample Docker Command Line

First build the Docker image, using the Dockerfile, and tag it agdrone/transformer-plotclip:2.2. 
Read about the [docker build](https://docs.docker.com/engine/reference/commandline/build/) command if needed.

```bash
docker build -t agdrone/transformer-plotclip:2.2 ./
```

Below is a sample command line that shows how the plot clip image could be run.
An explanation of the command line options used follows.
Be sure to read up on the [docker run](https://docs.docker.com/engine/reference/run/) command line for more information.

The files that are used in this example are available for downloading [plotclip_sample_data.tar.gz](https://de.cyverse.org/dl/d/4AF1E272-6A28-400D-B102-E0F6F168BA10/plotclip_sample_data.tar.gz).
The following can be used to download and extract the contents of the test data archive:
```bash
mkdir test_data
curl -X GET https://de.cyverse.org/dl/d/4AF1E272-6A28-400D-B102-E0F6F168BA10/plotclip_sample_data.tar.gz -o test_data/plotclip_sample_data.tar.gz
tar -xzvf test_data/plotclip_sample_data.tar.gz -C test_data/
```

The following command runs the plotclip Docker image:
```bash
docker run --rm --mount "src=/home/test_data,target=/mnt,type=bind" agdrone/transformer-plotclip:2.2 --working_space /mnt --metadata /mnt/experiment.yaml /mnt/plots.json /mnt/orthomosaic.tif
```

This example command line assumes the source files are located in the `/home/test_data` folder of the local machine.
The name of the Docker image to run is `agdrone/transformer-plotclip:2.2`.

We are using the same folder for the source files and the output files.
By using multiple `--mount` options, the source and output files can be separated.

**Docker commands** \
Everything between 'docker' and the name of the image are docker commands.

- `run` indicates we want to run an image
- `--rm` automatically delete the image instance after it's run
- `--mount "src=/home/test,target=/mnt,type=bind"` mounts the `/home/test_data` folder to the `/mnt` folder of the running image

We mount the `/home/test_data` folder to the running image to make files available to the software in the image.

**Image's commands** \
The command line parameters after the image name are passed to the software inside the image.
Note that the paths provided are relative to the running image (see the --mount option specified above).

- `--working_space "/mnt"` specifies the folder to use as a workspace
- `--metadata "/mnt/experiment.yaml"` is the name of the source metadata
- `/mnt/plots.json` the name of the GeoJSON file containing the plot geometries
- `/mnt/orthomosaic.tif` the GeoTIFF or LAS file to split by plot (in this example an TIFF file is specified) 

## Acceptance Testing

There are automated test suites that are run via [GitHub Actions](https://docs.github.com/en/actions).
In this section we provide details on these tests so that they can be run locally as well.

These tests are run when a [Pull Request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) or [push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository) occurs on the `develop` or `master` branches.
There may be other instances when these tests are automatically run, but these are considered the mandatory events and branches.

### PyLint and PyTest

These tests are run against any Python scripts that are in the repository.

[PyLint](https://www.pylint.org/) is used to both check that Python code conforms to the recommended coding style, and checks for syntax errors.
The default behavior of PyLint is modified by the `pylint.rc` file in the [Organization-info](https://github.com/AgPipeline/Organization-info) repository.
Please also refer to our [Coding Standards](https://github.com/AgPipeline/Organization-info#python) for information on how we use [pylint](https://www.pylint.org/).

The following command can be used to fetch the `pylint.rc` file:
```bash
wget https://raw.githubusercontent.com/AgPipeline/Organization-info/master/pylint.rc
```

Assuming the `pylint.rc` file is in the current folder, the following command can be used against the `plotclip.py` file:
```bash
# Assumes Python3.7+ is default Python version
python -m pylint --rcfile ./pylint.rc plotclip.py
``` 

In the `tests` folder there are testing scripts and their supporting files.
The tests are designed to be run with [Pytest](https://docs.pytest.org/en/stable/).
When running the tests, the root of the repository is expected to be the starting directory.

The command line for running the tests is as follows:
```bash
# Assumes Python3.7+ is default Python version
python -m pytest -rpP
```

If [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/) is installed, it can be used to generate a code coverage report as part of running PyTest.
The code coverage report shows how much of the code has been tested; it doesn't indicate **how well** that code has been tested.
The modified PyTest command line including coverage is:
```bash
# Assumes Python3.7+ is default Python version
python -m pytest --cov=. -rpP 
```

### Docker Testing

The Docker testing Workflow replicate the examples in this document to ensure they continue to work.

## Previous Version's Discontinued Features

- 12/22/2020: removed previously required sensor parameter since it wasn't used
- 4/1/2020: defaults to clipping RGB to the only the plot-image intersection by default; the `--full_plot_fill` command line flag restores previous default behavior
- version 2.0 merged clipped LAS files into a single file
