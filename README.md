# Transformer Plot Clip

Clip GeoTIFF or LAS files according to plots without merging.

## Authors

* Christophe Schnaufer, University of Arizona, Tucson, AZ
* Max Burnette, National Supercomputing Applications, Urbana, Il

## Sample Docker Command Line
Below is a sample command line that shows how the plot clip image could be run.
An explanation of the command line options used follows.
Be sure to read up on the [docker run](https://docs.docker.com/engine/reference/run/) command line for more information.

The files that are used in this example are available through Google Drive: [plotclip_sample_data.tar.gz](https://drive.google.com/file/d/1AHx6surpXV7izII2hn3c2wn_qnj5dvCb/view?usp=sharing).

```docker run --rm --mount "src=/home/test,target=/mnt,type=bind" agdrone/transformer-plotclip:2.2 --working_space /mnt --metadata /mnt/experiment.yaml stereoTop /mnt/plots.json /mnt/orthomosaic.tif```

This example command line assumes the source files are located in the `/home/test` folder of the local machine.
The name of the Docker image to run is `agdrone/transformer-plotclip:2.2`.

We are using the same folder for the source files and the output files.
By using multiple `--mount` options, the source and output files can be separated.

**Docker commands** \
Everything between 'docker' and the name of the image are docker commands.

- `run` indicates we want to run an image
- `--rm` automatically delete the image instance after it's run
- `--mount "src=/home/test,target=/mnt,type=bind"` mounts the `/home/test` folder to the `/mnt` folder of the running image

We mount the `/home/test` folder to the running image to make files available to the software in the image.

**Image's commands** \
The command line parameters after the image name are passed to the software inside the image.
Note that the paths provided are relative to the running image (see the --mount option specified above).

- `--working_space "/mnt"` specifies the folder to use as a workspace
- `--metadata "/mnt/experiment.yaml"` is the name of the source metadata
- `stereoTop` the name of the sensor associated with the source files
- `/mnt/plots.json` the name of the GeoJSON file containing the plot geometries
- `/mnt/orthomosaic.tif` the GeoTIFF or LAS file to split by plot (in this example an TIFF file is specified) 

## Testing Source Code

Please also refer to our [Coding Standards](https://github.com/AgPipeline/Organization-info#python) for information on how we use [pylint](https://www.pylint.org/).
A pylint command line is:
```bash
# Assumes Python3.7+ is default Python version
python -m pylint --rcfile ~/agpipeline/Organization-info/pylint.rc plotclip.py
``` 

In the `tests` folder there are testing scripts and their supporting files.
The tests are designed to be run with [Pytest](https://docs.pytest.org/en/stable/).
When running the tests, the root of the repository is expected to be the starting directory.

The command line for running the tests is as follows:
```bash
# Assumes Python3.7+ is default Python version
python -m pytest -rpP
```

If test coverage reporting is desired, we suggest using [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/).
After installing this tool, the following command line will include a coverage report in the output:
```bash
# Assumes Python3.7+ is default Python version
python -m pytest --cov=. -rpP 
```

## Previous Version's Discontinued Features

- 4/1/2020: defaults to clipping RGB to the only the plot-image intersection by default; the `--full_plot_fill` command line flag restores previous default behavior
- version 2.0 merged clipped LAS files into a single file
