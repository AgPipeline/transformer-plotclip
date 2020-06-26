"""Contains transformer configuration information
"""
from agpypeline.configuration import Configuration

class ConfigurationPlotclip(Configuration):
    # Silence this error until we have public methods
    # pylint: disable=too-few-public-methods
    # The version number of the transformer
    transformer_version = '2.0'

    # The transformer description
    transformer_description = 'TERRA-REF plot clipper'

    # Short name of the transformer
    transformer_name = 'terra.plotclipper'

    # The sensor associated with the transformer
    transformer_sensor = 'scanner3DTop'

    # The transformer type (eg: 'rgbmask', 'plotclipper')
    transformer_type = 'plotclipper'

    # The name of the author of the extractor
    author_name = 'Chris Schnaufer'

    # The email of the author of the extractor
    author_email = 'schnaufer@email.arizona.edu'

    # Contributors to this transformer
    contributors = ["Max Burnette"]

    # Reposity URI of where the source code lives
    repository = 'https://github.com/AgPipeline/transformer-plotclip.git'
