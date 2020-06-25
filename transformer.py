"""Performs plot-level clipping on geo referenced files
"""

import argparse
import copy
import datetime
import json
import logging
import os
import re
import subprocess
from typing import Optional
import numpy as np
from osgeo import gdal, ogr, osr
import liblas

import configuration
import transformer_class


class __internal__:
    """Class for internal use only functions
    """

    def __init__(self):
        """Initializes class instance
        """

    @staticmethod
    def get_sr_from_crs(crs: dict) -> osr.SpatialReference:
        """Returns a SpatialReference determined from the crs GeoJSON passed in
        Arguments:
            crs: the associated crs in GeoJson format
        Exceptions:
            Raises RuntimeError if a problem with the CRS is found
        """
        ogc_crs_prefix = 'urn:ogc:def:crs:'

        if 'type' not in crs:
            raise RuntimeError('GeoJSON crs type not specified "%s"' % str(crs))
        if crs['type'] == 'name':
            # eg: { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::32612" } }
            if 'properties' not in crs or 'name' not in crs['properties']:
                raise RuntimeError('Unknown named format specification in GeoJSON crs "%s"' % str(crs))
            sr_name = crs['properties']['name'][len(ogc_crs_prefix):]
            sr_type, sr_def = sr_name.split(':', 1)
            sr_type = sr_type.lower()
            sr_def = sr_def.strip(':')
        else:
            # eg: { "type": "EPSG", "properties": { "code": "4326" } }
            if 'properties' not in crs or 'code' not in crs['properties']:
                raise RuntimeError('Unknown EPSG format specification in GeoJSON crs "%s"' % str(crs))
            sr_type = crs['type'].lower()
            sr_def = crs['properties']['code']

        return_sr = osr.SpatialReference()
        if sr_type == 'epsg':
            if return_sr.ImportFromEPSG(int(sr_def)) != ogr.OGRERR_NONE:
                raise RuntimeError('Unable to load EPSG code "%s" from GeoJSON crs %s' % (sr_def, str(crs)))
        elif sr_type == 'ogc':
            if return_sr.ImportFromEPSG(4326) != ogr.OGRERR_NONE:
                raise RuntimeError('Unable to load default OGC from GeoJSON crs %s' % str(crs))
        else:
            raise RuntimeError('Spatial reference type is currently not supported: %s' % str(sr_type))

        return return_sr

    @staticmethod
    def get_geojson_file_sr(geojson: dict) -> Optional[osr.SpatialReference]:
        """Returns the spatial reference loaded from the GeoJSON
        Arguments:
            geojson: the loaded GeoJSON to find the CRS of
        Return:
            Returns a spatial reference of any specified CRS, or the default EPGS:4236 spatial reference
        Notes:
            Refer to the GeoJSON specification for information on CRS and assumed coordinate systems
        """
        if 'crs' in geojson:
            file_sr = __internal__.get_sr_from_crs(geojson['crs'])
            if not file_sr:
                logging.error('Unable to load CRS "%s"', str(geojson['crs']))
                return None
        else:
            file_sr = osr.SpatialReference()
            if file_sr.ImportFromEPSG(4326) != ogr.OGRERR_NONE:
                logging.error('Unable to load default SpatialReference 4326')
                return None

        return file_sr

    @staticmethod
    def get_plot_key_name(properties: dict, default_key: str = None) -> Optional[tuple]:
        """ Attempts to find the plot name from the set of properties
        Arguments:
            properties: a dictionary that's searched for well known plot name keys
            default_key: an optional key to check for first; if found it's returned otherwise a search is performed
        Return:
            The found plot key name and current plot name, or None if nothing useable was found
        """
        # Check if the key exists
        if default_key and default_key in properties:
            logging.debug('[get_plot_key_name] Default key "%s": "%s"', default_key, properties[default_key])
            return default_key, properties[default_key]

        # Search the dictionary
        if 'observationUnitName' in properties:
            logging.debug('[get_plot_key_name] observationUnitName: "%s"', default_key, properties['observationUnitName'])
            return 'observationUnitName', properties['observationUnitName']

        best_fit = None
        for one_key in properties:
            if (one_key.lower().find('plot') >= 0) and ((one_key.lower().find('name') >= 0) or one_key.lower().find('id') >= 0):
                logging.debug('[get_plot_key_name]  1 best fit "%s"', one_key)
                best_fit = one_key
            elif one_key.lower() == 'id' and not best_fit:
                logging.debug('[get_plot_key_name]  2 best fit "%s"', one_key)
                best_fit = one_key

        if best_fit:
            logging.debug('[get_plot_key_name] %s: "%s"', best_fit, properties[best_fit])
            return best_fit, properties[best_fit]

        return None

    @staticmethod
    def load_plotfile(plot_file: str, plot_column: str = None) -> dict:
        """Loads the GeoJSON plot file and returns a dict with plot names and geometries as key, value pairs
        Arguments:
            plot_file: the path of the GeoJSON file to load
            plot_column: optional parameter specifying the name of column containing the plot names
        Return:
            A dict containing the loaded plots with their associated geometries as ogr.Geometry
        Exceptions:
            Raises RuntimeError if a problem is found with the loaded JSON
            Other exceptions may be raised when loading and parsing the JSON
        """
        plots = {}

        # Load the file contents and check them
        with open(plot_file, 'r') as in_file:
            geojson = json.load(in_file)
            if not geojson:
                raise RuntimeError('No JSON was found in file: "%s"' % plot_file)
            for req_key in ['type', 'features']:
                if req_key not in geojson:
                    raise RuntimeError('Missing GeoJSON key "%s" detected in file: "%s"' % (req_key, plot_file))

        # Determine if there's a CRS and create a osr.SpatialReference for it.
        # Otherwise create a EPSG:4326 osr.SpatialReference to use
        file_sr = __internal__.get_geojson_file_sr(geojson)
        if not file_sr:
            raise RuntimeError('Unable to load CRS for file "%s"' % plot_file)

        # Loop through the features
        feature_idx = 0
        plot_key = None
        logging.debug("Have %s features", str(len(geojson['features'])))
        for one_feature in geojson['features']:
            # Initialize for each pass
            feature_idx += 1
            plot_name = None
            plot_key = plot_column if plot_column else plot_key

            # Look for values we need and get those values
            if 'type' not in one_feature or one_feature['type'] != 'Feature' or 'geometry' not in one_feature:
                logging.info('Skipping unknown feature at index %s: "%s"', str(feature_idx), str(plot_file))
                continue
            if 'properties' in one_feature:
                plot_key, plot_name = __internal__.get_plot_key_name(one_feature['properties'], plot_key)
            if not plot_name:
                plot_name = 'Plot ' + str(feature_idx)

            # Create the geometry
            plot_geom = ogr.CreateGeometryFromJson(json.dumps(one_feature['geometry']))
            if not plot_geom:
                raise RuntimeError('Unable to create geometry from JSON at index %s: "%s"' % (str(feature_idx), plot_file))

            # Assign the spatial reference if needed
            geom_sr = plot_geom.GetSpatialReference()
            if not geom_sr or not geom_sr.IsSame(file_sr):
                plot_geom.AssignSpatialReference(file_sr)

            # Store the plot
            logging.debug("Plot name: '%s'", plot_name)
            plots[plot_name] = plot_geom

        return plots

    @staticmethod
    def convert_geometry(geometry: ogr.Geometry, new_spatialreference: osr.SpatialReference) -> ogr.Geometry:
        """Converts the geometry to the new spatial reference if possible
        Arguments:
            geometry - The geometry to transform
            new_spatialreference - The spatial reference to change to
        Return:
            The transformed geometry or the original geometry. If either the
            new Spatial Reference parameter is None, or the geometry doesn't
            have a spatial reference, then the original geometry is returned.
        """
        if not new_spatialreference or not geometry:
            return geometry

        return_geometry = geometry
        try:
            geom_sr = geometry.GetSpatialReference()
            if geom_sr and not new_spatialreference.IsSame(geom_sr):
                transform = osr.CreateCoordinateTransformation(geom_sr, new_spatialreference)
                new_geom = geometry.Clone()
                if new_geom:
                    new_geom.Transform(transform)
                    return_geometry = new_geom
        except Exception as ex:
            logging.warning("Exception caught while transforming geometries: %s", str(ex))
            logging.warning("    Returning original geometry")

        return return_geometry

    @staticmethod
    def find_plots_intersect_boundingbox(bounding_box: ogr.Geometry, all_plots: dict) -> dict:
        """Take a list of plots and return only those overlapping bounding box.
        Arguments:
            bounding_box: the geometry of the bounding box
            all_plots: the dictionary of all available plots
        Return:
            A dictionary of all intersecting plots
        """
        bb_sr = bounding_box.GetSpatialReference()
        intersecting_plots = {}
        logging.debug("[find_plots_intersect_boundingbox] Bounding box %s %s", str(bb_sr), str(bounding_box))

        for plot_name in all_plots:
            current_poly = all_plots[plot_name]

            # Check for a need to convert coordinate systems
            check_poly = current_poly
            if bb_sr:
                poly_sr = current_poly.GetSpatialReference()
                if poly_sr and not bb_sr.IsSame(poly_sr):
                    # We need to convert to the same coordinate system before an intersection
                    check_poly = __internal__.convert_geometry(current_poly, bb_sr)

            logging.debug("[find_plots_intersect_boundingbox] Intersection with %s", str(check_poly))
            intersection_with_bounding_box = bounding_box.Intersection(check_poly)

            if intersection_with_bounding_box is not None:
                intersection = json.loads(intersection_with_bounding_box.ExportToJson())
                if 'coordinates' in intersection and len(intersection['coordinates']) > 0:
                    intersecting_plots[str(plot_name)] = current_poly

        return intersecting_plots

    @staticmethod
    def image_get_geobounds(filename: str) -> list:
        """Uses gdal functionality to retrieve recilinear boundaries from the file

        Args:
            filename(str): path of the file to get the boundaries from

        Returns:
            The upper-left and calculated lower-right boundaries of the image in a list upon success.
            The values are returned in following order: min_y, max_y, min_x, max_x. A list of numpy.nan
            in each position is returned if the boundaries can't be determined
        """
        try:
            src = gdal.Open(filename)
            ulx, xres, _, uly, _, yres = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)

            min_y = min(uly, lry)
            max_y = max(uly, lry)
            min_x = min(ulx, lrx)
            max_x = max(ulx, lrx)

            return [min_y, max_y, min_x, max_x]
        except Exception as ex:
            logging.info("[image_get_geobounds] Exception caught: %s", str(ex))
            if logging.getLogger().level == logging.DEBUG:
                logging.exception("[image_get_geobounds] Exception")

        return [np.nan, np.nan, np.nan, np.nan]

    @staticmethod
    def clip_raster(rast_path: str, bounds: tuple, out_path: str = None, compress: bool = False) ->\
            Optional[np.ndarray]:
        """Clip raster to polygon.
        Arguments:
            rast_path: path to raster file
            bounds: (min_y, max_y, min_x, max_x)
            out_path: if provided, where to save as output file
            compress: set to True to compress the image and False to leave it uncompressed
        Returns:
            The array of clipped pixels
        Notes:
            Oddly, the "features path" can be either a filename
            OR a geojson string. GDAL seems to figure it out and do
            the right thing.
            From http://karthur.org/2015/clipping-rasters-in-python.html
        """
        if not out_path:
            out_path = "temp.tif"

        # Clip raster to GDAL and read it to numpy array
        coords = "%s %s %s %s" % (bounds[2], bounds[1], bounds[3], bounds[0])
        if compress:
            cmd = 'gdal_translate -projwin %s "%s" "%s"' % (coords, rast_path, out_path)
        else:
            cmd = 'gdal_translate -co COMPRESS=LZW -projwin %s "%s" "%s"' % (coords, rast_path, out_path)
        subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        out_px = np.array(gdal.Open(out_path).ReadAsArray())

        if np.count_nonzero(out_px) > 0:
            if out_path == "temp.tif":
                os.remove(out_path)
            return out_px

        os.remove(out_path)
        return None

    @staticmethod
    def get_epsg(filename: str) -> Optional[str]:
        """Returns the EPSG of the georeferenced image file
        Args:
            filename(str): path of the file to retrieve the EPSG code from
        Return:
            Returns the found EPSG code, or None if it's not found or an error ocurred
        """
        logger = logging.getLogger(__name__)

        try:
            src = gdal.Open(filename)

            proj = osr.SpatialReference(wkt=src.GetProjection())

            return proj.GetAttrValue('AUTHORITY', 1)
        except Exception as ex:
            logger.warning("[get_epsg] Exception caught: %s", str(ex))
            if logging.getLogger().level == logging.DEBUG:
                logging.exception("[get_epsg] Exception")

        return None

    @staticmethod
    def get_las_epsg_from_header(header: liblas.header.Header) -> str:
        """Returns the found EPSG code from the LAS header
        Arguments:
            header: the loaded LAS header to find the SRID in
        Return:
            Returns the SRID as a string if found, None is returned otherwise
        """
        epsg = None
        search_terms_ordered = ['DATUM', 'AUTHORITY', '"EPSG"', ',']
        try:
            # Get the WKT from the header, find the DATUM, then finally the EPSG code
            srs = header.get_srs()
            wkt = srs.get_wkt().decode('UTF-8')
            idx = -1
            for term in search_terms_ordered:
                idx = wkt.find(term)
                if idx < 0:
                    break
            if idx >= 0:
                epsg = re.search(r'\d+', wkt[idx:])[0]
        except Exception as ex:
            logging.debug("Unable to find EPSG in LAS file header")
            logging.debug("    exception caught: %s", str(ex))

        return epsg

    @staticmethod
    def get_las_extents(file_path: str, default_epsg: int = None) -> Optional[ogr.Geometry]:
        """Calculate the extent of the given las file and return as a ogr.Geometry.
        Arguments:
            file_path: path to the file from which to load the bounds
            default_epsg: the default EPSG to assume if a file has a boundary but not a coordinate system
        Return:
            Returns the geometry representing the image boundary, or None if the
            bounds could not be loaded
        Notes:
            If a file doesn't have a coordinate system and a default epsg is specified, the
            return geometry will use the default_epsg.
            If a file doesn't have a coordinate system and there isn't a default epsg specified, the boundary
            of the image is not returned (None) and a warning is logged.
        """
        # Get the bounds and the EPSG code
        las_info = liblas.file.File(file_path, mode='r')
        min_bound = las_info.header.min
        max_bound = las_info.header.max
        epsg = __internal__.get_las_epsg_from_header(las_info.header)
        if epsg is None:
            if default_epsg is not None:
                epsg = default_epsg
            else:
                logging.warning("Unable to find EPSG and not default is specified for file '%s'", file_path)
                return None

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(min_bound[1], min_bound[0])  # Upper left
        ring.AddPoint(min_bound[1], max_bound[0])  # Upper right
        ring.AddPoint(max_bound[1], max_bound[0])  # lower right
        ring.AddPoint(max_bound[1], min_bound[0])  # lower left
        ring.AddPoint(min_bound[1], min_bound[0])  # Closing the polygon

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        ref_sys = osr.SpatialReference()
        if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
            poly.AssignSpatialReference(ref_sys)
            return poly

        logging.error("Failed to import EPSG %s for las file %s", str(epsg), file_path)
        return None

    @staticmethod
    def clip_las(las_path: str, clip_tuple: tuple, out_path: str) -> None:
        """Clip LAS file to polygon.
        Arguments:
          las_path: path to point cloud file
          clip_tuple: tuple containing (minX, maxX, minY, maxY) of clip bounds
          out_path: output file to write
        Notes:
            The clip_tuple is assumed to be in the correct coordinate system for the point cloud file
        """
        bounds_str = "([%s, %s], [%s, %s])" % (clip_tuple[0], clip_tuple[1], clip_tuple[2], clip_tuple[3])

        pdal_dtm = out_path.replace(".las", "_dtm.json")
        with open(pdal_dtm, 'w') as dtm:
            dtm_data = """{
                "pipeline": [
                    "%s",
                    {
                        "type": "filters.crop",
                        "bounds": "%s"
                    },
                    {
                        "type": "writers.las",
                        "filename": "%s"
                    }
                ]
            }""" % (las_path, bounds_str, out_path)
            logging.debug("Writing dtm file contents: %s", str(dtm_data))
            dtm.write(dtm_data)

        cmd = 'pdal pipeline "%s"' % pdal_dtm
        logging.debug("Running pipeline command: %s", cmd)
        subprocess.call([cmd], shell=True)
        os.remove(pdal_dtm)

    @staticmethod
    def get_image_bounds(file_path: str, default_epsg: int = None) -> Optional[ogr.Geometry]:
        """Loads the boundaries of the image file and returns the ogr.Geometry
           representing the bounds (including EPSG code)
        Arguments:
            file_path: path to the file from which to load the bounds
            default_epsg: the default EPSG to assume if a file has a boundary but not a coordinate system
        Return:
            Returns the geometry representing the image boundary, or None if the
            bounds could not be loaded
        Notes:
            If a file doesn't have a coordinate system and a default epsg is specified, the
            return geometry will use the default_epsg.
            If a file doesn't have a coordinate system and there isn't a default epsg specified, the boundary
            of the image is not returned (None) and a warning is logged.
        """
        # Get the bounds (if they exist)
        bounds = __internal__.image_get_geobounds(file_path)
        if bounds[0] == np.nan:
            return None

        epsg = __internal__.get_epsg(file_path)
        if epsg is None:
            if default_epsg:
                epsg = default_epsg
            else:
                logging.warning("Files does not have a coordinate system defined and no default was specified: '%s'",
                                file_path)
                return None

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(bounds[2], bounds[1])  # Upper left
        ring.AddPoint(bounds[3], bounds[1])  # Upper right
        ring.AddPoint(bounds[3], bounds[0])  # lower right
        ring.AddPoint(bounds[2], bounds[0])  # lower left
        ring.AddPoint(bounds[2], bounds[1])  # Closing the polygon

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        ref_sys = osr.SpatialReference()
        if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
            poly.AssignSpatialReference(ref_sys)
            return poly

        logging.error("Failed to import EPSG %s for image file %s", str(epsg), file_path)
        return None

    @staticmethod
    def get_files_to_process(file_list: list, sensor: str, default_epsg: int = None) -> dict:
        """Returns a dictionary of georeferenced files to process
        Arguments:
            file_list: the list of file paths to process
            sensor: the name of the sensor associated with the files
            default_epsg: the default EPSG value to use if a file is missing one
        Return:
            Returns a dictionary with the file names as keys. Each key's value is another dictionary containing
            the file path, file bounds (as ogr.Geometry), and the sensor name
        """
        files_to_process = {}
        for one_file in file_list:
            filename = os.path.basename(one_file)
            if filename in files_to_process:
                continue
            if not os.path.exists(one_file):
                logging.warning("Skipping file that does not exist: '%s'", one_file)
                continue

            if one_file.endswith('.tif'):
                files_to_process[filename] = {
                    'path': one_file,
                    'bounds': __internal__.get_image_bounds(one_file, default_epsg),
                    'sensor_name': sensor
                }
            elif one_file.endswith(".las"):
                files_to_process[filename] = {
                    'path': one_file,
                    'bounds': __internal__.get_las_extents(one_file, default_epsg),
                    'sensor_name': sensor
                }
        return files_to_process

    @staticmethod
    def geometry_to_tuples(geom: ogr.Geometry) -> tuple:
        """Returns the bounds of the shape
        Arguments:
            geom: the geometry to return the bounds of
        Return:
            A tuple containing the bounds in (min Y, max Y, min X, max X) order
        """
        current_env = geom.GetEnvelope()

        return current_env[2], current_env[3], current_env[0], current_env[1]

    @staticmethod
    def calculate_overlap_percent(check_bounds: ogr.Geometry, other_bounds: ogr.Geometry) -> float:
        """Calculates and returns the percentage overlap between the two boundaries.
           The calculation determines the overlap shape between the two parameters and
           then calculates the percentage by dividing the overlap area by the checking
           bounds area, and returns that value.
        Args:
            check_bounds: geometry of boundary to check
            other_bounds: geometry of boundary to check against
        Return:
            The calculated overlap percent (0.0 - 1.0) or 0.0 if there is no overlap.
            If an exception is detected, a warning message is logged and 0.0 is returned.
        """
        try:
            if check_bounds and other_bounds:
                intersection = other_bounds.Intersection(check_bounds)
                if intersection:
                    return intersection.Area() / check_bounds.Area()
        except Exception as ex:
            logging.warning("Exception caught while calculating shape overlap: %s", str(ex))

        return 0.0

    @staticmethod
    def clip_raster_intersection(file_path: str, file_bounds: ogr.Geometry, plot_bounds: ogr.Geometry, out_file: str) ->\
            Optional[int]:
        """Clips the raster to the intersection of the file bounds and plot bounds
        Arguments:
            file_path: the path to the source file
            file_bounds: the geometric boundary of the source file
            plot_bounds: the geometric boundary of the plot to clip to
            out_file: the path to store the clipped image
        Return:
            The number of pixels in the new image, or None if no pixels were saved
        Notes:
            Assumes the boundaries are in the same coordinate system
        Exceptions:
            Raises RuntimeError if the polygons are invalid
        """
        logging.debug("Clip to intersect of plot boundary: File: '%s' '%s' Plot: '%s'", file_path, str(file_bounds), str(plot_bounds))
        try:
            if not file_bounds or not plot_bounds:
                logging.error("Invalid polygon specified for clip_raster_intersection: File: '%s' plot: '%s'",
                              str(file_bounds), str(plot_bounds))
                raise RuntimeError("One or more invalid polygons specified when clipping raster")

            intersection = file_bounds.Intersection(plot_bounds)
            if not intersection or not intersection.Area():
                logging.info("File does not intersect plot boundary: %s", file_path)
                return None

            # Make sure we pass a multipolygon down to the tuple converter
            if intersection.GetGeometryName().startswith('MULTI'):
                multi_polygon = intersection
            else:
                multi_polygon = ogr.Geometry(ogr.wkbMultiPolygon)
                multi_polygon.AddGeometry(intersection)

            # Proceed to clip to the intersection
            tuples = __internal__.geometry_to_tuples(multi_polygon)
            return __internal__.clip_raster(file_path, tuples, out_path=out_file, compress=True)

        except Exception as ex:
            logging.exception("Exception caught while clipping image to plot intersection")
            raise ex

    @staticmethod
    def cleanup_request_md(source_md: dict) -> dict:
        """Makes a copy of the source metadata and cleans it up for use as plot-level information
        Arguments:
            source_md: the source metadata to clone and clean up
        Returns:
            returns the cleaned up metadata
        """
        if not source_md:
            return {}

        new_md = copy.deepcopy(source_md)
        new_md.pop('list_files', None)
        new_md.pop('context_md', None)
        new_md.pop('working_folder', None)

        return new_md

    @staticmethod
    def prepare_container_md(plot_name: str, plot_md: dict, sensor: str, source_file: str, result_files: list) -> dict:
        """Prepares the metadata for a single container
        Arguments:
            plot_name: the name of the container
            plot_md: the metadata associated with this container
            sensor: the name of the related sensor
            source_file: the name of the source file
            result_files: list of files to add to container metadata
        Return:
            The formatted metadata
        Notes:
            The files in result_files are checked for existence before being added to the metadata
        """
        cur_md = {
            'name': plot_name,
            'metadata': {
                'replace': True,
                'data': plot_md
            },
            'file': []
        }
        for one_file in result_files:
            if os.path.exists(one_file):
                cur_md['file'].append({
                    'path': one_file,
                    'key': sensor,
                    'metadata': {
                        'source': source_file,
                        'transformer': configuration.TRANSFORMER_NAME,
                        'version': configuration.TRANSFORMER_VERSION,
                        'timestamp': datetime.datetime.utcnow().isoformat(),
                        'plot_name': plot_name
                    }
                })
        return cur_md

    @staticmethod
    def merge_container_md(dest_md: list, new_md: dict) -> list:
        """Merges container level metadata ensuring there aren't any plot-level
           duplicates or duplicate file entries for a plot entry
        Arguments:
            dest_md: the list of current metadata to merge into
            new_md: the new metadata to merge
        Return:
            Returns a new list of metadata with the new metadata merged into it
        """
        # Return something meaningful if we have missing or empty dict
        if not dest_md:
            if new_md:
                return [new_md]
            return []

        # Copy the metadata and look for a match
        match_idx = -1
        md_len = len(dest_md)
        for idx in range(0, md_len):
            if dest_md[idx]['name'] == new_md['name']:
                match_idx = idx
                break

        # If no match found, add and return
        if match_idx == -1:
            dest_md.append(new_md)
            return dest_md

        # Merge the metadata
        working_md = dest_md[match_idx]
        if 'files' in new_md:
            if 'files' in working_md:
                # Only add files that aren't included in the destination metadata already
                for one_file in new_md['files']:
                    file_match_found = False
                    for match_file in working_md['files']:
                        if one_file['path'] == match_file['path']:
                            file_match_found = True
                            break
                    if not file_match_found:
                        dest_md[match_idx]['files'].append(one_file)
            else:
                # Target metadata doesn't have a 'files' entry
                dest_md[match_idx]['files'] = new_md['files']

        return dest_md


def add_parameters(parser: argparse.ArgumentParser) -> None:
    """Adds parameters
    Arguments:
        parser: instance of argparse
    """
    parser.add_argument('--epsg', type=int, nargs='?',
                        help='default epsg code to use if a file doesn\'t have a coordinate system')
    parser.add_argument('--full_plot_fill', action='store_true',
                        help='clipped images will be color filled to match the plot dimensions (outside the original image boundaries)')
    parser.add_argument('--plot_column', type=str,
                        help='the name of the column in the plot geometry file containing plot names')
    parser.add_argument('sensor', type=str, help='the name of the sensor associated with the source files')
    parser.add_argument('plot_file', type=str, help='the path of the GeoJSON file to use for plot boundaries')


def perform_process(transformer: transformer_class.Transformer, check_md: dict, transformer_md: dict, full_md: list) -> dict:
    """Performs the processing of the data
    Arguments:
        transformer: instance of transformer class
        check_md: metadata associated with this request
        transformer_md: metadata associated with this transformer
        full_md: the full set of metadata
    Return:
        Returns a dictionary with the results of processing
    """
    # pylint: disable=unused-argument
    # loop through the available files and clip data into plot-level files
    processed_files = 0
    processed_plots = 0
    start_timestamp = datetime.datetime.now()
    file_list = check_md['list_files']()
    files_to_process = __internal__.get_files_to_process(file_list, transformer.args.sensor, transformer.args.epsg)
    logging.info("Found %s files to process", str(len(files_to_process)))

    container_md = []
    if files_to_process:
        # Get all the possible plots
        logging.debug("Plots file: '%s' column: '%s'", str(transformer.args.plot_file), str(transformer.args.plot_column))
        all_plots = __internal__.load_plotfile(transformer.args.plot_file, transformer.args.plot_column)
        logging.debug("Loaded %s plots", str(len(all_plots)))

        for filename in files_to_process:
            processed_files += 1
            file_path = files_to_process[filename]['path']
            file_bounds = files_to_process[filename]['bounds']
            sensor = files_to_process[filename]['sensor_name']
            logging.debug("File bounds: %s", str(file_bounds))

            overlap_plots = __internal__.find_plots_intersect_boundingbox(file_bounds, all_plots)
            logging.info("Have %s plots intersecting file '%s'", str(len(overlap_plots)), filename)

            file_spatial_ref = file_bounds.GetSpatialReference()
            for plot_name in overlap_plots:
                processed_plots += 1
                plot_bounds = __internal__.convert_geometry(overlap_plots[plot_name], file_spatial_ref)
                logging.debug("Clipping out plot '%s': %s", str(plot_name), str(plot_bounds))
                if __internal__.calculate_overlap_percent(plot_bounds, file_bounds) < 0.10:
                    logging.info("Skipping plot with too small overlap: %s", plot_name)
                    continue
                tuples = __internal__.geometry_to_tuples(plot_bounds)

                plot_md = __internal__.cleanup_request_md(check_md)
                plot_md['plot_name'] = plot_name

                if filename.endswith('.tif'):
                    # If file is a geoTIFF, simply clip it
                    out_path = os.path.join(check_md['working_folder'], plot_name)
                    out_file = os.path.join(out_path, filename)
                    if not os.path.exists(out_path):
                        os.makedirs(out_path)

                    if not transformer.args.full_plot_fill:
                        __internal__.clip_raster_intersection(file_path, file_bounds, plot_bounds, out_file)
                    else:
                        logging.info("Clipping image to plot boundary with fill")
                        __internal__.clip_raster(file_path, tuples, out_path=out_file, compress=True)

                    cur_md = __internal__.prepare_container_md(plot_name, plot_md, sensor, file_path, [out_file])
                    container_md = __internal__.merge_container_md(container_md, cur_md)

                elif filename.endswith('.las'):
                    out_path = os.path.join(check_md['working_folder'], plot_name)
                    out_file = os.path.join(out_path, filename)
                    if not os.path.exists(out_path):
                        os.makedirs(out_path)

                    __internal__.clip_las(file_path, tuples, out_path=out_file)

                    cur_md = __internal__.prepare_container_md(plot_name, plot_md, sensor, file_path, [out_file])
                    container_md = __internal__.merge_container_md(container_md, cur_md)

    return {
        'code': 0,
        'container': container_md,
        configuration.TRANSFORMER_NAME:
        {
            'utc_timestamp': datetime.datetime.utcnow().isoformat(),
            'processing_time': str(datetime.datetime.now() - start_timestamp),
            'total_file_count': len(file_list),
            'processed_file_count': processed_files,
            'total_plots_processed': processed_plots,
            'sensor': transformer.args.sensor
        }
    }
