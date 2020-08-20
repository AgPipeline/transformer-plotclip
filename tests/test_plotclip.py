#!/usr/bin/env python3
"""Tests plotclip.py
"""
import os
import re
from subprocess import getstatusoutput
import pytest

# The name of the source file to test and it's path
SOURCE_FILE = 'plotclip.py'
SOURCE_PATH = os.path.abspath(os.path.join('.', SOURCE_FILE))


def test_exists():
    """Asserts that the source file is available"""
    assert os.path.isfile(SOURCE_PATH)


def test_usage():
    """Program prints a "usage" statement when requested"""
    for flag in ['-h', '--help']:
        ret_val, out = getstatusoutput(f'{SOURCE_PATH} {flag}')
        assert ret_val == 0
        assert re.match('usage', out, re.IGNORECASE)


def test_fail_get_sr_from_crs():
    """Passes bad GeoJSON to the function to test failure conditions"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'nothing': None, 'nada': 42})

    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'type': None})

    # Invalid or missing keys for CRS and EPSG
    for type_name in ['name', 'EPSG']:
        with pytest.raises(RuntimeError):
            pc.__internal__.get_sr_from_crs({'properties': None, 'type': type_name})
        with pytest.raises(RuntimeError):
            pc.__internal__.get_sr_from_crs({'properties': {}, 'type': type_name})
        with pytest.raises(RuntimeError):
            pc.__internal__.get_sr_from_crs({'properties': {'silly': 'putty'}, 'type': type_name})

    # Incomplete CRS tests
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'name': None}, 'type': 'name'})
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'name': 'urn:ogc:def:crs'}, 'type': 'name'})
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'name': 'urn:ogc:def:crs:'}, 'type': 'name'})
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'name': 'urn:ogc:def:crs:EPSG'}, 'type': 'name'})

    # Incomplete EPSG tests
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'code': None}, 'type': 'EPSG'})
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'name': 'invalid code'}, 'type': 'EPSG'})

    # Invalid type check
    with pytest.raises(RuntimeError):
        pc.__internal__.get_sr_from_crs({'properties': {'name': 'invalid code'}, 'type': 'invalid'})


def test_get_sr_from_crs():
    """Checks retrieving the various CRS"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    # Check normal CRS
    res = pc.__internal__.get_sr_from_crs({"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::32612"}})
    assert res is not None

    # Check normal EPSG
    res = pc.__internal__.get_sr_from_crs({"type": "EPSG", "properties": {"code": "4326"}})
    assert res is not None

    # Check special GeoJSON CRS
    res = pc.__internal__.get_sr_from_crs({"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC::CRS84"}})
    assert res is not None


def test_fail_get_geojson_file_sr():
    """Tests getting an incorrect GeoJSON result from a file or elsewhere"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    with pytest.raises(RuntimeError):
        pc.__internal__.get_geojson_file_sr({'crs': {'properties': {'name': 'invalid code'}, 'type': 'invalid'}})


def test_get_geojson_file_sr():
    """Tests getting an a CRS from a simulated image file GeoJson"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    res = pc.__internal__.get_geojson_file_sr({'crs':
                                                   {'properties': {'name': 'urn:ogc:def:crs:EPSG::32612'},
                                                    'type': 'name'}})
    assert res is not None

    # The default GeoJSON CRS should be returned
    res = pc.__internal__.get_geojson_file_sr({})
    assert res is not None


def test_fail_get_plot_key_name():
    """Tests failure to get a plot name"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    # Variations for not having a valid key
    res = pc.__internal__.get_plot_key_name({})
    assert res is None

    res = pc.__internal__.get_plot_key_name({'nothing': 'here'})
    assert res is None

    # Doesn't check beyond the first level of the dict
    res = pc.__internal__.get_plot_key_name({'nothing': 'here', 'something': {'id': 'here'}})
    assert res is None


def test_get_plot_key_name():
    """Tests getting a plot key"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    res = pc.__internal__.get_plot_key_name({'observationUnitName': 'the_id'})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'observationUnitName'
    assert res[1] == 'the_id'

    # Check combinations of "plot" and "name" with "id" as well as case insensitivity
    res = pc.__internal__.get_plot_key_name({'plot_id': 'the_id'})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'plot_id'
    assert res[1] == 'the_id'

    res = pc.__internal__.get_plot_key_name({'nameid': 'the_id'})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'nameid'
    assert res[1] == 'the_id'

    res = pc.__internal__.get_plot_key_name({'Name_Id': 'the_id'})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'Name_Id'
    assert res[1] == 'the_id'

    # Check for plain "id"
    res = pc.__internal__.get_plot_key_name({'ID': 'the_id'})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'ID'
    assert res[1] == 'the_id'

    # Check for embedded "id" (last chance before failing)
    res = pc.__internal__.get_plot_key_name({'planetoID': 'the_id'})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'planetoID'
    assert res[1] == 'the_id'

    # Check for a returned None
    res = pc.__internal__.get_plot_key_name({'observationUnitName': None})
    assert res is not None
    assert len(res) == 2
    assert res[0] == 'observationUnitName'
    assert res[1] is None
