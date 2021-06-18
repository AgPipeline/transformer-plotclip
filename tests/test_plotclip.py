#!/usr/bin/env python3
"""Tests plotclip.py
"""
import os
import re
import json
from subprocess import getstatusoutput
import pytest

# The name of the source file to test and it's path
SOURCE_FILE = 'plotclip.py'
SOURCE_PATH = os.path.abspath(os.path.join('.', SOURCE_FILE))

# Path relative to the current folder where the testing JSON file are
TESTING_JSON_FILE_PATH = './tests'

def test_exists():
    """Asserts that the source file is available"""
    assert os.path.isfile(SOURCE_PATH)


def test_usage():
    """Program prints a "usage" statement when requested"""
    for flag in ['-h', '--help']:
        ret_val, out = getstatusoutput(f'{SOURCE_PATH} {flag}')
        assert re.match('usage', out, re.IGNORECASE)
        assert ret_val == 0


def test_fail_get_sr_from_crs():
    """Passes bad GeoJSON to the function to test failure conditions"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'fail_get_sr_from_crs.json'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test_data = json.load(in_file)
        for test in test_data:
            with pytest.raises(RuntimeError):
                pc.__internal__.get_sr_from_crs(test)


def test_get_sr_from_crs():
    """Checks retrieving the various CRS"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'get_sr_from_crs.json'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test_data = json.load(in_file)
        for test in test_data:
            res = pc.__internal__.get_sr_from_crs(test)
            assert res is not None


def test_fail_get_geojson_file_sr():
    """Tests getting an incorrect GeoJSON result from a file or elsewhere"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'fail_get_geojson_file_sr.json'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test_data = json.load(in_file)
        for test in test_data:
            with pytest.raises(RuntimeError):
                pc.__internal__.get_geojson_file_sr(test)


def test_get_geojson_file_sr():
    """Tests getting an a CRS from a simulated image file GeoJson"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'get_geojson_file_sr.json'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test_data = json.load(in_file)
        for test in test_data:
            res = pc.__internal__.get_geojson_file_sr(test)
            assert res is not None


def test_fail_get_plot_key_name():
    """Tests failure to get a plot name"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'fail_get_plot_key_name.json'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test_data = json.load(in_file)
        for test in test_data:
            res = pc.__internal__.get_plot_key_name(test)
            assert res == (None, None)


def test_get_plot_key_name():
    """Tests getting a plot key"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'get_plot_key_name.json'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test_data = json.load(in_file)
        for test in test_data:
            res = pc.__internal__.get_plot_key_name(test)
            assert res != (None, None)
            assert len(res) == 2
            assert res[0] in test
            assert res[1] == test[res[0]]


def test_fail_load_plot_file():
    """Tests the failure cases for loading JSON files"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    # This test builds off the other tests that parse the JSON.
    # We don't try all combinations of failure, just the ones directly in the function
    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'fail_load_plot_file.txt'))
    assert os.path.exists(data_file_name)

    with open(data_file_name, 'r') as in_file:
        test = in_file.readline()
        test_index = 1
        while test:
            json_file_name = '%d_fail_load_plot.json' % test_index
            with open(json_file_name, 'w') as out_file:
                out_file.write(test)
            with pytest.raises(RuntimeError):
                pc.__internal__.load_plot_file(json_file_name)
            os.unlink(json_file_name)
            test_index += 1
            test = in_file.readline()


def test_load_plot_file():
    """Tests the successful loading of a JSON file"""
    # pylint: disable=import-outside-toplevel
    import plotclip as pc

    # We just test the successful loading of a valid file
    data_file_name = os.path.realpath(os.path.join(TESTING_JSON_FILE_PATH, 'plots.geojson'))
    assert os.path.exists(data_file_name)

    # Testing with different plot ID columns
    for column_name in [None, 'Id', 'purple']:
        plots = pc.__internal__.load_plot_file(data_file_name, column_name)
        assert len(plots) == 2
        plot_keys = plots.keys()
        assert 1 in plot_keys
        assert 'Plot 2' in plot_keys
        for one_key in plot_keys:
            assert plots[one_key] is not None
