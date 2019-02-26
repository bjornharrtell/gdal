#!/usr/bin/env pytest
# -*- coding: utf-8 -*-
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  FlatGeobuf driver test suite.
# Author:   Björn Harrtell <bjorn@wololo.org>
#
###############################################################################
# Copyright (c) 2018-2019, Björn Harrtell <bjorn@wololo.org>
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
###############################################################################

import json
import math
import os


from osgeo import osr
from osgeo import ogr
from osgeo import gdal

import gdaltest
import ogrtest
import pytest

### utils

def verify_flatgeobuf_copy(name, fids, names):

    if gdaltest.fgbpoint_feat is None:
        print('Missing features collection')
        return False

    fname = os.path.join('tmp', name + '.fgb')
    ds = ogr.Open(fname)
    if ds is None:
        print('Can not open \'' + fname + '\'')
        return False

    lyr = ds.GetLayer(0)
    if lyr is None:
        print('Missing layer')
        return False

    ######################################################
    # Test attributes
    ret = ogrtest.check_features_against_list(lyr, 'FID', fids)
    if ret != 1:
        print('Wrong values in \'FID\' field')
        return False

    lyr.ResetReading()
    ret = ogrtest.check_features_against_list(lyr, 'NAME', names)
    if ret != 1:
        print('Wrong values in \'NAME\' field')
        return False

    ######################################################
    # Test geometries
    lyr.ResetReading()
    for i in range(len(gdaltest.fgbpoint_feat)):

        orig_feat = gdaltest.fgbpoint_feat[i]
        feat = lyr.GetNextFeature()

        if feat is None:
            print('Failed trying to read feature')
            return False

        if ogrtest.check_feature_geometry(feat, orig_feat.GetGeometryRef(),
                                          max_error=0.001) != 0:
            print('Geometry test failed')
            gdaltest.fgbpoint_feat = None
            return False

    gdaltest.fgbpoint_feat = None

    lyr = None

    return True


def copy_shape_to_flatgeobuf(name, compress=None):
    if gdaltest.flatgeobuf_drv is None:
        return False

    if compress is not None:
        if compress[0:5] == '/vsig':
            dst_name = os.path.join('/vsigzip/', 'tmp', name + '.fgb' + '.gz')
        elif compress[0:4] == '/vsiz':
            dst_name = os.path.join('/vsizip/', 'tmp', name + '.fgb' + '.zip')
        elif compress == '/vsistdout/':
            dst_name = compress
        else:
            return False
    else:
        dst_name = os.path.join('tmp', name + '.fgb')

    print dst_name

    ds = gdaltest.flatgeobuf_drv.CreateDataSource(dst_name)
    if ds is None:
        return False

    ######################################################
    # Create layer
    lyr = ds.CreateLayer(name, None, ogr.wkbPoint)
    if lyr is None:
        return False

    ######################################################
    # Setup schema (all test shapefiles use common schmea)
    ogrtest.quick_create_layer_def(lyr,
                                   [('FID', ogr.OFTReal),
                                    ('NAME', ogr.OFTString)])

    ######################################################
    # Copy in gjpoint.shp

    dst_feat = ogr.Feature(feature_def=lyr.GetLayerDefn())

    src_name = os.path.join('data', name + '.shp')
    shp_ds = ogr.Open(src_name)
    shp_lyr = shp_ds.GetLayer(0)

    feat = shp_lyr.GetNextFeature()
    gdaltest.fgbpoint_feat = []

    while feat is not None:
        gdaltest.fgbpoint_feat.append(feat)

        dst_feat.SetFrom(feat)
        lyr.CreateFeature(dst_feat)

        feat = shp_lyr.GetNextFeature()

    shp_lyr = None
    lyr = None

    ds = None

    return True

### tests

def test_ogr_flatgeobuf_1():

    gdaltest.flatgeobuf_drv = ogr.GetDriverByName('FlatGeobuf')

    if gdaltest.flatgeobuf_drv is not None:
        return
    pytest.fail()

def test_ogr_flatgeobuf_9():
    if gdaltest.flatgeobuf_drv is None:
        pytest.skip()

    gdaltest.tests = [
        ['gjpoint', [1], ['Point 1']],
        #['gjline', [1], ['Line 1']],
        #['gjpoly', [1], ['Polygon 1']],
        #['fgbmultipoint', [1], ['MultiPoint 1']],
        #['fgbmultiline', [2], ['MultiLine 1']],
        #['fgbmultipoly', [2], ['MultiPoly 1']]
    ]

    for i in range(len(gdaltest.tests)):
        test = gdaltest.tests[i]

        rc = copy_shape_to_flatgeobuf(test[0])
        assert rc, ('Failed making copy of ' + test[0] + '.shp')

        rc = verify_flatgeobuf_copy(test[0], test[1], test[2])
        assert rc, ('Verification of copy of ' + test[0] + '.shp failed')
