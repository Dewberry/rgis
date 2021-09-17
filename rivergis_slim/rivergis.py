# -*- coding: utf-8 -*-

"""
/***************************************************************************
Name                 : RiverGIS
Description          : HEC-RAS tools for QGIS
Date                 : December, 2015
copyright            : (C) 2015 by RiverGIS Group
email                : rpasiok@gmail.com, damnback333@gmail.com
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from __future__ import absolute_import

# from builtins import range
# import os
# import json
import logging

# from . import river_database as rivdb
# from . import hecobjects as heco
# from . import ras2dFunctions as r2d

import river_database as rivdb
import hecobjects as heco
import ras2dFunctions as r2d


class RiverGIS:
    def __init__(self, dbname, host, port, user, password, schema, srid):

        self.rdb = rivdb.RiverDatabase(dbname, host, port, user, password, schema, srid)

    def rasCreateRdbTables(self):

        tables = [
            heco.BreakLines2d,
            heco.BreakPoints2d,
            heco.FlowAreas2d,
            heco.MeshPoints2d,
        ]
        tables.sort(key=lambda x: x().order)

        for hecobj_class in tables:
            obj = self.rdb.process_hecobject(hecobj_class, "pg_create_table")
            try:
                logging.info("{0} OK".format(obj.name))
            except:
                logging.error("{0} - failure!".format(obj))

    def run(
        self,
        out_ras_geom_file,
        flow_area_2d_file,
        breaklines_file=None,
        mesh_preview_file=None,
        mesh_pts_preview_file=None,
    ):
        try:

            # SETUP
            self.rdb.connect_pg()
            self.rdb.create_schema(self.rdb.SCHEMA)
            self.rdb.create_spatial_index()
            self.rasCreateRdbTables()

            # LOAD GIS DATA IN
            self.rdb.load_gis_file(flow_area_2d_file, heco.FlowAreas2d().name)
            if breaklines_file:
                self.rdb.load_gis_file(breaklines_file, heco.BreakLines2d().name)

            # CREATE RAS DATA
            r2d.ras2dCreate2dPoints(self)

            r2d.ras2dPreviewMesh(self, output_filename=mesh_preview_file)

            if mesh_pts_preview_file:
                r2d.ras2dPreviewMeshPoints(self, mesh_pts_preview_file)

            r2d.ras2dSaveMeshPtsToGeometry(self, out_ras_geom_file)

        except Exception as e:
            logging.error(e, exc_info=1)
        finally:
            # CLEANUP
            self.rdb.drop_schema(self.rdb.SCHEMA, cascade=True)
            self.rdb.disconnect_pg()


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s -- [ %(filename)s:%(funcName)s:%(lineno)d ]  ::  %(message)s"
    )

    rgis = RiverGIS("gis", "localhost", "5432", "docker", "docker", "testschema", 2284)
    rgis.run(
        "/home/abrazeau/workbench/repos/rgis/data/test_prj.g01",
        "/home/abrazeau/workbench/repos/rgis/data/input_area.shp",
        "/home/abrazeau/workbench/repos/rgis/data/breaklines.geojson",
        "/home/abrazeau/workbench/repos/rgis/data/mesh_preview.geojson",
        "/home/abrazeau/workbench/repos/rgis/data/mesh_pts_preview.geojson",
    )
