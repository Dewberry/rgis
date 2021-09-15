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
    def __init__(self, dbname, host, port, user, password):

        self.rdb = rivdb.RiverDatabase(dbname, host, port, user, password)
        self.schema = self.rdb.SCHEMA

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
                logging.info("  {0} OK".format(obj.name))
            except:
                logging.error("  {0} - failure!".format(obj))

    # def rasImportLayersIntoRdbTables(self):
    #     """Import chosen layers into PostGIS database."""
    #     from .dlg_rasImportDataIntoRasTables import DlgImportDataIntoRasTables

    #     self.addInfo("<br><b>Import data into RAS PostGIS tables...</b>")
    #     if not self.curConnName or not self.schema:
    #         self.addInfo("No PostGIS database or schema selected. Choose a connection and schema.")
    #         return
    #     importData = DlgImportDataIntoRasTables(self)
    #     importData.exec_()

    def run(self, flow_area_2d_file, out_ras_geom_file):
        try:
            self.rdb.connect_pg()
            self.rdb.create_schema(self.rdb.SCHEMA)
            self.rdb.create_spatial_index()
            self.rasCreateRdbTables()
            self.rdb.load_gis_file(flow_area_2d_file, heco.FlowAreas2d().name)
            r2d.ras2dCreate2dPoints(self)
            # r2d.ras2dPreviewMesh(self)
            r2d.ras2dSaveMeshPtsToGeometry(self, out_ras_geom_file)
        except Exception as e:
            logging.error(e, exc_info=1)
        finally:
            self.rdb.disconnect_pg()


if __name__ == "__main__":

    rgis = RiverGIS("gis", "localhost", "5432", "docker", "docker")
    rgis.run(
        "/home/abrazeau/workbench/repos/rgis/data/input_area.shp",
        "/home/abrazeau/workbench/repos/rgis/data/test_prj.g01",
    )
