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
from builtins import range
import os
import logging
from datetime import datetime

# from . import hecobjects as heco
import hecobjects as heco
from math import floor
from psycopg2 import ProgrammingError
import geopandas as gpd
from shapely.ops import voronoi_diagram


def ras2dCreate2dPoints(rgis):
    """
    Create 2D computational points for each 2D flow area.
    Points are regularly spaced (based on CellSize attribute of the FlowArea2D table) except for breaklines, where they are aligned to form a cell face exactly at a breakline.
    Points spacing along and across a breakline is read from CellSizeAlong and CellSizeAcross attributes of BreakLines2D table, respectively. A number of cells rows to align with a beakline can be given.
    Create breakpoints at locations where a cell face is needed (on a breakline).
    """
    logging.info("Creating computational points for 2D flow areas")

    # and create breaklines with a linear measure
    qry = 'SELECT * FROM "{0}"."{1}"'.format(rgis.rdb.SCHEMA, heco.FlowAreas2d().name)
    chk2dAreas = rgis.rdb.run_query(qry, fetch=True)
    if not chk2dAreas:
        logging.info(
            "No 2d flow area in the database. Import or create it before generating 2d computational points. Cancelling..."
        )
        return

    logging.info("Creating regular mesh points...")

    # create regular mesh points
    # and delete points located too close to the 2D area boundary
    rgis.rdb.process_hecobject(heco.MeshPoints2d, "pg_create_table")
    rgis.rdb.process_hecobject(heco.MeshPoints2d, "pg_create_mesh")

    # find which breakline line belongs to which 2d flow area
    # and create breaklines with a linear measure
    rgis.rdb.process_hecobject(heco.BreakLines2d, "pg_flow_to_breakline")
    rgis.rdb.process_hecobject(heco.BreakLines2d, "pg_breaklines_m")
    rgis.rdb.process_hecobject(heco.BreakLines2d, "pg_drop_by_buffer")

    logging.info("Creating mesh points along structures...")

    # check if breaklines and breakpoints exist in the database
    bls_exist = False
    bps_exist = False
    for t in rgis.rdb.list_tables():
        if t == heco.BreakLines2d().name:
            bls_exist = True
        if t == heco.BreakPoints2d().name:
            bps_exist = True

    if bls_exist:
        # find measures of breakpoints along breaklines
        # there was a change in the alg name between PostGIS 2.0 and 2.1
        # ST_Line_Locate_Point -> ST_LineLocatePoint
        # qry = "SELECT PostGIS_Full_Version() AS ver;"
        # postgisVersion = rgis.rdb.run_query(qry, True)[0]["ver"].split('"')[1][:5]
        # pgMajV = int(postgisVersion[:1])
        # pgMinV = int(postgisVersion[2:3])
        # if pgMajV < 2:
        #     locate = "ST_Line_Locate_Point"
        # elif pgMajV >= 2 and pgMinV == 0:
        #     locate = "ST_Line_Locate_Point"
        # else:
        #     locate = "ST_LineLocatePoint"

        # find breakline that a breakpoint is located on ( tolerance = 10 [map units] )
        if bps_exist:
            breakPtsLocTol = 10
            rgis.rdb.process_hecobject(
                heco.BreakPoints2d, "pg_bpoints_along_blines", tolerance=breakPtsLocTol, func_name="ST_LineLocatePoint"
            )
        # find breaklines with measures
        qry = """
        SELECT
            "BLmID",
            "AreaID",
            "CellSizeAlong" AS csx,
            "CellSizeAcross" AS csy,
            ST_Length(geom) AS len,
            "RowsAligned" AS rows
        FROM
            "{0}"."breaklines2d_m";
        """
        qry = qry.format(rgis.rdb.SCHEMA)
        bls = rgis.rdb.run_query(qry, True)

        for line in bls:
            if not line["csx"] or not line["csy"] or not line["rows"]:
                logging.info("Empty BreakLines2d attribute! Cancelling...")
                return
            dist_x = float(line["csx"])
            width = float(line["csy"])
            id = line["BLmID"]
            leng = float(line["len"])
            rows = int(line["rows"])
            imax = int(leng / (dist_x))

            # check if breakpoints exist on the breakline
            qry = """
            SELECT
                bp."BPID"
            FROM
                "{0}"."breaklines2d_m" AS bl,
                "{0}"."breakpoints2d" AS bp
            WHERE
                bl."BLmID" = {1} AND
                bp."BLmID" = bl."BLmID";
            """
            try:
                qry = qry.format(rgis.rdb.SCHEMA, id)
                bp_on_bl = rgis.rdb.run_query(qry, True)
                logging.debug("Breakline BLmID={0}, {1}".format(id, bp_on_bl))
            except ProgrammingError:
                bp_on_bl = None

            if not bp_on_bl:
                # no BreakPoints2d: create aligned mesh at regular interval = CellSizeAlong
                logging.debug("Creating regular points for breakline BLmID={0} (no breakpoints)".format(id))
                for i in range(0, imax + 1):
                    dist = i * dist_x
                    for j in range(0, rows):
                        rgis.rdb.process_hecobject(
                            heco.MeshPoints2d,
                            "pg_aligned_mesh",
                            cellsize=dist_x,
                            measure=dist,
                            offset=j * width + width / 2,
                            blid=id,
                        )

            # create cell faces at breakline's breakpoints
            else:
                qry = """
                SELECT DISTINCT
                    p."Fraction"
                FROM
                    "{0}"."breakpoints2d" AS p
                WHERE
                    p."BLmID" = {1};
                """
                qry = qry.format(rgis.rdb.SCHEMA, id)
                ms = rgis.rdb.run_query(qry, True)

                logging.debug("Creating breakpoints for structure id={0} (with breakpoints)".format(id))
                sm_param = 4
                db_min = 10.0 ** 9
                # breakpoints m list (linear locations on current structure)
                bm = []
                # linear measures of mesh points to be created
                mpts = []

                for m in ms:
                    bm.append(float(m["Fraction"]))
                    logging.debug("BreakPoint2d fraction: {0}".format(float(m["Fraction"])))

                # sort the list
                bm.sort()

                for i, m in enumerate(bm):
                    # calculate minimal distance between breakpoints
                    if i > 0:
                        db_min = min(bm[i] - bm[i - 1], db_min)
                logging.debug("Min dist between breakpoints db_min={0}".format(db_min))
                # create 2 mesh points on both sides of a breakpoint at a distance db_min / sm_param
                dist_min = min(db_min / sm_param, 0.5 * dist_x / leng)
                cs_min = dist_min * leng
                logging.debug("dist_min={0}, cs_min={1}".format(dist_min, cs_min))
                for m in bm:
                    mpts.append(max(0.0001, m - dist_min))
                    mpts.append(min(m + dist_min, 0.9999))

                # find gaps between points along a breakline longer than 3 * dist_min
                gaps = []
                for i, m in enumerate(mpts):
                    logging.debug("m={0}".format(m))
                    if i > 0:
                        dist = m - mpts[i - 1]
                        if dist > 3 * dist_min:
                            gaps.append([m, dist])

                # create mesh points filling the gaps
                for g in gaps:
                    m, dist = g
                    # how many points to insert?
                    k = int(floor(dist / (2 * dist_min)))
                    # distance between new points
                    cs = dist / k
                    for j in range(1, k):
                        mpts.append(m - j * cs)
                        logging.debug("gap: dist={0}, m={1}".format(cs, m - j * cs))

                # insert aligned mesh points into table
                for m in sorted(mpts):
                    for j in range(0, rows):
                        rgis.rdb.process_hecobject(
                            heco.MeshPoints2d,
                            "pg_aligned_mesh",
                            cellsize=cs_min,
                            measure=m * leng,
                            offset=j * width + width / 2,
                            blid=id,
                        )

    logging.info("Deleting mesh points located too close to each other or outside the 2D area...")
    rgis.rdb.process_hecobject(heco.MeshPoints2d, "pg_clean_points")
    logging.info("Done")


def ras2dPreviewMesh(rgis, output_filename):
    """Build and load Voronoi polygons for the mesh points"""

    logging.info("Creating mesh_preview")

    areas = rgis.rdb.table_to_gdf("flowareas2d")
    mesh_pts = rgis.rdb.table_to_gdf("meshpoints2d")
    assert areas.crs == mesh_pts.crs, "flowareas2d and meshpoints2d have different coordinate systems!"

    multipts = mesh_pts.geometry.unary_union
    voronoi = voronoi_diagram(multipts).geoms

    new_feats = []
    for area in areas.geometry:
        for item in voronoi:
            poly_cut = item.intersection(area)
            new_feats.append(poly_cut)

    data = list(zip(range(len(new_feats)), new_feats))

    out_mesh_preview = gpd.GeoDataFrame(
        gpd.pd.DataFrame(data, columns=["id", "geometry"]), crs=areas.crs, geometry="geometry"
    )

    too_many_faces = out_mesh_preview[
        out_mesh_preview.geometry.apply(lambda x: True if len(x.exterior.coords) > 8 else False)
    ]

    for i in too_many_faces["id"]:
        logging.warning("Mesh Cell {0} has too many faces!".format(i))

    if len(too_many_faces) > 0:
        for g in too_many_faces["geometry"]:
            rgis.rdb.subdivide_poly_update_meshpts2d(g.wkt)
        logging.warning("Recursing to eliminate cells with too many faces...")
        ras2dPreviewMesh(rgis, output_filename)

    else:
        if output_filename:
            logging.info("Saving mesh_preview to: {}".format(output_filename))
            out_mesh_preview.to_file(output_filename, driver="GeoJSON")
        logging.info("Number of faces has been enforced")

    return


def ras2dPreviewMeshPoints(rgis, output_filename):
    logging.info("Saving meshpoints2d to: {}".format(output_filename))
    rgis.rdb.table_to_gdf("meshpoints2d").to_file(output_filename, driver="GeoJSON")
    return


def ras2dSaveMeshPtsToGeometry(rgis, geoFileName):
    """Saves mesh points from current schema and table 'mesh_pts' to HEC-RAS geometry file"""

    logging.info("Saving 2D Flow Area to HEC-RAS geometry file {}...".format(geoFileName))

    # get mesh points extent
    qry = """
    SELECT
        ST_XMin(ST_Collect(geom)) AS xmin,
        ST_XMax(ST_collect(geom)) AS xmax,
        ST_YMin(ST_collect(geom)) AS ymin,
        ST_YMax(ST_collect(geom)) AS ymax
    FROM
        "{0}"."{1}";
    """
    qry = qry.format(rgis.rdb.SCHEMA, heco.MeshPoints2d().name)
    pExt = rgis.rdb.run_query(qry, True)[0]
    xmin, xmax, ymin, ymax = [pExt["xmin"], pExt["xmax"], pExt["ymin"], pExt["ymax"]]
    buf = max(0.2 * (xmax - xmin), 0.2 * (ymax - ymin))
    pExtStr = "{:.2f}, {:.2f}, {:.2f}, {:.2f}".format(xmin - buf, xmax + buf, ymax + buf, ymin - buf)

    # get list of mesh areas
    qry = """
    SELECT
        "AreaID",
        "Name",
        ST_X(ST_Centroid(geom)) AS x,
        ST_Y(ST_Centroid(geom)) AS y,
        ST_NPoints(geom) AS ptsnr
    FROM
        "{0}"."{1}";
    """
    qry = qry.format(rgis.rdb.SCHEMA, heco.FlowAreas2d().name)
    t = ""

    for area in rgis.rdb.run_query(qry, True):
        qry = """
        SELECT
            ST_AsText(geom) AS geom
        FROM
            "{0}"."{1}"
        WHERE
            "AreaID" = {2};
        """
        qry = qry.format(rgis.rdb.SCHEMA, heco.FlowAreas2d().name, area["AreaID"])
        res = rgis.rdb.run_query(qry, True)[0]["geom"]
        ptsList = res[9:-2].split(",")
        ptsTxt = ""
        for pt in ptsList:
            x, y = [float(c) for c in pt.split(" ")]
            ptsTxt += "{:>16.4f}{:>16.4f}\n".format(x, y)
        t += """

Storage Area={0:<14},{1:14},{2:14}
Storage Area Surface Line= {3:d}
{4}
Storage Area Type= 0
Storage Area Area=
Storage Area Min Elev=
Storage Area Is2D=-1
Storage Area Point Generation Data=,,,
""".format(
            area["Name"], area["x"], area["y"], area["ptsnr"], ptsTxt
        )

        qry = """
        SELECT
            ST_X(geom) AS x,
            ST_Y(geom) AS y
        FROM
            "{0}"."{1}"
        WHERE
            "AreaID" = {2};
        """
        qry = qry.format(rgis.rdb.SCHEMA, heco.MeshPoints2d().name, area["AreaID"])
        pkty = rgis.rdb.run_query(qry, True)

        coords = ""
        for i, pt in enumerate(pkty):
            if i % 2 == 0:
                coords += "{:16.2f}{:16.2f}".format(float(pt["x"]), float(pt["y"]))
            else:
                coords += "{:16.2f}{:16.2f}\n".format(float(pt["x"]), float(pt["y"]))

        t += """Storage Area 2D Points= {0}
{1}
Storage Area 2D PointsPerimeterTime={2}
Storage Area Mannings=0.06
2D Cell Volume Filter Tolerance=0.003
2D Face Profile Filter Tolerance=0.003
2D Face Area Elevation Profile Filter Tolerance=0.003
2D Face Area Elevation Conveyance Ratio=0.02

""".format(
            len(pkty), coords, datetime.now().strftime("%d%b%Y %H:%M:%S")
        )

    qry = """
    SELECT 
        'Breakline' || LPAD("BLID"::TEXT, 3, '0') AS breakline_name,
        ST_NPoints(geom) AS num_coords,
        SUBSTRING(ST_AsText(geom), '\((.+)\)') AS coords
    FROM {0}.breaklines2d;
    """.format(
        rgis.rdb.SCHEMA
    )
    for i, breakline in enumerate(rgis.rdb.run_query(qry, True)):
        ptsList = breakline["coords"].split(",")
        ptsTxt = ""
        for i, pt in enumerate(ptsList):
            x, y = [float(c) for c in pt.split(" ")]
            if i % 2 == 0:
                ptsTxt += "{:>16.4f}{:>16.4f}".format(x, y)
            else:
                ptsTxt += "{:>16.4f}{:>16.4f}\n".format(x, y)

        t += """BreakLine Name={0}
BreakLine CellSize Min=
BreakLine CellSize Max=
BreakLine Near Repeats=3
BreakLine Protection Radius=-1
BreakLine Polyline= {1} 
{2}""".format(
            breakline["breakline_name"], breakline["num_coords"], ptsTxt
        )

    qry = """
    SELECT
        'BoundLine' || LPAD("AreaID"::TEXT, 3, '0')  AS bc_line_name,
        "Name" AS bc_line_storage_area_name,
        ST_X(ST_LineInterpolatePoint(
            ST_ExteriorRing(ST_ConvexHull(ST_Buffer(ST_Union(geom), MAX("CellSize") * 2))),
            0.50)) AS middle_x,
        ST_Y(ST_LineInterpolatePoint(
            ST_ExteriorRing(ST_ConvexHull(ST_Buffer(ST_Union(geom), MAX("CellSize") * 2))),
            0.50)) AS middle_Y,
        ST_NPoints(ST_Union(geom)) AS num_coords,
        SUBSTRING(ST_AsText(ST_ExteriorRing(ST_ConvexHull(ST_Buffer(ST_Union(geom), MAX("CellSize") * 2)))), '\((.+)\)') AS coords
    FROM {0}.flowareas2d
    GROUP BY "AreaID", "Name";""".format(
        rgis.rdb.SCHEMA
    )
    for i, bc_line in enumerate(rgis.rdb.run_query(qry, True)):
        ptsList = bc_line["coords"].split(",")
        ptsTxt = ""
        for i, pt in enumerate(ptsList):
            x, y = [float(c) for c in pt.split(" ")]
            if i == 0:
                start_x, start_y = x, y
            if i == len(ptsList) - 1:
                end_x, end_y = x, y
            if i % 2 == 0:
                ptsTxt += "{:>16.4f}{:>16.4f}".format(x, y)
            else:
                ptsTxt += "{:>16.4f}{:>16.4f}\n".format(x, y)

        t += """BC Line Name={0}                     
BC Line Storage Area={1}          
BC Line Start Position= {2} , {3} 
BC Line Middle Position= {4} , {5} 
BC Line End Position= {6} , {7}
BC Line Arc= {8} 
{9}
BC Line Text Position= 1.79769313486232E+308 , 1.79769313486232E+308
""".format(
            bc_line["bc_line_name"],
            bc_line["bc_line_storage_area_name"],
            start_x,
            start_y,
            bc_line["middle_x"],
            bc_line["middle_y"],
            end_x,
            end_y,
            bc_line["num_coords"],
            ptsTxt,
        )

    if not os.path.isfile(geoFileName):
        createNewGeometry(geoFileName, pExtStr)

    geoFile = open(geoFileName, "r")
    geoLines = geoFile.readlines()
    geoFile.close()

    geoFile = open(geoFileName, "w")
    geo = ""
    for line in geoLines:
        if not line.startswith("Chan Stop Cuts"):
            geo += line
        else:
            geo += t
            geo += line
    geoFile.write(geo)
    geoFile.close()

    logging.info("Saved to: {}".format(geoFileName))


def createNewGeometry(filename, extent):
    t = """Geom Title=Import from RiverGIS
Program Version=5.00
Viewing Rectangle= {0}


Chan Stop Cuts=-1



Use User Specified Reach Order=0
GIS Ratio Cuts To Invert=-1
GIS Limit At Bridges=0
Composite Channel Slope=5
""".format(
        extent
    )
    geoFile = open(filename, "w")
    geoFile.write(t)
    geoFile.close()
