Risk-based Ranking of Highway Bridges
==========================

Setup 1: GIS data to PostGIS database and CSV data
-----
GIS data in shapefile can be downloaded from OpenStreetMap (https://www.openstreetmap.org)

TAZ file for LA can be downloaded from http://gisdata.scag.ca.gov/SitePages/GIS%20Library.aspx

CTPP ODs can be downloaded from http://www.fhwa.dot.gov/planning/census_issues/ctpp

Bridge data can be downloaded from http://www.fhwa.dot.gov/bridge/nbi/ascii.cfm

1. Installation of QGIS and PostGIS can be found on https://github.com/megacell/user-equilibrium-synthetic

2. Instructions for importing shapefiles to PostGIS database and export to CSV data can be found on https://github.com/megacell/user-equilibrium-synthetic


Setup 2: create nodes and links files (db tables and csv files)
-----
1. Import highways (using filter) and taz shapefiles (and clip rectangular if
   applicable) into QGIS

2. Clip highway with either taz shapefile or clip rectangular
   (QGIS->Vector->Geoprocessing->Clip)

3. Add point layer *nodes* and line layer *motorways*

4. Create nodes and motorways (single line). Open Snap for better selection (QGIS->Snapping
   Options)
  1. nodes should have x and y attributes as longitude and latitude respectively

5. Import nodes and motorways to PostGIS

6. Run motorway2link.sql to create a table of *links*

7. Create a layer of voronoi polygons of nodes (QGIS->Geometry tools->Voronoi
   polygons), import the shapefile to PostGIS

8. Import bridge data into QGIS, two numerical attributes corresponding to
   longitudes and latitudes of bridges should be firstly converted from text data and added to the database

9. Clip bridge data with taz or clipping rectangular and filter all the bridges
   with the following criteria, save as shapefile and import the shapefile into
   PostGIS as *bridges*

   ```sql
   -- recond_type_005a = 1 for bridges ON highways
   -- route_prefix_005b = 1 or 2 for interstate and numbered highways
   -- super- or sub-structure condition rating should be equal to less than 6
   "recond_type_005a" = '1' and
   "route_prefix_005b" in ('1', '2') and
   ("superstructure_cond_059" in ('0', '1', '2', '3', '4', '5', '6') or
   "substructure_cond_060" in ('0', '1', '2', '3', '4', '5', '6'))
   ```

10. Update table bridges with *onlink* data  by running bridgeonlink.sql

Setup 3: update python files with directory paths
-----
Use ```bash grep -nr -exclude-dir ./.git 'xxx' ./``` to make sure all directory
is correct


Setup 4: run python scripts
-----
Run **ue_LA_county.py** for traffic assignment without bridges

Run **bridge_ranking_par.py** for Monte Carlo simulation. Make sure that the
process number and sample number have been correctly set
