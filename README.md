Risk-based Ranking of Highway Bridges
==========================

Setup 1: GIS data to PostGIS database and CSV data
-----
GIS data in shapefile can be downloaded from OpenStreetMap (https://www.openstreetmap.org)

TAZ file for LA can be downloaded from http://gisdata.scag.ca.gov/SitePages/GIS%20Library.aspx

CTPP ODs can be downloaded from http://www.fhwa.dot.gov/planning/census_issues/ctpp

1. Installation of QGIS and PostGIS can be found on https://github.com/megacell/user-equilibrium-synthetic

2. Instructions for importing shapefiles to PostGIS database and export to CSV data can be found on https://github.com/megacell/user-equilibrium-synthetic


Setup 2: create nodes and links files (db tables and csv files)
-----
1. Import highways (using filter) and taz shapefiles (and clip rectangular if
   applicable) into QGIS

2. Clip highway with either taz shapefile or clip rectangular
   (QGIS->Vector->Geoprocessing->Clip)

3. Add point layer __nodes and line layer __motorways

4. Create nodes and motorways (single line). Open Snap for better selection (QGIS->Snapping
   Options)

5. 
