# Krillmod
## Overview of files
krillmod.py: Main python file for working with other modules.\
get_trajectory.py: Module containing functions for pre-processing data.\
plot_trajectory.py: Module with plotting functionality.

## SSMU polygons:

[comment]: <> (This is a comment, it will not be included)
[comment]: <> (The shapefiles from ccamlr)
```python
# Make sure the .shx file is in the same folder as .shp
import geopandas as gpd
file="ssmusPolygon.shp" 
df = gpd.read_file(file)
print(df)
```
```
 GAR_ID  ...                                           geometry
0    95561  ...  POLYGON ((-50.00000 -64.00000, -50.00000 -65.0...
1    95562  ...  POLYGON ((-61.11353 -64.32835, -61.11216 -64.3...
2    95563  ...  POLYGON ((-60.53607 -61.73011, -60.53600 -61.7...
3    95564  ...  POLYGON ((-56.37107 -61.54253, -56.38300 -61.5...
4    95588  ...  POLYGON ((-58.24770 -63.45529, -58.25221 -63.4...
5    95589  ...  POLYGON ((-55.50100 -62.20903, -55.50140 -62.2...
6    95590  ...  POLYGON ((-55.50100 -62.20903, -56.37107 -61.5...
7    95591  ...  POLYGON ((-30.00000 -64.00000, -30.95238 -64.0...
8    95592  ...  POLYGON ((-44.76102 -59.89825, -45.96612 -60.5...
9    95593  ...  POLYGON ((-43.00918 -60.91173, -44.39301 -60.7...
10   95594  ...  POLYGON ((-43.00918 -60.91173, -43.03899 -60.9...
11   95595  ...  POLYGON ((-30.00000 -57.00000, -30.95238 -57.0...
12   95596  ...  POLYGON ((-37.12200 -52.83920, -37.12200 -52.8...
13   95597  ...  POLYGON ((-37.13500 -55.84870, -37.13500 -55.8...
14   95598  ...  POLYGON ((-30.00000 -64.00000, -30.00000 -57.0...
15   95599  ...  POLYGON ((-26.64942 -60.17018, -26.85006 -60.1...
16   95604  ...  POLYGON ((-57.03671 -63.35067, -54.19000 -62.6...
[17 rows x 15 columns]
```
