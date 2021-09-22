"""
This program is free software. You can use, copy, distribute, modify and/or redistribute
it under the terms of the MIT/Expat License. See the file LICENSE for more details.

(c) 2021 Sandro Sousa

If you use this software please consider citing the original paper:

S. Sousa, V. Nicosia
"Quantifying ethnic segregation in cities through random walks"


-----------------------------------------------------------------------------

Reads a shapefile and population data, find the neighbouring nodes in the
giant component and save results to local files.
Inputs:
<shape>:    Shapefile of the census tracts at country extent;
<code_name>:The loop uk code (GISJOIN) string to match the tracts;
<pop_data>: CSV file containing the tracts within the CSA;
[fileout]:  Optional filename to save the nodes properties of the resulting
connected component.

Returns: a undirected edge list in the format (i, j) and
the equivalent prop_file format for the random walk. Considers pop_data IDs as
input for polygon filtering.

-----------------------------------------------------------------------------
"""

import pandas as pd
from shapely.geometry import shape
import fiona
import sys
import itertools
import os


if len(sys.argv) < 4:
    print("Usage: %s <shape> <code_name> <pop_data> [fileout]\n" % sys.argv[0])
    exit(1)

arg_shape = sys.argv[1]
arg_code = str(sys.argv[2])
pop_data = pd.read_csv(sys.argv[3], sep=',', dtype=str)
edges = []

# create edge_list using shape file features ID
with fiona.open(arg_shape) as shp:
    feat_dict = {feat['properties'][arg_code]:feat for feat in shp}
shp.close()

for i, j in itertools.combinations(pop_data.iloc[:, 0], 2):
        geom1 = shape(feat_dict[i]['geometry'])
        geom2 = shape(feat_dict[j]['geometry'])
        if geom1.touches(geom2):
            idx_i = pop_data[pop_data.iloc[:,0] == i].index[0]
            idx_j = pop_data[pop_data.iloc[:,0] == j].index[0]
            print(i, j)
            edges.append("%s %s\n" %(idx_i,idx_j))
            # print(idx_j, idx_i)
            # for validation remove comment
            # print(idx_i, idx_j, i, j)
            # print(idx_j, idx_i, j, i)

# save prop file if argument is passed
# ignore area code and sum column
if len(sys.argv) > 4:
    pop_data.iloc[:, 2:].to_csv(sys.argv[4], sep=' ', header=False, index=True)

# save edge list with
edge_ids = "edges_ids_%s.txt" %os.path.basename(sys.argv[3])[:-4]
with open(edge_ids, "w") as out:
    for line in edges:
        out.write(line)
