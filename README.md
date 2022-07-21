# ethnic-segregation-rw

This repository contains software to compute the segregation indices
based on Class Coverage Times (CCT), as explained in the paper:  

> S. Sousa, V. Nicosia
> "Quantifying ethnic segregation in cities through random walks".
> arXiv: https://arxiv.org/abs/2010.10462

The code can be copied, used, modified, and redistributed under the
terms of the MIT/Expat License. Please see the file LICENSE for
additional details.


Software requirements:
\* *The list bellow indicates the versions at which the experiments were conducted.*
- anaconda 4.8.2 (recommended for operating system compatibility)
- Python 3.7.6
- numpy 1.18.1
- pandas 1.0.1
- shapely 1.6.4
- fiona 1.8.6


# Obtaining the edge-list and node properties for CCT data  

The script `get_adjacent_polygons.py` returns the edge-list and the nodes
properties (ethnicity distribution). It takes the following inputs:

> shape:     Shapefile of the census tracts at country extent  
> code_name: The loop uk code (GISJOIN) string to match the tracts  
> pop_data:  CSV file containing the tracts within the CSA  
> fileout:   Optional filename to save the nodes properties of the resulting  

Constructing the individual CSV file for each city is trivial and can be
obtained either by filtering the data by `name` using your preferred editor or
using a library such as pandas in python to group rows by the column. The CSV
files from this process are available in the folder
[census-original](census-original/). Note that the largest component is not
automatically provided by the script, a version that provides the largest
connected component is available
[here](https://mygit.katolaz.net/covid_19_ethnicity/rw-segregation/src/branch/master/cct/get_adjacencies_connected.py).

Note that the `code_name` parameter is different for each country, the Shapefile
for US systems can be downloaded at [US Census
Tiger](https://www2.census.gov/geo/tiger/TIGER2010DP1/Tract_2010Census_DP1.zip)
and they use the `GEOID10` for joining data with census tables. The Shapefile for
the UK systems can be downloaded at the [UK Data
Service](https://statistics.ukdataservice.ac.uk/dataset/2011-census-geography-boundaries-wards-and-electoral-divisions)
where the `geo_code` is used to merge data.

US example:  
`python get_adjacent_polygons.py Tract_2010Census_DP1.shp GEOID10 data/census-original/us/SF1P9_Census_Atlanta.csv nodes_ethnics_Census_Atlanta.txt > edge_list_Census_Atlanta.txt`

UK example:  
`python get_adjacent_polygons.py infuse_ward_lyr_2011_merged.shp geo_code data/census-original/uk/QS211_Wards_Bristol.csv nodes_ethnics_Wards_Bristol.txt > edge_list_Wards_Bristol.txt`

A second edge-list file will be created automatically with following the naming
convention:
`edges_ids_` + `pop_data filename` = `edges_ids_QS211_Wards_Bristol.txt`

Edge-list details:
* Node IDs are defined by the row index of the population table
* Each edge is reported once


# Running the CCT random walk  

The random walk on the adjacency graph of real systems can be simulated by
running the python script `rw_cct_fnt.py` and the following input must be
provided:

> edge:    Edge-list file;  
> prop:    Node properties file with the frequency of each ethnicity;  
> num:     Number of walk repetitions from the same node;  
> idx:     range of nodes or single node ID to run the walk from, e.g.: 1, 0-10.  

Example:  
`python rw_cct_fnt.py data/census-clean/edge_list_Wards_Bristol.txt data/census-clean/nodes_ethnics_distrib_Wards_Bristol.txt 1000 0-143 >> [output_filename]`

All the files provided in [data/census-clean](data/census-clean) are already in
the input format needed for the script. The script will print, for each node
passed to the <idx> parameter, the following results:
```
"Node ID" "[list with the CCT for each fraction c]"
```
where each line corresponds to the average over all repetitions of the random
walk from node i. The runtime depends on the number of nodes in the graph and
the number of repetitions per node, which is expected to take a few seconds
for a single node with 100 repetitions. The runtime can vary significantly per
repetition since the trajectory of the walker is a pure stochastic process.

The synthetic systems use the script `rw_cct_fnt_synthetic.py` which needs the
following input:

> edge:    Edge-list file;  
> C:       number of classes to construct the MxC matrix;  
> num:     Number of walk repetitions from the same node;  
> idx:     range of nodes or single node ID to run the walk from, e.g.: 1, 0-10.  
> mode:    the spatial pattern of the population distribution

Example:  
`python rw_cct_fnt_synthetic.py data/synthetic/edge_list_16x16_grid.txt 5 1000 0-256 544 >> [output_filename]`

*For details about each pattern, please consult the inline documentation in the script.*


# Description of data sets

**Census-clean:**

- `edge_list_Census_[city].txt`: Edge-list for US systems with the assigned numeric ID.
- `edge_list_Wards_[city].txt`: Edge-list for UK systems with the assigned numeric ID.
- `nodes_ethnics_Census_[city].txt`: Ethnic data associated with the census tract.
- `nodes_ethnics_distrib_Wards_[city].txt`: Ethnic data associated with the wards.

**Census-original:**
- `QS211_[scale]_[city].csv`: Census ethnic data with headers at the scales
LSOA, OA and Wards.
- `SF1P9_Census_[city].csv`: Census ethnic data with headers.

**dfa:**  
- `edges_un_[scale]_[city].txt`: Edge-list for UK systems with the assigned
numeric ID \*.
- `nodes_[scale]_entropy_[city].txt`: Entropy of the population distribution for
UK systems with the assigned numeric ID \*.

\* *All files available at three spatial scales LSOA, OA and Wards.*

**Socio-economic-US:**
- `ACSDP5Y2011.DP03_data_with_overlays_2021-02-21T040944.csv`: Socio-Economic Indicators at `CSA` level.
- `ACSDP5Y2011.DP03_metadata_2021-02-21T040944.csv`: Metadata information.
- `ACSDP5Y2011.DP03_table_title_2021-02-21T040944.txt`: Description of the fields in the `CSV` table.

**Synthetic:**
- `edge_list_[size]_[boundary].txt`: Edge-list for synthetic systems for
distinct lattice size and boundary conditions.
