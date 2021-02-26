import networkx as nx
import matplotlib.pyplot as plt
import osmnx as ox
import geopandas as gpd
import pandas as pd
import shapefile
from shapely.geometry import Polygon,MultiPolygon
from shapely.ops import unary_union
from shapely.geometry import box
from shapely.validation import explain_validity
from shapely.geos import TopologicalError

''' This file recreates the road polygons of an area. Starting by some given road polygons for some area, we download the graph (nodes and edges).
Then based on the intersection nodes we recreate the shape of those polygons in order to model the intersections separetly.
We dissolve the initial road polygons to a Multipolygon, then we create a buffer around each node and we find the intersection of this buffer with the dissolved polygon
this creates new polygons.

2 things need to change for the different datasets. 1) The EPSG line 54 (and 127) +  2) The distance of the buffers around the intersection nodes!'''

#read a polygons file, create a new field same for all entries and dissolve everything to a big Multipolygon
#save the dissolved Multipoligon to a new shp file
p = gpd.read_file("polygons_to_dissolve.shp")
p['d_field'] = 1
p_new = p.dissolve(by='d_field')
p_new.to_file("diss_polygons.shp")

#Before generate a graph with osmnx we need to project shapefile to wgs84 lat lon!! Later we ll project the generated nodes to the CRS of the specific shapefile
nn = gpd.read_file("diss_polygons.shp")
print(nn.crs)
nn = nn.to_crs("EPSG:4326")
nn.to_file("wgs84_dissolved_polygons.shp")

d_polygons = gpd.read_file("wgs84_dissolved_polygons.shp")
pl = d_polygons.loc[0,'geometry']
print(pl.bounds)
b = box(pl.bounds[0],pl.bounds[1],pl.bounds[2],pl.bounds[3])
print(b)

#Generate the graph WGS84
G = ox.graph_from_polygon(b, network_type='drive')
ox.plot_graph(G)
ox.io.save_graph_shapefile(G,"graph")

nodes_shapefile = "C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/git2/graph/nodes.shp"
nod = gpd.read_file(nodes_shapefile)

deg_higher_3 = nod[(nod.street_cou > 2)]
deg3 = deg_higher_3.to_file("int_nodes.shp")

#open the intersection nodes file and create a buffer for each node
#int_nodes is projected at WGS84 --> Change to correspondig crs
df_multipolygon = gpd.read_file("diss_polygons.shp")
df_nodes = gpd.read_file("int_nodes.shp")
df_nodes = df_nodes.to_crs("EPSG:4326")
df_nodes.to_file("int_nodes_projected.shp")
df_final_nodes = gpd.read_file("int_nodes_projected.shp")
print(df_final_nodes.crs)
print(df_multipolygon.crs)
all_polygons = df_multipolygon.loc[0,'geometry']
v = all_polygons.is_valid
print(v)
if v == False:
    print(explain_validity(all_polygons))
    all_polygons = all_polygons.buffer(0)
    print(explain_validity(all_polygons))

#first create a small buffer for each node in order to find the nodes that are close together (in the same polygon or really close)
buffers_small = []
for i in range(len(df_final_nodes)-1):
    geomm = df_final_nodes.loc[i,'geometry']
    bf = geomm.buffer(10)
    buffers_small.append([bf,geomm])
#find which nodes are intersecting with the initial dissolved polygon and keep only them
bufs = []
only_nodes2 = []
for b in buffers_small:
    if b[0].intersects(all_polygons):
        bufs.append(b)
        only_nodes2.append(b[1])

for i in bufs:
    print(b[0])
    print(b[1])

print("check1")
#based on if the small buffers are intersecting find the nodes that are close together (remove them from only_nodes2 list)
only_buffers = []
only_nodes = []
for bb in bufs:
    for bbb in bufs:
        if bb[0]!=bbb[0] and bb[0].intersects(bbb[0]):
            only_nodes.append(bb[1])
            only_buffers.append(bb[0])
            if bb[1] in only_nodes2:
                only_nodes2.remove(bb[1])
#make a union of buffers of the nodes that are close together, and find the centroid of this union buffer
op = []
for i in only_nodes:
    print(i)
    test = []
    s = i.buffer(52)
    print(s)
    for b in only_buffers:
        if s.intersects(b):
            test.append(b)
    cu = unary_union(test)
    op.append(cu)

print("check2")

opop = []
for mult in op:
    cent = mult.centroid
    opop.append(cent)

#all_pnts is a list which contains all the unique (isolate) nodes from the only_nodes2 list (the nodes that are not too close with other nodes) +
#+ the centroids of the union buffer of the nodes that are close together
all_pnts = []
for i in only_nodes2:
    all_pnts.append(i)

for ii in opop:
    all_pnts.append(ii)


#int_polygons a list that contains all the intersection polygons (based on the intersection of the buffer with the initial dissolved layer of road polygons)
int_polygons = []
for i in all_pnts:
    bfr = i.buffer(27)
    check = bfr.intersection(all_polygons)
    c = check.simplify(0.05, preserve_topology=False)
    int_polygons.append(c)

#dissolve the separate intersection polygons to one Multipolygon(use it later to find the symmetric difference with the dissolced intial polygon)
cu2 = unary_union(int_polygons)
#create a file with all the intersection polygons
d = {'geometry':[]}
for p in cu2:
    d['geometry'].append(p)

gdf = gpd.GeoDataFrame(d, crs="EPSG:3879")
gdf.to_file("intersection_polygons.shp")

'''
int_polygons2 = []
for i in cu2:
    print(i)
    s= i.simplify(0.05, preserve_topology=False)
    ss = s.symmetric_difference(all_polygons)
    print("A")
    print(ss)
    #int_polygons2.append(s)'''

rest_polygons_multi = cu2.symmetric_difference(all_polygons)
rest_polygons = list(rest_polygons_multi)



new_road_polygons = []
for i in int_polygons:
    new_road_polygons.append(i)

for ii in rest_polygons:
    new_road_polygons.append(ii)

#for f in new_road_polygons:
    #print(f)
