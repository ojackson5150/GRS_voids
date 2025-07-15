"""
@author: Josepf Amador

Determines wether a set of points falls inside a polygon defined by another
set of points. In the context of this research, determines if a celestial
object falls in the void catalog footprint"""

from descartes import PolygonPatch # comes with bug. See NOTE in perimeter.py and step 1 in Part 2 of README.
from shapely.geometry import Point
from shapely.geometry import Polygon, MultiPolygon
import pandas as pd
from matplotlib import pyplot as plt
from shapely.ops import unary_union
import numpy as np
import alphashape



def filter_by_footprint(cel_obj, foot_print_fn, voids, ret_mask = False):
    #edited to filter by the MultiPolygons because multiple regions
    
    """Filter the celestial object table by the polygon defined by the points in 
     the foot_print table.
     input: cel_obj: pandas.DataFrame of celestial objects
     Returns: Pandas dataframe filtered by the footprint
     optional: Return the mask that determines wether or not it is in the footprint (Does not filter it), 
     otherwise, just returns the pandas DataFrame"""
    
    # Read in data
    # cel_obj = pd.read_excel(cel_obj_fn)
    footprint = pd.read_excel(foot_print_fn)

    
    # for SUTTER void catalog catalog
    # convert footprint data to format readable by polygon
    
    footprint_list = []
    for x, y in zip(footprint.RAdeg, footprint.DEdeg):
        footprint_list.append((x,y))
        
    footprint_polygon = Polygon(footprint_list)

    # Represent the celestial objects as Point objects
    cel_obj_Points = []
    for ra, de in zip(cel_obj.RAdeg, cel_obj.DEdeg):
        cel_obj_Points.append(Point(ra, de))

    # Do the roar (filtering)
    in_foot_print = [None]*len(cel_obj_Points)
    for i, point in enumerate(cel_obj_Points):
        in_foot_print[i] = footprint_polygon.contains(point)

    if ret_mask:
        cel_obj['foot_print_mask'] = in_foot_print
        return cel_obj
    else:
        cel_obj = cel_obj[in_foot_print]
        fig, ax = plt.subplots()
        ax.scatter(voids.RAdeg, voids.DEdeg, marker='.')
        ax.scatter(footprint.RAdeg, footprint.DEdeg, marker='.')
        if isinstance(footprint_polygon, Polygon):
            ax.add_patch(PolygonPatch(footprint_polygon, alpha=0.5))
        elif isinstance(footprint_polygon, MultiPolygon):
            for polygon in footprint_polygon.geoms:
                #print(polygon)
                ax.add_patch(PolygonPatch(polygon, alpha=0.5))
                #print('success')
        ax.scatter(cel_obj.RAdeg, cel_obj.DEdeg, marker='+', s=250, color='red')
        fig.legend(['Void Centers','Bounding Points', 'Footprint','Gamma-Ray Sources'])
        ax.set_ylabel('Declination (Deg)')
        ax.set_xlabel('Right Ascencion (Deg)')
        ax.set_title('Cosmic Void Centers and Perimeter')

        plt.show()
        return cel_obj
