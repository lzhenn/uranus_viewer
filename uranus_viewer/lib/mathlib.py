#/usr/bin/env python3
import numpy as np
import xarray as xr


def windspeed(var1,var2):
    
    ws_value=np.sqrt(var1.values**2+var2.values**2)
    return ws_value

def locate_position(lat0, lon0, lat_mesh, lon_mesh):
    # Calculate the index of the closest latitude and longitude points
    lat_index = np.abs(lat_mesh[:,0] - lat0).argmin()
    lon_index = np.abs(lon_mesh[0,:] - lon0).argmin()
    latmin, latmax=lat_mesh[:,0].min(), lat_mesh[:,0].max()
    lonmin, lonmax=lon_mesh[0,:].min(), lon_mesh[0,:].max()
    if lat0>latmax or lat0<latmin:
        lat_index=-1
    if lon0>lonmax or lon0<lonmin:
        lon_index=-1
   # Return the matrix index
    return lat_index, lon_index

def search_tc_center(slp, ws10, idx, idy, r=50):
    # slice the xarray variable to select the box of interest
    box = slp[max(0, idx-r):min(slp.shape[0], idx+r+1), 
              max(0, idy-r):min(slp.shape[1], idy+r+1)]
    # find the minimum value and its position in the box
    
    min_pos = np.unravel_index(box.argmin(), box.shape)
    min_idx, min_idy = min_pos
    try:
        minslp=box[min_idx, min_idy].values
    except AttributeError:
        minslp=box[min_idx, min_idy]
    box=ws10[max(0, idx-r):min(ws10.shape[0], idx+r+1),
             max(0, idy-r):min(ws10.shape[1], idy+r+1)]
    max_pos = np.unravel_index(box.argmax(), box.shape)
    try:
        maxws=box[max_pos[0], max_pos[1]].values
    except AttributeError:
        maxws=box[max_pos[0], max_pos[1]]
    # convert the position in the box to the position in the original xarray variable
    min_x = min_idx + max(0, idx-r)
    min_y = min_idy + max(0, idy-r)
    
    return min_x, min_y, minslp, maxws