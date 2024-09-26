#!/usr/bin/env python3
"""specific module for IO"""
# ---imports---
import xarray as xr
import wrf
import pandas as pd
from . import const
from netCDF4 import Dataset
# ---Module regime consts and variables---
print_prefix='lib.io>>'

def get_ibtrack(tcname,yr,resample_frq='',src='HKO'):
    ibtracs = pd.read_csv(
        const.machine_dic['IBTRCK_FN'], header=0, skiprows=[1],na_values={"NaN":' '})
    df = ibtracs.loc[(ibtracs.SEASON==int(yr))]
    df = df.loc[(df.TRACK_TYPE=='main')]
    df = df.loc[(df.NAME==tcname)]
    try:
        df = df[['ISO_TIME',f'{src}_LAT',f'{src}_LON',f'{src}_CAT', f'{src}_WIND', f'{src}_PRES']]
    except KeyError:
        df = df[['ISO_TIME',f'{src}_LAT',f'{src}_LON', f'{src}_WIND', f'{src}_PRES']]
        df0 = df[~df.apply(lambda row: row.str.isspace()).any(axis=1)]
        if df0.empty:
            df = df[['ISO_TIME',f'{src}_LAT',f'{src}_LON', f'{src}_PRES']]
            df0= df[~df.apply(lambda row: row.str.isspace()).any(axis=1)]
        if df0.empty:
            df = df[['ISO_TIME',f'{src}_LAT',f'{src}_LON', f'{src}_WIND']]
            df0= df[~df.apply(lambda row: row.str.isspace()).any(axis=1)]
    df =df[~df.apply(lambda row: row.str.isspace()).any(axis=1)]
    
    df['ISO_TIME'] = pd.to_datetime(df['ISO_TIME'])
    for itm in [f'{src}_LAT',f'{src}_LON', f'{src}_WIND', f'{src}_PRES']:
        df[itm]=pd.to_numeric(df[itm])

    df.set_index('ISO_TIME', inplace=True)
    # drop the rows in df that have empty cells
    if not(resample_frq==''):
        df = df.resample(resample_frq).ffill()
    return df 
# ---Classes and Functions---
def get_file_tf(fn, source='roms'):
    '''get file timeframe'''
    if source=='roms':
        ds=xr.open_dataset(fn)
        tf=ds['ocean_time']
        ds.close()
    return tf
def get_single_var2d(fn, varname, source='roms',lay=0):
    '''get single var 2D single frame
        fn: roms file
        varname: variable name
        ifrm: frame index
    '''
    if source=='roms':
        ds=xr.open_dataset(fn)
        var=ds[varname]
        ds.close()
    elif source=='wrf':
        wrf_file=Dataset(fn)
        var = wrf.getvar(wrf_file, varname, timeidx=0)
    if len(var.shape)==3:
        var=var[lay,:,:]
    return var

def get_var2d_from_3d(fn, varname, source='wrf',plv=500):
    '''get single var 2D single frame
        fn: roms file
        varname: variable name
        ifrm: frame index
    '''
    if source=='roms':
        ds=xr.open_dataset(fn)
        var=ds[varname]
        ds.close()
    elif source=='wrf':
        wrf_file=Dataset(fn)
        var = wrf.getvar(wrf_file, varname, timeidx=0)
        # Extract the pressure levels
        p = wrf.getvar(wrf_file, "pressure", timeidx=0)
        var = wrf.interplevel(var, p, plv)
    return var



# ---Unit test---
if __name__ == '__main__':
    pass