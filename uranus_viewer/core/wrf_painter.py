#!/usr/bin/env python3
'''
Base Class: Painter
    Derived Class: WRFPainter
'''

from . import painter
from ..lib import utils, io, const, cfgparser, mathlib

from netCDF4 import Dataset
import cmaps
import numpy as np
import pandas as pd
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import shapely.geometry as sgeom
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
matplotlib.use('Agg') 

import sys, os, subprocess, glob
import wrf

# -------Global Envrionmental Vars--------
print_prefix='lib.painter>>'

CWD=sys.path[0]

# Constants
BIGFONT=22
MIDFONT=18
SMFONT=14
# Province level
SHP_LV1='cnhimap.dbf'
# City level
SHP_LV2='gadm41_CHN_2.dbf'
# County level
SHP_LV3='county_2004.dbf'
# For D03 Coastline
SHP_LV4='china_coastline.dbf'
# -------Global Envrionmental Vars--------


class WRFPainter(painter.Painter):
    """WRF Painter"""
    def update(self):
        '''update WRF specific files'''
        wrfcfg=self.cfg['WRF']
        self.domid=wrfcfg['domid']
        self.utc_flag=wrfcfg.getboolean('wrf_flag')
        self.max_workers=int(wrfcfg['max_workers'])
        
        self.spatial_draw=wrfcfg.getboolean('spatial_draw')
        self.debug=wrfcfg.getboolean('debug')
        self.met_em_flag=wrfcfg.getboolean('met_em_flag')
        self.varlist=cfgparser.cfg_get_varlist(self.cfg, 'WRF', 'varlist')
        self.form_anim=wrfcfg.getboolean('form_animation')
        
        self.run_sp=wrfcfg.getboolean('run_special')
        if self.run_sp:
            self.spjob=wrfcfg['special_job'].split(':')
    def load_data(self):
        '''load WRF files according to domain ID'''
        if self.met_em_flag:
            keywd='met_em'
        else:
            keywd='wrfout'
        utils.write_log(f'get {keywd}...')
        fn_stream=subprocess.check_output(
                'ls '+self.arch_root+'/'+keywd+'?'+self.domid+'*', 
                shell=True).decode('utf-8')
        fn_list=fn_stream.split()
        self.file_list=fn_list
        if self.debug:
            self.file_list=self.file_list[:self.max_workers]
        wrf_sample=Dataset(fn_list[0])
        lsmask = wrf.getvar(wrf_sample, 'LANDMASK')
        # Get the latitude and longitude points
        self.lats, self.lons = wrf.latlon_coords(lsmask)
        self.nrow, self.ncol = (
            self.lats.sizes['south_north'],
            self.lats.sizes['west_east'])
        utils.write_log(f'fecthed %3d {keywd} files' % len(fn_list))
    def drive(self):
        import concurrent.futures
        self.update()
        self.load_data()
        fn_lst=self.file_list
        if self.spatial_draw:
            if self.cfg.has_option('WRF','varplv'):
                self.plv=int(self.cfg['WRF']['varplv'])
            else:
                self.plv=0
            # overlap tc track
            tc=self.cfg['WRF']['overlap_tc'].split(':')
            para_dic={}
            if len(tc)==2:
                obv_trck=io.get_ibtrack(tc[0], tc[1], resample_frq='1H')
                sim_trck=pd.read_csv(
                    f'{self.arch_root}/tc_track.csv', parse_dates=['time'],index_col=['time'])
                para_dic.update({'obv_trck':obv_trck, 'sim_trck':sim_trck})
            for varname in self.varlist:
                para_dic.update({
                    'fig_root':self.fig_root, 
                    'lats':self.lats.values, 
                    'lons':self.lons.values,'add_shp':self.add_shp, 
                    'plv':self.plv,
                    'file_list':self.file_list,'met_em_flag':self.met_em_flag,
                    'varname':varname, 'diff_flag':self.diff_flag})
                if self.diff_flag:
                    para_dic['file_list']=self.file_list2
                with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                    futures = [executor.submit(
                        draw2d_map,  idx, file_path, para_dic) for idx, file_path in enumerate(fn_lst)] 
                        
                    concurrent.futures.wait(futures)
                    for idx, f in enumerate(futures):
                        if f.exception(): 
                            print(f.result(),idx)
                            #exit()
            # end of loop vars
        # end of spatial draw
        
        if self.form_anim:
            for varname in self.varlist:
                if self.met_em_flag:
                    prefix='met_em'
                else:
                    prefix='wrf'
                if self.diff_flag:
                    fnpattern=f'{prefix}.{varname}-diff.'
                else:
                    fnpattern=f'{prefix}.{varname}.'
                anim_fn=f'{fnpattern}mp4'
                if self.cfg.has_option('anim','fnpattern'):
                    fnpattern=self.cfg['anim']['fnpattern']
                if self.cfg.has_option('anim','outfn'):
                    anim_fn=self.cfg['anim']['outfn']
                fps=self.cfg['anim']['fps']
                utils.write_log(f'form animation@${fps}fps:{fnpattern}.mp4')
                '''
                This command takes a sequence of PNG image files with a frame rate of {fps} frames per second, 
                converts them to an MP4 video file with a frame rate of 8 frames per second, 
                and applies the libx264 codec and yuv420p pixel format to the output.
                '''
                spell=f'cd {self.fig_root};ffmpeg -y -r {fps} -pattern_type glob -i \'{fnpattern}*.png\''+\
                    f' -c:v libx264 -vf "fps={fps},format=yuv420p" {anim_fn}'
                utils.write_log(spell)
                subprocess.call(spell, shell=True)
        # end of form animation
        
        if self.run_sp:
            if self.spjob[0]=='loc_track':
                obv_trck=io.get_ibtrack(self.spjob[1], self.spjob[2], resample_frq='1H')
                sim_trck = pd.DataFrame(columns=['time', 'lat', 'lon', 'wsmax', 'slp'])
                para_dic={
                    'lats':self.lats.values, 
                    'lons':self.lons.values, 
                    'obv_trck':obv_trck, 'met_em_flag':self.met_em_flag}
            
                with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                    futures = [executor.submit(
                        find_tc_center, idx, file_path, para_dic) for idx, file_path in enumerate(fn_lst)] 
                        
                    concurrent.futures.wait(futures)
                    for idx, f in enumerate(futures):
                        if f.exception(): 
                            print(f.result(),idx)
                            exit()
                        if f.result() is not None:
                            sim_trck=pd.concat(
                                [sim_trck, pd.DataFrame(f.result(), index=[0])], ignore_index=True)
                sim_trck.to_csv(f'{self.arch_root}/tc_track.csv')
def find_tc_center(idx, fn, para_dic):
    '''find tc center'''
    lats,lons=para_dic['lats'],para_dic['lons']
    obv_trck=para_dic['obv_trck']
    met_em_flag=para_dic['met_em_flag']
    if met_em_flag:
        tstmp=fn.split('/')[-1][11:-3]
    else:
        tstmp=fn.split('/')[-1][11:]
    tfrm = datetime.strptime(tstmp,'%Y-%m-%d_%H:%M:%S')
    utils.write_log(f'try to find tc center @ {tstmp}')
    try:
        match_row = obv_trck.loc[tfrm]
        obv_lat, obv_lon=float(match_row['HKO_LAT']),float(match_row['HKO_LON'])
        idx,idy=mathlib.locate_position(obv_lat, obv_lon, lats, lons)
        if -1 in [idx,idy]: 
            utils.write_log(f'TC not in domain @ {tstmp}')
            return
    except KeyError:
        utils.write_log(f'no match row found for {tstmp}')
        return
    utils.write_log(f'TC in domain @ {tstmp}, deal with slp...')
    if met_em_flag:
        slp=io.get_single_var2d(fn, 'PMSL', source='wrf')
        u_srf=io.get_single_var2d(fn, 'UU', source='wrf')[:,1:]
        v_srf=io.get_single_var2d(fn, 'VV', source='wrf')[1:,:]
        ws10=mathlib.windspeed(u_srf, v_srf)
    else:
        slp=io.get_single_var2d(fn, 'slp', source='wrf')
        ws10=io.get_single_var2d(fn, 'wspd_wdir10', source='wrf')
    sim_idx, sim_idy, minslp, maxws=mathlib.search_tc_center(slp, ws10, idx, idy)
    rt_dic={'time':tfrm, 'lat':lats[sim_idx, sim_idy], 'lon':lons[sim_idx, sim_idy], 'wsmax':maxws, 'slp':minslp}
    return rt_dic 
def draw2d_map(idx, fn, para_dic):
    '''
    draw 2d spatial 
    '''
    # read var
    add_shp=para_dic['add_shp'] 
    varname, lats, lons=para_dic['varname'],para_dic['lats'], para_dic['lons']
    met_em_flag=para_dic['met_em_flag']
    utils.write_log(f'{print_prefix}ploting {varname} from {fn}')
    try:
        att_dic=const.WRF_VAR_DIC[varname]    
    except KeyError:
        att_dic=const.WRF_VAR_DIC['precip']
    if met_em_flag:
        tstmp=fn.split('/')[-1][11:-3]
    else:
        tstmp=fn.split('/')[-1][11:]
    tfrm = datetime.strptime(tstmp,'%Y-%m-%d_%H:%M:%S') 
    if 'obv_trck' in para_dic:
        obv_trck,sim_trck=para_dic['obv_trck'],para_dic['sim_trck']
        strt_loc=obv_trck.index.get_loc(sim_trck.index[0])
        try:
            curr_loc=obv_trck.index.get_loc(tfrm)
        except KeyError:
            exit() 
    if varname in ['ua', 'va', 'geopt']:
        var2d=io.get_var2d_from_3d(fn, varname, source='wrf',plv=para_dic['plv'])
    elif varname == 'uv10':
        if met_em_flag:
            label='ReAnalysis'
            var2d=io.get_single_var2d(fn, 'LANDMASK', source='wrf')
            u_srf=io.get_single_var2d(fn, 'UU', source='wrf')[:,1:]
            v_srf=io.get_single_var2d(fn, 'VV', source='wrf')[1:,:]
            var2d.data=mathlib.windspeed(u_srf, v_srf)
        else:
            label='WRF'
            var2d=io.get_single_var2d(
                fn, 'wspd_wdir10', source='wrf', lay=0) 
        
    elif 'precip' in varname:
        diff_frm=int(varname.split('_')[-1])
        if idx<diff_frm:
            return
        else:
            rainc0 = io.get_single_var2d(fn, 'RAINC', source='wrf')
            rainnc0 = io.get_single_var2d(fn, 'RAINNC', source='wrf')
            bkrc0=io.get_single_var2d(fn, 'I_RAINC', source='wrf')
            bkrnc0=io.get_single_var2d(fn, 'I_RAINNC', source='wrf')
            fn_ref=para_dic['file_list'][idx-diff_frm]
            rainc_hrs = io.get_single_var2d(fn_ref, 'RAINC', source='wrf')
            rainnc_hrs = io.get_single_var2d(fn_ref, 'RAINNC', source='wrf')
            bkrc_hrs=io.get_single_var2d(fn_ref, 'I_RAINC', source='wrf')
            bkrnc_hrs=io.get_single_var2d(fn_ref, 'I_RAINNC', source='wrf')
            var2d=rainc0
            var2d.data = (wrf.to_np(rainc0) + wrf.to_np(rainnc0) + 100*wrf.to_np(bkrc0) + 100*wrf.to_np(bkrnc0) - wrf.to_np(rainc_hrs)
                - wrf.to_np(rainnc_hrs)-100*wrf.to_np(bkrc_hrs)-100*wrf.to_np(bkrnc_hrs))
    else:
        var2d=io.get_single_var2d(fn, varname, source='wrf')
    # prepare painting
    '''
    tfrm = pd.to_datetime(
        np.datetime_as_string(var2d['Time'],unit='h'),
        format='%Y-%m-%dT%H')
    '''
    cmap=att_dic['cmap']
    levels=att_dic['clvl']
    
    if para_dic['plv'] == 0:
        title_txt=f'{varname} ({att_dic["units"]})'
    else:
        title_txt=f'{varname} ({att_dic["units"]}) at {para_dic["plv"]}hPa'
    if add_shp:
        fig = plt.figure(figsize=[8, 8],layout='constrained')
        proj=wrf.get_cartopy(var2d)
        # Set projection and plot the main figure
        ax = fig.add_axes([0.1, 0.05, 0.9, 0.85], projection=proj)
        # add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                color='grey', linestyle='--', linewidth=1)
        ax.coastlines()
    else:
        fig, ax = plt.subplots(figsize=[10, 8])
    if para_dic['diff_flag']:
        var2d_ref=io.get_single_var2d(
            para_dic['file_list2'][idx], varname, source='wrf')
        var2d.values=var2d.values-var2d_ref.values
        cmap=att_dic['diff_cmap']
        levels=att_dic['diff_clvl']
        title_txt=varname+'_Diff'
        arr=var2d.values
        # Mask the nan values in the array
        plt.pcolormesh(
            lons, lats, arr,
            vmax=levels[-1], vmin=levels[0],
            cmap=cmap)
    else:
        if add_shp:
            try:
                if att_dic['norm']==True:
                    norm=matplotlib.colors.BoundaryNorm(
                        boundaries=levels,ncolors=cmap.N,extend='both')
                    plt.pcolormesh(
                        lons, lats, var2d.values,
                        norm=norm,
                        transform=ccrs.PlateCarree(), cmap=cmap)
            except KeyError:
                plt.pcolormesh(
                    lons, lats, var2d.values,
                    vmax=levels[-1], vmin=levels[0],
                    transform=ccrs.PlateCarree(), cmap=cmap)
        else:
            plt.pcolormesh(
                lons, lats, var2d.values,
                vmax=levels[-1], vmin=levels[0],
                cmap=cmap)
        if 'obv_trck' in para_dic:
            if add_shp:
                ax.plot(obv_trck[strt_loc:curr_loc]['HKO_LON'], obv_trck[strt_loc:curr_loc]['HKO_LAT'], 
                        color='black', transform=ccrs.PlateCarree(), 
                        marker='^', linewidth=2, linestyle='solid', markersize=4,label='OBV')
            else:
                ax.plot(obv_trck[strt_loc:curr_loc]['HKO_LON'], obv_trck[strt_loc:curr_loc]['HKO_LAT'], color='black', 
                    marker='^', linewidth=2, linestyle='solid', markersize=4,label='OBV')
            try:
                curr_loc_sim=sim_trck.index.get_loc(tfrm)
            except:
                return
            if add_shp:
                ax.plot(sim_trck[:curr_loc_sim]['lon'], sim_trck[:curr_loc_sim]['lat'], color='blue', 
                    marker='*', linewidth=1, linestyle='dashed', markersize=3,label=label,transform=ccrs.PlateCarree())
            else:
                ax.plot(sim_trck[:curr_loc_sim]['lon'], sim_trck[:curr_loc_sim]['lat'], color='blue', 
                    marker='*', linewidth=1, linestyle='dashed', markersize=3,label=label)
            ax.legend(loc='upper left')
        '''
        plt.contourf(
            lons, lats, 
            var2d.values,
            levels=levels, extend='both', 
            cmap=cmap)
        '''
    tfrm_str=tfrm.strftime('%Y-%m-%d %H:%M:%S')
    title_txt=title_txt+'@'+tfrm_str
    
    plt.title(title_txt)

    # Add a color bar
    cbar=plt.colorbar(ax=ax)
    
    # Adjust the position of the colorbar
    cbar.ax.set_position([0.9, 0.1, 0.05, 0.7])
    #plt.tight_layout() 
    savefig(
        para_dic, varname, tfrm.strftime('%Y%m%d%H%M'))
    plt.close()

def savefig(para_dic, varname, time_str):
    '''save fig portal'''
    save_dir=para_dic['fig_root']
    diff_flag=para_dic['diff_flag']
    met_flag=para_dic['met_em_flag']
    plv=para_dic['plv']
    
    if met_flag:
        prefix='met_em'
    else:
        prefix='wrf'
    if plv ==0:
        if diff_flag:
            fn=f'{prefix}.{varname}-diff.{time_str}.png'
        else:
            fn=f'{prefix}.{varname}.{time_str}.png'
    else:
        if diff_flag:
            fn=f'{prefix}.{varname}-diff.{str(plv)}hpa.{time_str}.png'
        else:
            fn=f'{prefix}.{varname}.{str(plv)}hpa.{time_str}.png'
 
    
    utils.write_log(f'saving {fn} to {save_dir}') 
    if not(os.path.isdir(save_dir)):
            os.mkdir(save_dir)        
    plt.savefig(os.path.join(save_dir,fn), dpi=100)
