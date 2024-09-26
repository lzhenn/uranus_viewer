#!/usr/bin/env python3
'''
Base Class: Painter
    Derived Class: ROMSPainter
'''
from . import painter
from . import driver
from ..lib import utils, io, const, cfgparser, mathlib, station

import cmaps
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg') 

import sys, os, subprocess


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
SHP_LV2='gadm41_CHN_1.dbf'
# County level
SHP_LV3='county_2004.dbf'
# For D03 Coastline
SHP_LV4='china_coastline.dbf'

SHP_SELECT=SHP_LV4
# -------Global Envrionmental Vars--------


class ROMSPainter(painter.Painter):
    """ROMS Painter"""
    def update(self):
        '''update ROMS specific files'''
        romscfg=self.cfg['ROMS']
        self.domid=romscfg['domid']
        self.varlist=cfgparser.cfg_get_varlist(self.cfg, 'ROMS', 'varlist')
        self.max_workers=int(romscfg['max_workers'])
        self.spatial_draw=romscfg.getboolean('spatial_draw')
        if self.spatial_draw:
            self.kwargs=romscfg['spatial_kwargs']
        self.form_anim=romscfg.getboolean('form_animation')
        self.run_sp=romscfg.getboolean('run_special')
        if self.run_sp:
            self.spjob=romscfg['special_job'].split(':')

    def load_data(self):
        '''load ROMS files according to domain ID'''
        utils.write_log('get roms hist file...')
        fn_stream=subprocess.check_output(
            'ls '+self.arch_root+'/roms_his_'+self.domid+'*', 
            shell=True).decode('utf-8')
           
        fn_list=fn_stream.split()
        self.file_list=fn_list
        
        if self.diff_flag:
            utils.write_log('get roms ref file...')
            fn_stream=subprocess.check_output(
                'ls '+self.arch_root2+'/roms_his_'+self.domid+'*', 
                shell=True).decode('utf-8')
            self.file_list2=fn_stream.split()
        self.roms_num_file=len(fn_list)
        utils.write_log('read %3d roms_his files' % len(fn_list))
        ds_roms=xr.load_dataset(fn_list[0])
        self.ds_roms=ds_roms
        # land sea mask
        self.lsmask =ds_roms['mask_rho']
        # vertical coordinate
        s_rho=ds_roms['s_rho'][:]
        self.isurf=s_rho.size-1
        # Get the latitude and longitude points
        self.lats, self.lons = ds_roms['lat_rho'].values, ds_roms['lon_rho'].values
        self.nrow, self.ncol = self.lats.shape
    def drive(self):
   
        import concurrent.futures
        self.update()
        #draw2d_map(fn_lst[0], self.fig_root, self.lats, self.lons, varname='Hwave')
        if self.spatial_draw:
            
            self.load_data()
            fn_lst=self.file_list
            
            for varname in self.varlist:
                vardic=const.ROMS_VAR_DIC[varname]
                coord=vardic['coord']
                para_dic={
                    'fig_root':self.fig_root, 
                    'lats':self.ds_roms[f'lat_{coord}'].values, 
                    'lons':self.ds_roms[f'lon_{coord}'].values, 
                    'add_shp':self.add_shp,'domid':self.domid, 
                    'varname':varname, 'watermark':self.watermark, 
                    'diff_flag':self.diff_flag,'kwargs':self.kwargs}
                if self.diff_flag:
                    para_dic['file_list2']=self.file_list2
                if varname in const.ROMS_STATIC_VAR:
                    fn=fn_lst[0]
                    draw2d_map(0, fn, para_dic)
                else:
                    with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                        futures = [executor.submit(
                            draw2d_map,  idx, file_path, para_dic) for idx, file_path in enumerate(fn_lst)] 
                        
                        concurrent.futures.wait(futures)
                        for idx, f in enumerate(futures):
                            if f.exception(): 
                                print(f.result(),idx)
                                exit()
        if self.form_anim:
            for varname in self.varlist:
                domid=self.domid
                if self.diff_flag:
                    fnpattern=f'roms.{domid}.{varname}-diff.'
                else:
                    fnpattern=f'roms.{domid}.{varname}.'
                fps=self.cfg['anim']['fps']
                utils.write_log(f'form animation@${fps}fps:{fnpattern}.mp4')
                '''
                This command takes a sequence of PNG image files with a frame rate of {fps} frames per second, 
                converts them to an MP4 video file with a frame rate of 8 frames per second, 
                and applies the libx264 codec and yuv420p pixel format to the output.
                '''
                spell=f'cd {self.fig_root};ffmpeg -y -r {fps} -pattern_type glob -i \'{fnpattern}*.png\''+\
                    f' -c:v libx264 -vf "fps=8,format=yuv420p" {fnpattern}mp4'
                utils.write_log(spell)
                subprocess.call(spell, shell=True)
                
        if self.run_sp:
            if self.spjob[0]=='extract_ts':
                self.load_data()
                site_path=os.path.join(
                    const.machine_dic['SRC_ROOT'],'db',self.spjob[1]+'.csv')
                stas=station.construct_stas(site_path)
                for sta in stas:
                    sta.ix,sta.iy=mathlib.locate_position(
                        sta.lat, sta.lon, self.lats, self.lons)
                fn_lst=self.file_list
                for varname in self.varlist:
                    para_dic={
                        'varname':varname,
                        'stations':stas}
                    utils.write_log(f'{print_prefix}extracting time series for {varname}')
                    stas_ts = pd.DataFrame()
                    with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                        futures = [executor.submit(
                            extract_ts, file_path, para_dic) for file_path in fn_lst] 
                            
                        concurrent.futures.wait(futures)
                        for idx, f in enumerate(futures):
                            if f.exception(): 
                                print(f.result(),idx)
                                exit()
                            if f.result() is not None:
                                stas_ts=pd.concat(
                                    [stas_ts, pd.DataFrame(f.result(), index=[0])], ignore_index=True)
                        utils.write_log(f'save to {self.arch_root}/stas_{varname}_ts.csv')
                        stas_ts.to_csv(f'{self.arch_root}/stas_{varname}_ts.csv')
            
            
            elif self.spjob[0]=='2d_max_draw':
                self.load_data()
                for var in self.varlist:
                    for fn in self.file_list:
                        utils.write_log(f'{print_prefix}aggregating 2d max {var} from {fn}')
                        ds = xr.open_dataset(fn)
                        try:
                            max_var.values = np.maximum(ds[var].values, max_var.values)
                        except UnboundLocalError:
                            max_var = ds[var]
                        ds.close()
                    max_var=max_var
                    att_dic=const.ROMS_VAR_DIC[var]
                    cmap=att_dic['cmap']
                    levels=att_dic['clvl']
                     
                    data_fn=self.arch_root+'/roms_max_'+var+'_'+self.domid+'.nc'
                    utils.write_log(f'max values save to {data_fn}')
                    max_var.to_netcdf(data_fn)
                    
                    if self.add_shp:
                        fig, ax = plt.subplots(figsize=[10, 8])
                        proj=ccrs.PlateCarree()
                        # Set projection and plot the main figure
                        ax = fig.add_axes([0.1, 0.05, 0.9, 0.85], projection=proj)
                        # add gridlines
                        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                                color='grey', linestyle='--', linewidth=0.5)
                        ax.coastlines()
                        ax.add_feature(cfeature.BORDERS, color='grey', linestyle='--')
                        try:
                            if att_dic['norm']==True:
                                norm=matplotlib.colors.BoundaryNorm(
                                    boundaries=levels,ncolors=cmap.N,extend='both')
                            plt.pcolormesh(
                                self.lons, self.lats, max_var.values[0,:,:],
                                norm=norm,
                                transform=ccrs.PlateCarree(), cmap=cmap)
                        except KeyError:
                            plt.pcolormesh(
                                self.lons, self.lats, max_var.values[0,:,:],
                                vmax=levels[-1], vmin=levels[0],
                                transform=ccrs.PlateCarree(), cmap=cmap)
                        plt.gca().set_facecolor("tan") 
                    else:
                        fig, ax = plt.subplots(figsize=[10, 8],frameon=True,constrained_layout=True)
                        plt.pcolormesh(
                            self.lons, self.lats, max_var.values[0,:,:],
                            vmax=levels[-1], vmin=levels[0],
                            cmap=cmap)
                        plt.gca().set_facecolor("black") 
                        
                    if self.cfg.has_option('special', '2d_max_draw_contour_lv'):
                        lv=self.cfg.getfloat('special', '2d_max_draw_contour_lv')
                        if self.add_shp:
                            contour = ax.contour(
                                self.lons, self.lats, max_var.values[0,:,:], levels=[lv], 
                                colors='white', linewidths=1.5, transform=ccrs.PlateCarree()) 
                        else:
                            contour = ax.contour(
                                self.lons, self.lats, max_var.values[0,:,:], levels=[lv], 
                                colors='white', linewidths=1.5)
                    plt.title(f'Maximum {var}', fontsize=MIDFONT)
                    plt.xticks( fontsize=SMFONT)
                    plt.yticks( fontsize=SMFONT)
                    cbar=plt.colorbar(ax=ax, shrink=0.8)
                    cbar.ax.tick_params(labelsize=SMFONT)
                    # Remove default matplotlib frame
                    plt.subplots_adjust(left=0, right=1, top=1, bottom=0) 
                    # Remove the frame around the plot
                    para_dic={
                    'fig_root':self.fig_root, 'domid':self.domid, 
                    'varname':var, 'watermark':self.watermark, 'diff_flag':self.diff_flag}
 
                    savefig(para_dic, var, 'dom_max')
                    plt.close()
                    del max_var
def extract_ts(fn, para_dic):
    varname,stas=para_dic['varname'],para_dic['stations']
    tfrm=io.get_file_tf(fn)
    tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values[0]
    rt_dic={
        'time':tfrm_str}
    
    for sta in stas:
        ix, iy=sta.ix, sta.iy
        var2d=io.get_single_var2d(fn, varname)
        rt_dic[sta.name]=var2d.values[ix, iy]
    return rt_dic 

def draw2d_map(idx, fn, para_dic):
    '''
    draw 2d spatial sigwave height map
    kwargs include:
        time:201809120000,201809130000 (yyyymmddHHMMSS)
        sub_dom:114.17,114.33,22.37,22.51 (lonmin,lonmax,latmin,latmax)
        cmaps:blre
        clvl:0,3,31
        rm_areamean
        rm_val:22.45,114.25 or 9.0 (lat,lon/val) 
    '''
    # read var 
    add_shp=para_dic['add_shp']
    varname, lats, lons=para_dic['varname'],para_dic['lats'], para_dic['lons']
    kwargs=para_dic['kwargs']
    if 'time' in kwargs and not(varname in const.ROMS_STATIC_VAR):
        [stime, etime]=utils.parse_kwargs(kwargs, 'time','strlist') 
        stime,etime=pd.to_datetime(stime),pd.to_datetime(etime)
        tfrm=io.get_file_tf(fn)
        if stime>tfrm or etime<tfrm:
            print(f'skip {tfrm.dt.strftime("%Y-%m-%d %H:%M:%S").values[0]}')
            return
    att_dic=const.ROMS_VAR_DIC[varname]
    
    utils.write_log(f'{print_prefix}ploting {varname} from {fn}')
    
    var2d=io.get_single_var2d(fn, varname)
    if 'sub_dom' in kwargs:
        [lonmin, lonmax, latmin, latmax]=utils.parse_kwargs(
            kwargs, 'sub_dom', 'float')
        var2d_subdom = xr.where(
            (lons>= lonmin) & (lons <= lonmax) &
            (lats>= latmin) & (lats <= latmax), var2d, np.nan)
    if 'rm_areamean' in kwargs:
        if 'sub_dom' in kwargs:
            var2d=var2d-var2d_subdom.mean()
        else:
            var2d=var2d-var2d.mean()
    if 'rm_val' in kwargs:
        rt_vals=utils.parse_kwargs(kwargs, 'rm_val', 'float')
        if len(rt_vals)==2:
            lat0, lon0=rt_vals[0], rt_vals[1]
            # Calculate the distance to each point
            distance = np.sqrt((lats - lat0)**2 + (lons - lon0)**2)
            # Find the indices of the minimum distance
            eta, xi = np.unravel_index(np.argmin(distance), distance.shape)
            var2d=var2d-var2d.sel(eta_rho=eta, xi_rho=xi)
        else:
            var2d=var2d-rt_vals[0]
    if varname in const.ROMS_STATIC_VAR:
        lmask=io.get_single_var2d(fn, 'mask_rho')
        var2d=xr.where(lmask>0, var2d, np.nan) 
    
    cmap=att_dic['cmap']
    levels=att_dic['clvl']
    if 'cmap' in kwargs:
        cmapname=utils.parse_kwargs(kwargs, 'cmap')
        cmap=const.CMAPS[cmapname]
    if 'clvl' in kwargs:
        largs=utils.parse_kwargs(kwargs, 'clvl','float')
        levels=np.linspace(largs[0], largs[1], int(largs[2]))
   
    title_txt=f'{varname} ({att_dic["units"]})'
 
    if add_shp:
        fig= plt.figure(
            figsize=[10, 8])
        proj=ccrs.PlateCarree()
        # Set projection and plot the main figure
        ax = fig.add_axes([0.1, 0.05, 0.9, 0.85], projection=proj)
        # add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                color='grey', linestyle='--', linewidth=0.5)
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=':')
    else:
        fig, ax = plt.subplots(figsize=[10, 8])
   
    if para_dic['diff_flag']:
        var2d_ref=io.get_single_var2d(para_dic['file_list2'][idx], varname)
        var2d.values=var2d.values-var2d_ref.values
        cmap=att_dic['diff_cmap']
        levels=att_dic['diff_clvl']
        title_txt=varname+'_Diff'
        #arr=var2d.values[0,:,:]
        arr=var2d.values
        # Mask the nan values in the array
        plt.contourf(
            lons, lats, arr,
            levels=levels, extend='both', 
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
        
            plt.gca().set_facecolor("tan") 
        else:
        #ax=set_canvas_common(lats, lons, watermark=para_dic['watermark'])
            plt.pcolormesh(
                lons, lats, var2d.values,
                vmax=levels[-1], vmin=levels[0],
                cmap=cmap)
            plt.xticks( fontsize=SMFONT)
            plt.yticks( fontsize=SMFONT)

            plt.gca().set_facecolor("black")  

    if 'sub_dom' in kwargs:
        plt.xlim(lonmin, lonmax)
        plt.ylim(latmin, latmax) 
    # Add a color bar
    cbar=plt.colorbar(ax=ax, shrink=0.9, extendfrac='auto')
    cbar.ax.tick_params(labelsize=SMFONT) 
    
    if varname in const.ROMS_STATIC_VAR: 
        plt.title(title_txt, fontsize=MIDFONT)
        savefig(
            para_dic, varname, 'static')
    else:
        tfrm=var2d['ocean_time']
        #tfrm=var2d['ocean_time'][0]
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=title_txt+'@'+tfrm_str
        plt.title(title_txt, fontsize=MIDFONT)
        savefig(
            para_dic, varname, tfrm.dt.strftime('%Y%m%d%H%M').values)
    
    plt.close()



def savefig(para_dic, varname, time_str,bbox_inches=''):
    '''save fig portal'''
    save_dir=para_dic['fig_root']
    diff_flag=para_dic['diff_flag']
    domid=para_dic['domid']
    if diff_flag:
        fn=f'roms.{domid}.{varname}-diff.{time_str}.png'
    else:
        fn=f'roms.{domid}.{varname}.{time_str}.png'
    utils.write_log(f'saving {fn} to {save_dir}') 
    if not(os.path.isdir(save_dir)):
            os.mkdir(save_dir)        
    if bbox_inches=='':
        plt.savefig(os.path.join(save_dir,fn), dpi=100)
    else:
        plt.savefig(os.path.join(save_dir,fn), dpi=100, bbox_inches=bbox_inches)


def set_canvas_common(lats, lons, watermark=False):
    '''set common properties of the canvas'''
    lonmin, lonmax=lons.min(), lons.max()
    latmin, latmax=lats.min(), lats.max()

    proj = ccrs.Mercator(
        central_longitude=115., min_latitude=-80.0, max_latitude=84.0, 
        globe=None, latitude_true_scale=23.0, false_easting=0.0, 
        false_northing=0.0, scale_factor=None) 
    
    fig = plt.figure(figsize=[10.24, 7.68])
    
    # Set projection and plot the main figure
    ax = fig.add_axes([0.08, 0.01, 0.8, 0.94], projection=proj)

    # Set figure extent
    ax.set_extent(
        [lonmin, lonmax, latmin, latmax],
        crs=ccrs.PlateCarree())
    if watermark:
        ax.text(0.5, 0.5, const.WATERMARK_TXT, transform=ax.transAxes,
        fontsize=12, color='gray', alpha=0.35,
        ha='center', va='center', rotation=30) 
    # add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
        color='grey', linestyle='--', linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False
    

    #amdn_shp=shpreader.Reader(const.SHP_ROOT+SHP_SELECT).geometries()
    #amdn_shp_outer=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
    #ax.add_geometries(
    #    amdn_shp_outer, ccrs.PlateCarree(),
    #    facecolor='none', edgecolor='black',linewidth=1, zorder = 1)
    #ax.coastlines()
    
    # plot shp boundaries
    #ax.add_geometries(
    #    amdn_shp, ccrs.PlateCarree(),
    #    facecolor='none', edgecolor='black',linewidth=.5, zorder = 1)

    return ax    
