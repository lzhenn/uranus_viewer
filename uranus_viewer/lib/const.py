import cmaps 
import numpy as np
PKG_NAME = 'uranus'
machine_dic={
    'SRC_ROOT':'/home/lzhenn/array74/workspace/uranus_viewer/uranus_viewer/',
    'SHP_ROOT':'/home/lzhenn/array74/workspace/njord_pipeline/lib/shp/',
    'IBTRCK_FN':'/home/lzhenn/array74/data/ecf/ibtracs.WP.list.v04r00.csv',
}
WATERMARK_TXT='Created by Zhenning LI (zhenningli91@gmail.com)'

CMAPS={
    'blre':cmaps.BlWhRe,'blvi200':cmaps.BlGrYeOrReVi200,
    'drywet':cmaps.GMT_drywet,'wbgyr':cmaps.WhiteBlueGreenYellowRed,
    'precip':cmaps.precip2_17lev,'cmocean':cmaps.cmocean_deep,
    'nice':cmaps.nice_gfdl}
WRF_VAR_DIC={
    'rh2':{
        'units':'%',
        'clvl':np.linspace(40,100,121),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['drywet'],
        'diff_cmap':CMAPS['drywet']
    },
    'SST':{
        'units':'K',
        'clvl':np.linspace(295,305,51),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },
    'TSK':{
        'units':'K',
        'clvl':np.linspace(285,315,121),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },
    'U10':{
        'units':'m/s',
        'clvl':np.linspace(-20,20,81),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blre'],
        'diff_cmap':CMAPS['blvi200']
    },
    'uv10':{
        'units':'m/s',
        'clvl':np.linspace(0,40,81),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['wbgyr'],
        'diff_cmap':CMAPS['blvi200']
    },
    'ua':{
        'units':'m/s',
        'clvl':np.linspace(-20,20,81),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':'coolwarm',
        'diff_cmap':CMAPS['blvi200']
    },
    'Q2':{
        'units':'kg kg-1',
        'clvl':np.linspace(0.005,0.025,101),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['drywet'],
        'diff_cmap':CMAPS['blvi200']
    },  
    'slp':{
        'units':'hPa',
        'clvl':np.linspace(970,1020,101),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },
    'PMSL':{
        'units':'Pa',
        'clvl':np.linspace(97000,102000,101),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },  
    'LANDMASK':{
        'units':'',
        'clvl':np.linspace(0,1,1),
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },
    'precip_1':{
        'units':'mm/hr',
        'clvl': [0.1,0.5,1,2,3,5,7,10,15,20,30,40,50,60,70,80,100] ,
        'norm':True,
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['precip'],
        'diff_cmap':CMAPS['blvi200']
    },
    'precip':{
        'units':'mm',
        'clvl': [0.1,0.5,1,2,3,5,7,10,15,20,30,40,50,60,70,80,100] ,
        'norm':True,
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['precip'],
        'diff_cmap':CMAPS['blvi200']
    },
    'precip_24':{
        'units':'mm',
        'clvl': [0.1,0.5,1,3,5,10,15,20,30,40,60,80,100,120,150,200,250],
        'diff_clvl':np.linspace(-50,50,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },
}

WRF_VAR_DIC['va']=WRF_VAR_DIC['ua']

ROMS_VAR_DIC={
    'Hwave':{
        'coord':'rho',
        'units':'m',
        'clvl':np.linspace(0,4.0,81),
        'diff_clvl':np.linspace(-1.0,1.0,51),
        'cmap':CMAPS['wbgyr'],
        'diff_cmap':CMAPS['blre']
    },
    'zeta':{
        'coord':'rho','units':'m',
        'clvl':np.linspace(0,4.0,81),
        'diff_clvl':np.linspace(-0.5,0.5,51),
        'cmap':CMAPS['wbgyr'],
        'diff_cmap':'seismic'
    },
    'ubar':{
        'coord':'u','units':'m/s',
        'clvl':np.linspace(-1.0,1.0,101),
        'diff_clvl':np.linspace(-0.05,0.05,51),
        'cmap':CMAPS['blvi200'],
        'diff_cmap':CMAPS['blvi200']
    },
    'h':{
        'coord':'rho','units':'m',
        'clvl':np.linspace(0,40,81),
        'diff_clvl':np.linspace(-0.05,0.05,51),
        'cmap':CMAPS['wbgyr'],
        'diff_cmap':CMAPS['cmocean']
    },
}   

ROMS_STATIC_VAR=['h']