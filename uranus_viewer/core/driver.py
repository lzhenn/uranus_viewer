#/usr/bin/env python3
''' KERNEL MODULE '''

# ---imports---
from ..lib import utils
from . import roms_painter, wrf_painter, swan_painter
# ---Module regime consts and variables---


# ---Classes and Functions---
class Driver:
    def __init__(self, cfg):
        self.cfg=cfg
        self.run_wrf=cfg['WRF'].getboolean('wrf_flag')
        self.run_roms=cfg['ROMS'].getboolean('roms_flag')
        self.run_swan=cfg['SWAN'].getboolean('swan_flag')
    def drive(self):
        if self.run_wrf:
            utils.write_log('WRF Postprocessing...')
            painter=wrf_painter.WRFPainter(self.cfg)
            painter.drive()
        if self.run_roms:
            utils.write_log('ROMS Postprocessing...')
            painter=roms_painter.ROMSPainter(self.cfg)
            painter.drive() 
        if self.run_swan:
            utils.write_log('SWAN Postprocessing...')
            painter=swan_painter.SWANPainter(self.cfg)
            painter.drive()

# ---Unit test---
if __name__ == '__main__':
    pass
