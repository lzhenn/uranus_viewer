# conda create -n uranus_viewer python=3.10
# conda activate uranus_viewer
# conda install numpy pandas scipy netCDF4 xarray wrf-python cftime 
python -m uranus_viewer.run_viewer $1
