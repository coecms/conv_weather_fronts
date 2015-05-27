#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

lats = np.arange(-90.0, 90.0, 0.75)
longs = np.arange(0.0, 359.25, 0.75)

infile_name = 'rec_front_1979_01.v26.nc'
outfile_name = 'out.nc'

def create_maps(nc_var, nnumber, npoint, lats, longs):
    """
    Input: netCDF variable, timestep, number of fronts, number of points
        numpy array with latitudes, numpy array with longitudes
    Output: 2D numpy array of ints., 3x2D arrays of floats
    """
    front_map = np.ones((len(lats), len(longs)), dtype=np.int) * -9999
    t_map = np.ones((len(lats), len(longs)), dtype=np.float) * -9999.0
    u_map = t_map[:, :]
    v_map = t_map[:, :]
    for n in xrange(nnumber):
        for p in xrange(npoint):
            lat, lon, speed_t, speed_u, speed_v = nc_var[n, p, :]
            if lat < -1000 :
                break
            latidx = np.argmin(abs(lats - lat))
            lonidx = np.argmin(abs(longs - lon))
            if ((lon > longs[-1]) and
                (360.0+longs[0] - lon < lon - longs[-1])):
                lonidx = 0
            front_map[latidx, lonidx] = n
            t_map[latidx, lonidx] = speed_t
            u_map[latidx, lonidx] = speed_u
            v_map[latidx, lonidx] = speed_v
    return (front_map, t_map, u_map, v_map)


# input data file
infile = nc.Dataset(infile_name, 'r')
ntime = len(infile.dimensions['time'])
nnumber = len(infile.dimensions['number'])
npoint = len(infile.dimensions['point'])
ndata = len(infile.dimensions['data'])

# output data file
outfile = nc.Dataset(outfile_name, 'w')

# Create dimensions:
# Time, unlimited
outfile.createDimension('time', size=None)
timevar = outfile.createVariable(
    'time', infile.variables['time'].datatype, ('time',))
timevar.units = infile.variables['time'].units

# Latitude
outfile.createDimension('latitude', size=len(lats))
latsvar = outfile.createVariable(
    'latitude', np.float32, ('latitude',))
latsvar.long_name = 'Latitude'
latsvar.units = 'degree_north'

# Longitude
outfile.createDimension('longitude', size=len(longs))
longsvar = outfile.createVariable(
    'longitude', np.float32, ('longitude',))
longsvar.long_name = 'Longitude'
longsvar.units = 'degree_east'

# Create Variables
# Cold / Warm / Stationary Fronts
cf_var = outfile.createVariable(
    'cold_fronts', np.int16, ('time', 'latitude', 'longitude'),
    fill_value = -9999)
wf_var = outfile.createVariable(
    'warm_fronts', np.int16, ('time', 'latitude', 'longitude'),
    fill_value = -9999)
sf_var = outfile.createVariable(
    'stat_fronts', np.int16, ('time', 'latitude', 'longitude'),
    fill_value = -9999)

# Theta Gradient
t_var = outfile.createVariable(
    'thetaw_gradient', np.float32, ('time', 'latitude', 'longitude'),
    fill_value = -9999.0)

# Speed U and V
u_var = outfile.createVariable(
    'u_speed', np.float32, ('time', 'latitude', 'longitude'),
    fill_value = -9999.0)
v_var = outfile.createVariable(
    'v_speed', np.float32, ('time', 'latitude', 'longitude'),
    fill_value = -9999.0)

# Initialise the dimensions
timevar[:] = infile.variables['time'][:]
latsvar[:] = lats
longsvar[:] = longs

cf_data = np.zeros((nnumber, npoint, 5), dtype=np.float)
wf_data = np.zeros_like(cf_data)
sf_data = np.zeros_like(cf_data)

for t in xrange(ntime):
    cf_data[:, :, :] = infile.variables['cold_fronts'][t, :, :, :]
    wf_data[:, :, :] = infile.variables['warm_fronts'][t, :, :, :]
    sf_data[:, :, :] = infile.variables['stat_fronts'][t, :, :, :]
    cf_map, t_map_cf, u_map_cf, v_map_cf = \
              create_maps(cf_data, nnumber, npoint, lats, longs)
    wf_map, t_map_wf, u_map_wf, v_map_wf = \
              create_maps(wf_data, nnumber, npoint, lats, longs)
    sf_map, t_map_sf, u_map_sf, v_map_sf = \
              create_maps(sf_data, nnumber, npoint, lats, longs)

    cf_var[t, :, :] = cf_map[:, :]
    wf_var[t, :, :] = wf_map[:, :]
    sf_var[t, :, :] = sf_map[:, :]

    t_map = np.where(t_map_wf > -9998.0, t_map_wf,
                              np.where(t_map_cf != 0.0, t_map_cf,
                                       t_map_sf))
    u_map = np.where(u_map_wf > -9998.0, u_map_wf,
                              np.where(u_map_cf != 0.0, u_map_cf,
                                       u_map_sf))
    v_map = np.where(v_map_wf > -9998.0, v_map_wf,
                              np.where(v_map_cf != 0.0, v_map_cf,
                                       v_map_sf))
    t_var[t, :, :] = t_map[:, :]
    u_var[t, :, :] = u_map[:, :]
    v_var[t, :, :] = v_map[:, :]
    
    print "{}/{}".format(t, ntime)
    

outfile.close()
infile.close()
