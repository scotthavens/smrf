"""
Functions to output as a netCDF
"""

# from scipy import stats
import logging
import os
from datetime import datetime

import netCDF4 as nc
import numpy as np
from spatialnc.proj import add_proj, add_proj_from_file

from smrf import __version__
from smrf.utils import utils

# import pandas as pd


class output_netcdf():
    """
    Class output_netcdf() to output values to a netCDF file
    """

    type = 'netcdf'
    fmt = '%Y-%m-%d %H:%M:%S'
    cs = (6, 10, 10)

    def __init__(self, variable_list, topo, time, outConfig):
        """
        Initialize the output_netcdf() class

        Args:
            variable_list: list of dicts, one for each variable
            topo: loadTopo instance
        """

        self._logger = logging.getLogger(__name__)

        # go through the variable list and make full file names
        for v in variable_list:
            variable_list[v]['file_name'] = \
                str(variable_list[v]['out_location'] + '.nc')

        self.variable_list = variable_list

        # process the time section
        self.run_time_step = int(time['time_step'])
        self.out_frequency = int(outConfig['frequency'])
        self.outConfig = outConfig

        # determine the x,y vectors for the netCDF file
        x = topo.x
        y = topo.y
        self.mask = topo.mask

        dimensions = ('time', 'y', 'x')
        self.date_time = {}
        now_str = datetime.now().strftime(self.fmt)

        # Retrieve projection information from topo
        map_meta = add_proj_from_file(topo.topoConfig['filename'])

        for v in self.variable_list:

            f = self.variable_list[v]

            if os.path.isfile(f['file_name']):
                self._logger.warning('Opening {}, data may be overwritten!'
                                  .format(f['file_name']))

                # open in append mode
                s = nc.Dataset(f['file_name'], 'a')
                h = '[{}] Data added or updated'.format(now_str)
                setattr(s, 'last_modified', h)

                if 'projection' not in s.variables.keys():
                    s = add_proj(s, map_meta=map_meta)

            else:
                self._logger.debug('Creating %s' % f['file_name'])
                s = nc.Dataset(f['file_name'], 'w', format='NETCDF4',
                                                    clobber=True)

                # add dimensions
                s.createDimension(dimensions[0], None)
                s.createDimension(dimensions[1], y.shape[0])
                s.createDimension(dimensions[2], x.shape[0])

                # create the variables
                s.createVariable('time', 'f', (dimensions[0]))
                s.createVariable('y', 'f', dimensions[1])
                s.createVariable('x', 'f', dimensions[2])
                s.createVariable(f['variable'], 'f',
                                 (dimensions[0], dimensions[1], dimensions[2]),
                                 chunksizes=self.cs)

                # # define some attributes
                s.variables['time'].setncattr(
                        'units',
                        'hours since {}'.format(time['start_date']))
                s.variables['time'].setncattr('calendar', 'standard')
                s.variables['time'].setncattr('time_zone', time['time_zone'])
                s.variables['time'].setncattr('long_name', 'time')

                # the y variable attributes
                s.variables['y'].setncattr('units','meters')
                s.variables['y'].setncattr('description','UTM, north south')
                s.variables['y'].setncattr('long_name','y coordinate')

                # the x variable attributes
                s.variables['x'].setncattr('units','meters')
                s.variables['x'].setncattr('description','UTM, east west')
                s.variables['x'].setncattr('long_name','x coordinate')

                # the variable attributes
                s.variables[f['variable']].setncattr(
                        'module',
                        f['module'])
                s.variables[f['variable']].setncattr(
                        'units',
                        f['info']['units'])
                s.variables[f['variable']].setncattr(
                        'long_name',
                        f['info']['long_name'])


                # define some global attribute
                s.setncattr_string('Conventions', 'CF-1.6')
                s.setncattr_string('dateCreated', now_str)
                s.setncattr_string('title', 'Distributed {0} data from SMRF'
                                            ''.format(f['info']['long_name']))
                s.setncattr_string('history', ('[{}] Create netCDF4 file'
                                              '').format(now_str))
                s.setncattr_string('institution', ('USDA Agricultural Research'
                                                   ' Service, Northwest'
                                                   ' Watershed Research'
                                                   ' Center'))
                s = add_proj(s, map_meta=map_meta)

                s.setncattr_string('references',
                        'Online documentation smrf.readthedocs.io;'
                        ' https://doi.org/10.1016/j.cageo.2017.08.016')

                s.variables['y'][:] = y
                s.variables['x'][:] = x

            s.setncattr_string('source',
                    'SMRF {}'.format(utils.getgitinfo()))
            s.close()

    def output(self, variable, data, date_time):
        """
        Output a time step

        Args:
            variable: variable name that will index into variable list
            data: the variable data
            date_time: the date time object for the time step
        """

        self._logger.debug('{0} Writing variable {1} to netCDF'
                           .format(date_time, variable))

#         f = self.variable_list[variable]['nc_file']
        f = nc.Dataset(self.variable_list[variable]['file_name'],
                       'a', 'NETCDF4')

        # the current time integer
        times = f.variables['time']
        t = nc.date2num(date_time.replace(tzinfo=None),
                        times.units,
                        times.calendar)

        if len(times) != 0:
            index = np.where(times[:] == t)[0]
            if index.size == 0:
                index = len(times)
            else:
                index = index[0]
        else:
            index = len(times)

#         self.date_time[variable] = np.append(self.date_time[variable], pd.to_datetime(date_time))

            # insert the time
        f.variables['time'][index] = t

        # insert the data
        if self.outConfig['mask_output']:
            f.variables[variable][index, :] = data*self.mask
        else:
            f.variables[variable][index, :] = data

        # synce the data to disk to that it can be viewed immediately
#         f.sync()
        f.close()
