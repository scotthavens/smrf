
import os
import unittest
from copy import deepcopy
from glob import glob

import requests
from inicheck.tools import cast_all_variables
from inicheck.utilities import pcfg

from smrf.framework.model_framework import can_i_run_smrf, run_smrf
from tests.test_configurations import SMRFTestCase


class TestLoadGrid(SMRFTestCase):

    def compare_hrrr_gold(self, out_dir):
        """
        Compare the model results with the gold standard

        Args:
            out_dir: the output directory for the model run
        """

        output_dir = os.path.join(self.test_dir, out_dir)
        s = os.path.join(output_dir, '*.nc')
        file_names = glob(os.path.realpath(s))

        # path to the gold standard
        gold_path = os.path.realpath(
            os.path.join(self.test_dir, 'RME', 'gold_hrrr'))

        for file_name in file_names:
            nc_name = file_name.split('/')[-1]
            print('  Comparing {}'.format(nc_name))

            v = nc_name.split('.')[0]
            gold_file = os.path.join(gold_path, nc_name)

            self.compare_netcdf_files(gold_file, file_name)

    def test_grid_wrf(self):
        """ WRF NetCDF loading """

        config = deepcopy(self.base_config)
        del config.raw_cfg['csv']

        adj_config = {
            'gridded': {
                'data_type': 'wrf',
                'wrf_file': './RME/gridded/WRF_test.nc',
                'zone_number': '11',
                'zone_letter': 'N'
            },
            'system': {
                'threading': 'False',
                'log_file': './RME/output/log.txt'
            },
            'precip': {
                'station_adjust_for_undercatch': 'False'
            },
            'thermal': {
                'correct_cloud': 'False',
                'correct_veg': 'True'
            },
            'time': {
                'start_date': '2015-03-03 00:00',
                'end_date': '2015-03-03 04:00'
            }
        }

        config.raw_cfg.update(adj_config)

        # set the distribution to grid, thermal defaults will be fine
        for v in self.dist_variables:
            config.raw_cfg[v]['grid_mask'] = 'False'

        # fix the time to that of the WRF_test.nc
        config.apply_recipes()
        config = cast_all_variables(config, config.mcfg)

        # ensure that the recipes are used
        assert 'station_adjust_for_undercatch' not in config.cfg['precip'].keys(
        )
        self.assertTrue(config.cfg['thermal']['correct_cloud'] == False)
        self.assertTrue(config.cfg['thermal']['correct_veg'] == True)

        result = can_i_run_smrf(config)
        self.assertTrue(result)

    def test_grid_hrrr_local(self):
        """ HRRR grib2 loading with local elevation gradient """

        config = deepcopy(self.base_config)
        del config.raw_cfg['csv']

        adj_config = {
            'gridded': {
                'data_type': 'hrrr_grib',
                'hrrr_directory': './RME/gridded/hrrr_test/',
                'zone_number': 11,
                'zone_letter': 'N'
            },
            'time': {
                'start_date': '2018-07-22 16:00',
                'end_date': '2018-07-22 20:00'
            },
            'system': {
                'threading': False,
                'log_file': './RME/output/log.txt'
            },
            'air_temp': {
                'grid_local': True,
                'grid_local_n': 25
            },
            'vapor_pressure': {
                'grid_local': True,
                'grid_local_n': 25
            },
            'precip': {
                'grid_local': True,
                'grid_local_n': 25,
                'precip_temp_method': 'dew_point'
            },
            'thermal': {
                'correct_cloud': True,
                'correct_veg': True
            },
            'albedo': {
                'grain_size': 300.0,
                'max_grain': 2000.0
            }
        }
        config.raw_cfg.update(adj_config)

        # set the distribution to grid, thermal defaults will be fine
        for v in self.dist_variables:
            config.raw_cfg[v]['distribution'] = 'grid'
            config.raw_cfg[v]['grid_mask'] = 'False'

        config.apply_recipes()
        config = cast_all_variables(config, config.mcfg)

        # ensure that the recipes are used
        self.assertTrue(
            'station_adjust_for_undercatch' not in config.cfg['precip'].keys())
        self.assertTrue(config.cfg['thermal']['correct_cloud'] == True)
        self.assertTrue(config.cfg['thermal']['correct_veg'] == True)

        result = can_i_run_smrf(config)
        self.assertTrue(result)

        self.compare_hrrr_gold(config.raw_cfg['output']['out_location'][0])

    def test_grid_netcdf(self):
        """ Generic NetCDF loading """

        config = deepcopy(self.base_config)
        del config.raw_cfg['csv']

        generic_grid = {'data_type': 'netcdf',
                        'netcdf_file': './RME/gridded/netcdf_test.nc',
                        'zone_number': '11',
                        'zone_letter': 'N',
                        'air_temp': 'air_temp',
                        'vapor_pressure': 'vapor_pressure',
                        'precip': 'precip',
                        'wind_speed': 'wind_speed',
                        'wind_direction': 'wind_direction',
                        'thermal': 'thermal',
                        'cloud_factor': 'cloud_factor'}
        config.raw_cfg['gridded'] = generic_grid
        config.raw_cfg['system']['time_out'] = '25'
        config.raw_cfg['system']['queue_max_values'] = '2'
        # Doesn't work with true
        config.raw_cfg['system']['threading'] = 'False'

        # set the distribution to grid, thermal defaults will be fine
        for v in self.dist_variables:
            config.raw_cfg[v]['distribution'] = 'grid'
            config.raw_cfg[v]['grid_mask'] = 'False'

        config.raw_cfg['thermal']['correct_cloud'] = 'False'
        config.raw_cfg['thermal']['correct_veg'] = 'True'

        # fix the time to that of the WRF_test.nc
        config.raw_cfg['time']['start_date'] = '2015-03-03 00:00'
        config.raw_cfg['time']['end_date'] = '2015-03-03 04:00'

        config.apply_recipes()
        config = cast_all_variables(config, config.mcfg)

        # ensure that the recipes are used
        self.assertTrue(
            'station_adjust_for_undercatch' not in config.cfg['precip'].keys())
        self.assertTrue(config.cfg['thermal']['correct_cloud'] == False)
        self.assertTrue(config.cfg['thermal']['correct_veg'] == True)

        result = can_i_run_smrf(config)
        self.assertTrue(result)
