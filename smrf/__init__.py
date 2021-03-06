# -*- coding: utf-8 -*-
'''
SMRF - Snow Modeling Resource Framework
'''
__version__ = '0.9.1'

import os

__core_config__ = os.path.abspath(os.path.dirname(__file__) + '/framework/CoreConfig.ini')
__recipes__ = os.path.abspath(os.path.dirname(__file__) + '/framework/recipes.ini')
__config_changelog__ = os.path.abspath(os.path.dirname(__file__) + '/framework/changelog.ini')


__config_titles__ = {
                "topo":"Files for DEM and vegetation",
                "time": "Dates to run model",
                "stations":"Stations to use",
                "csv":"CSV section configurations",
                "mysql":"MySQL database",
                "gridded":"Gridded datasets configurations",
                "air_temp":"Air temperature distribution",
                "vapor_pressure":"Vapor pressure distribution",
                "wind": "Wind speed and wind direction distribution",
                "precip": "Precipitation distribution",
                "albedo":"Albedo distribution",
                "cloud_factor": "Cloud Factor - Fraction used to limit solar radiation Cloudy (0) - Sunny (1)",
                "solar":"Solar radiation",
                "thermal":"Incoming thermal radiation",
                "soil_temp":"Soil temperature",
                "output":"Output variables",
                "system":"System variables and Logging"
}

from . import data, distribute, envphys, framework, output, spatial, utils # isort:skip

__config_header__ = utils.utils.getConfigHeader()
__config_checkers__ = 'utils.utils'
# -*- coding: utf-8 -*-
