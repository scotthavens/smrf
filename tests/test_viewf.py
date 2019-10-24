from smrf.data import loadTopo
from smrf.utils.topo import hor1d

from tests.test_configurations import SMRFTestCase

import matplotlib.pyplot as plt


class TestViewf(SMRFTestCase):

    def test_hor1d(self):

        topo_config = {
            'basin_lon': -116.7547,
            'basin_lat': 43.067,
            'filename': 'tests/RME/topo/topo.nc',
            'type': 'netcdf',
            'threading': False
        }

        # IPW topo calc
        topo = loadTopo.topo(topo_config, calcInput=True,
                             tempDir='tests/RME/output')

        h = hor1d.hor1f_simple(topo.dem[0, :])

        h
