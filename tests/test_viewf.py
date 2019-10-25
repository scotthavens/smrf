from spatialnc import ipw
import matplotlib.pyplot as plt
import numpy as np

from smrf.data import loadTopo
from smrf.utils.topo import hor1d

from tests.test_configurations import SMRFTestCase


class TestViewf(SMRFTestCase):

    def test_hor1d(self):

        topo_config = {
            'basin_lon': -116.7547,
            'basin_lat': 43.067,
            'filename': 'tests/RME/topo/topo.nc',
            'type': 'netcdf',
            'threading': False
        }

        # ~/code/ipw/src/bin/topocalc/horizon/hor1d/hor1d -a 90 dem.ipw > hor1d.ipw
        # hor1d in the East direction
        hipw = ipw.IPW('tests/RME/hor1d.ipw')
        hipw = hipw.bands[0].data[0, :]

        # IPW topo calc
        topo = loadTopo.topo(topo_config, calcInput=True,
                             tempDir='tests/RME/output')

        dx = np.mean(np.diff(topo.x))
        h = hor1d.hor1f(topo.dem[0, :])
        hcos = hor1d.horval(topo.dem[0, :], dx, h)
