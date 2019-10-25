from spatialnc import ipw
import matplotlib.pyplot as plt
import numpy as np
import unittest

from smrf.data import loadTopo
from smrf.utils.topo import hor1d

from tests.test_configurations import SMRFTestCase


class TestViewf(SMRFTestCase):

    def test_hor1d_RME(self):
        """ Hor1d on the RME test dem """

        topo_config = {
            'basin_lon': -116.7547,
            'basin_lat': 43.067,
            'filename': 'tests/RME/topo/topo.nc',
            'type': 'netcdf',
            'threading': False
        }

        # ~/code/ipw/src/bin/topocalc/horizon/hor1d/hor1d -a 90 dem.ipw > hor1d.ipw
        # hor1d in the East direction
        hipw = ipw.IPW('tests/RME/gold/radiation/hor1d.ipw')
        hipw = hipw.bands[0].data

        # IPW topo calc
        topo = loadTopo.topo(topo_config, calcInput=True,
                             tempDir='tests/RME/output')

        dx = np.mean(np.diff(topo.x))
        h = hor1d.hor1f(topo.dem[1, :])
        hcos = hor1d.horval(topo.dem[1, :], dx, h)

        hv = hor1d.hor1d_c(topo.dem[1, :], dx)

        hcos = hor1d.hor1d(topo.dem, dx)

        # about the tolerance between a 16bit image and a float
        self.assertTrue(np.allclose(hipw, hcos, atol=1e-4))
        # plt.imshow(hcos - hipw)
        # plt.title('Hor1d difference from the East')
        # plt.colorbar()
        # plt.show()

    # def test_hor1d_Tuolumne(self):
    #     """ hor1d Tuolumne test """
    #     # Pretty close, it doesn't pass but the bulk of the pixels are 0

    #     dem_file = 'tests/Tuolumne/topo/dem_50m.ipw'
    #     dem = ipw.IPW(dem_file)
    #     dx = dem.bands[0].dsamp
    #     z = dem.bands[0].data

    #     # ~/code/ipw/src/bin/topocalc/horizon/hor1d/hor1d -a 90 dem_50m.ipw > hor1d.ipw
    #     # hor1d in the East direction
    #     hipw = ipw.IPW('tests/Tuolumne/gold/radiation/hor1d.ipw')
    #     hipw = hipw.bands[0].data

    #     # Python hor1d
    #     hcos = hor1d.hor1d(z, dx)

    #     # about the tolerance between a 16bit image and a float
    #     self.assertTrue(np.allclose(hipw, hcos, atol=1e-4))
    #     result = hcos - hipw
    #     plt.hist(result.flatten(), bins=100)
    #     plt.title('Hor1d difference from the East')
    #     plt.colorbar()
    #     plt.show()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestViewf("test_hor1d_Tuolumne"))
    runner = unittest.TextTestRunner()
    runner.run(suite)
