from spatialnc import ipw
import matplotlib.pyplot as plt
import numpy as np
import unittest

from smrf.data import loadTopo
from smrf.utils.topo import hor1d, viewf

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

        # compare the 1D horizon functions
        h = hor1d.hor1f(topo.dem[1, :])
        hcos = hor1d.horval(topo.dem[1, :], dx, h)

        # C version
        hv = hor1d.hor1d_c(topo.dem[1, :], dx)
        self.assertTrue(np.all(hv == hcos))

        # 2D horizon functions
        hcos = hor1d.hor1d(topo.dem, dx)
        hcos2 = hor1d.hor2d_c(topo.dem, dx)
        self.assertTrue(np.all(hcos == hcos2))

        # about the tolerance between a 16bit image and a float
        self.assertTrue(np.allclose(hipw, hcos, atol=1e-4))
        # plt.imshow(hcos - hipw)
        # plt.title('Hor1d difference from the East')
        # plt.colorbar()
        # plt.show()

    def test_hor1d_Tuolumne(self):
        """ hor1d Tuolumne test """
        # Pretty close, it doesn't pass but the bulk of the pixels are 0

        dem_file = 'tests/Tuolumne/topo/dem_50m.ipw'
        dem = ipw.IPW(dem_file)
        dx = dem.bands[0].dsamp
        z = dem.bands[0].data

        # convert to double
        z = z.astype(np.double)

        # ~/code/ipw/src/bin/topocalc/horizon/hor1d/hor1d -a 90 dem_50m.ipw > hor1d.ipw
        # hor1d in the East direction
        hipw = ipw.IPW('tests/Tuolumne/gold/radiation/hor1d.ipw')
        hipw = hipw.bands[0].data

        # Python hor1d
        hcos = hor1d.hor2d_c(z, dx)

        # result = hcos - hipw
        # plt.hist(result.flatten(), bins=100)
        # plt.title('Hor1d difference from the East')
        # # plt.colorbar()
        # plt.show()

    def test_viewf_RME(self):
        """ Test the view factor for RME """

        topo_config = {
            'basin_lon': -116.7547,
            'basin_lat': 43.067,
            'filename': 'tests/RME/topo/topo.nc',
            'type': 'netcdf',
            'threading': False
        }

        # ~/code/ipw/src/bin/topocalc/horizon/hor1d/hor1d -a 90 dem.ipw > hor1d.ipw
        # hor1d in the East direction
        # hipw = ipw.IPW('tests/RME/gold/radiation/hor1d.ipw')
        # hipw = hipw.bands[0].data

        # IPW topo calc for sky view factor
        topo = loadTopo.topo(topo_config, calcInput=True,
                             tempDir='tests/RME/output')
        dx = np.mean(np.diff(topo.x))
        ipw_svf = topo.stoporad_in.bands[3].data
        ipw_tcf = topo.stoporad_in.bands[3].data

        # Calcualte the sky view factor
        svf, tcf = viewf.viewf(
            topo.dem, dx, slope=topo.slope, aspect=topo.aspect)

        plt.imshow(ipw_svf - svf)
        plt.title('Sky view factor')
        plt.colorbar()
        plt.show()

        plt.imshow(ipw_tcf - tcf)
        plt.title('Terrain configuration factor')
        plt.colorbar()
        plt.show()

    def test_skew(self):
        """ Test the skew of an image """

        arr = np.tile(range(1000), (1000, 1))

        for angle in range(-45, 45, 5):
            # skew the initial array
            sarr = viewf.skew(arr, angle=angle)

            # skew it back to original
            sbarr = viewf.skew(sarr, angle=angle, fwd=False)

            # ensure that the skewed back is the same as the original
            self.assertTrue(np.all(arr == sbarr))

        # test the error
        self.assertRaises(ValueError, viewf.skew, arr, -100)
        self.assertRaises(ValueError, viewf.skew, arr, 100)
