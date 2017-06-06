"""
Created on March 14, 2017
Originally written by Scott Havens in 2015
@author: Micah Johnson

**Custom NASDE Models**:
    When creating a new NASDE model make sure you adhere to the following:

    1. Add a new method with the other models with a unique name.

    2. Add your method to the available_models dictionary using the format of the
       expected config file string followed by your function name.

    3. Update documentation and run smrf!

"""

__version__ = '0.2.1'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from distutils.command.check import check


def calc_phase_and_density(temperature, precipitation, nasde_model):
    """
    Uses various new accumulated snow density models to estimate the snow
    density of precipitation that falls during sub-zero conditions.
    The models all are based on the dew point temperature and the amount of
    precipitation, All models used here must return a dictionary containing the
    keywords pcs and rho_s for percent snow and snow density respectively.

    Args:
        temperature: a single timestep of the distributed dew point temperature

        precipitation: a numpy array of the distributed precipitation

        nasde_model: string value set in the configuration file representing the method
                    for estimating density of new snow that has just fallen.

                    Available Options:
                        * susong1999, see :func:`~smrf.envphys.snow.susong1999`
                        * piecewise_susong1999, see :func:`~smrf.envphys.snow.piecewise_susong1999`
                        * marks2017, see :func:`~smrf.envphys.snow.marks2017`

    Returns:
        tuple:
           - **snow_density** (*numpy.array*) - Snow density values in kg/m^3
            estimated by the selelected nasde_model.

           - **perc_snow** (*numpy.array*) - Percent of the precip that is snow
            in values 0.0-1.0, which are estimated by the selelected nasde_model.

    """
    # convert the inputs to numpy arrays
    precip = np.array(precipitation)
    temperature = np.array(temperature)

    snow_density = np.zeros(precip.shape)
    perc_snow = np.zeros(precip.shape)
    #New accumulated snow point models can go here.

    if nasde_model in available_models.keys():
        result = available_models[nasde_model](temperature,precip)
    else:
        raise ValueError("{0} is not an implemented new accumulated snow density model (NASDE)! Check the config file under precipitation".format(nasde_model))

    snow_density = result['rho_s']
    perc_snow = result['pcs']

    return snow_density,perc_snow


def calc_perc_snow(Tpp, Tmax=0.0, Tmin=-10.0):
    """
    Calculates the percent snow for the nasde_models piecewise_susong1999 and marks2017.

    Args:
        Tpp: A numpy array of temperature, use dew point temperature
        if available [degree C].

        Tmax: Max temperature that the percent snow is estimated. Default is 0.0 Degrees C.

        Tmin: Minimum temperature that percent snow is changed. Default is -10.0 Degrees C.

    Returns:
        numpy.array:
        A fraction of the precip at each pixel that is snow provided by Tpp.
    """

    #Coefficients for snow relationship
    Tr0 = 0.5
    Pcr0 = 0.25
    Pc0 = 0.75

    # this shouldn't have to be done again unless it's outside this module
    # We can set a flag later if it should be checked
    #Set a cap on temperature
#     Tpp, tsnow = check_temperature(Tpp, Tmax = Tmax, Tmin = Tmin)


    pcs = np.zeros(Tpp.shape)

    pcs[Tpp <= -0.5] = 1.0

    ind = (Tpp > -0.5) & (Tpp <= 0.0)
    if np.any(ind):
        pcs[ind] = (-Tpp[ind] / Tr0) * Pcr0 + Pc0

    ind = (Tpp > 0.0) & (Tpp <= Tmax +1.0)
    if np.any(ind):
        pcs[ind] = (-Tpp[ind] / (Tmax + 1.0)) * Pc0 + Pc0

    return pcs


def check_temperature(Tpp, Tmax = 0.0, Tmin = -10.0):
    """
    Sets  the precipitation temperature and snow temperature.

    Args:
        Tpp: A numpy array of temperature, use dew point temperature
        if available [degrees C].

        Tmax: Caps the temperature of the snow [degrees C].

        Tmin: Minimum temperature that the precipitation temperature [degrees C].

    Returns:
        tuple:
           - **Tpp** (*numpy.array*) - Modified precipitation temperature that
            is thresholded with a minimum set by tmin.

           - **tsnow** (*numpy.array*) - Temperature of the surface of the snow
            set by the precipitation temperature and thresholded by tmax where
            tsnow > tmax = tmax.
    """

    Tpp[Tpp < Tmin] = Tmin

    tsnow = Tpp.copy()
    tsnow[Tpp > Tmax] = Tmax

    return Tpp, tsnow


#BEGIN NASDE MODELS HERE AND BELOW

def susong1999(temperature, precipitation):
    """
    Follows the IPW command mkprecip

    The precipitation phase, or the amount of precipitation falling as rain or snow, can significantly
    alter the energy and mass balance of the snowpack, either leading to snow accumulation or inducing
    melt :cite:`Marks&al:1998` :cite:`Kormos&al:2014`. The precipitation phase and initial snow density are
    based on the precipitation temperature (the distributed dew point temperature) and are estimated
    after Susong et al (1999) :cite:`Susong&al:1999`. The table below shows the relationship to
    precipitation temperature:

    ========= ======== ============ ===============
    Min Temp  Max Temp Percent snow Snow density
    [deg C]   [deg C]  [%]          [kg/m^3]
    ========= ======== ============ ===============
    -Inf      -5       100          75
    -5        -3       100          100
    -3        -1.5     100          150
    -1.5      -0.5     100          175
    -0.5      0        75           200
    0         0.5      25           250
    0.5       Inf      0            0
    ========= ======== ============ ===============

    Args:
    precipitation - numpy array of precipitation values [mm]
    temperature - array of temperature values, use dew point temperature
        if available [degrees C]

    Returns:
        dictionary:
           - **perc_snow** (*numpy.array*) - Percent of the precipitation that is snow
            in values 0.0-1.0.

           - **rho_s** (*numpy.array*) - Snow density values in kg/m^3.

    """

    # create a list from the table above
    t = []
    t.append( {'temp_min': -1e309,  'temp_max': -5,     'snow': 1,    'density':75} )
    t.append( {'temp_min': -5,      'temp_max': -3,     'snow': 1,    'density':100} )
    t.append( {'temp_min': -3,      'temp_max': -1.5,    'snow': 1,    'density':150} )
    t.append( {'temp_min': -1.5,    'temp_max': -0.5,   'snow': 1,    'density':175} )
    t.append( {'temp_min': -0.5,    'temp_max': 0.0,    'snow': 0.75, 'density':200} )
    t.append( {'temp_min': 0.0,     'temp_max': 0.5,    'snow': 0.25, 'density':250} )
    t.append( {'temp_min': 0.5,     'temp_max': 1e309,  'snow': 0,    'density':0} )


    # preallocate the percent snow (ps) and snow density (sd)
    ps = np.zeros(precipitation.shape)
    sd = np.zeros(ps.shape)

    # if no precipitation return all zeros
    if np.sum(precipitation) == 0:
        return ps, sd

    # determine the indicies and allocate based on the table above
    for row in t:

        # get the values between the temperature ranges that have precip
        ind = [(temperature >= row['temp_min']) & (temperature < row['temp_max'])]

        # set the percent snow
        ps[ind] = row['snow']

        # set the density
        sd[ind] = row['density']


    # if there is no precipitation at a pixel, don't report a value
    # this may make isnobal crash, I'm not really sure
    ps[precipitation == 0] = 0
    sd[precipitation == 0] = 0

    return {'pcs':ps, 'rho_s':sd}


def piecewise_susong1999(Tpp, precip, Tmax = 0.0, Tmin = -10.0, check_temps=True):
    """
    Follows :func:`~smrf.envphys.snow.piecewise_susong1999` but is the piecewise form of table shown there.
    This model adds to the former by accounting for liquid water effect near 0.0 Degrees C.

    The table was estimated by Danny Marks in 2017 which resulted in the
    piecewise equations below:

    Percent Snow:

    .. math::

        \\%_{snow} = \\Bigg \\lbrace{
            \\frac{-T}{T_{r0}} P_{cr0} + P_{c0}, \\hfill -0.5^{\\circ} C \\leq T \\leq 0.0^{\\circ} C
             \\atop
             \\frac{-T_{pp}}{T_{max} + 1.0} P_{c0} + P_{c0}, \\hfill 0.0^\\circ C \\leq T \\leq T_{max}
            }

    Snow Density:
        .. math::
            \\rho_{s} = 50.0 + 1.7 * (T_{pp} + 15.0)^{ex}


            ex = \\Bigg \\lbrace{
                ex_min + \\frac{T_{range} + T_{snow} - T_{max}}{T_{range}} * ex_{r}, \\hfill ex < 1.75
                \\atop
                1.75, \\hfill, ex > 1.75
                }
    Args:
        Tpp: A numpy array of temperature, use dew point temperature
        if available [degree C].

        precip: A numpy array of precip in millimeters.

        Tmax: Max temperature that snow density is modeled. Default is 0.0 Degrees C.

        Tmin: Minimum temperature that snow density is changing. Default is -10.0 Degrees C.

        check_temps: A boolean value check to apply special temperature constraints,
            this is done using :func:`~smrf.envphys.snow.check_temperature`. Default is True.

    Returns:
        dictionary:
           - **pcs** (*numpy.array*) - Percent of the precipitation that is snow
            in values 0.0-1.0.

           - **rho_s** (*numpy.array*) - Density of the fresh snow in kg/m^3.

    """
    ex_max = 1.75
    exr = 0.75
    ex_min = 1.0

    Tz = 0.0

    # again, this shouldn't be needed in this context
    if check_temps:
        Tpp, tsnow = check_temperature(Tpp, Tmax=Tmax, Tmin=Tmin)

    pcs = calc_perc_snow(Tpp,Tmin=Tmin, Tmax = Tmax)

    # new snow density - no compaction
    Trange = Tmax - Tmin
    ex = ex_min + (((Trange + (tsnow - Tmax)) / Trange) * exr)

    ex[ex > ex_max] = ex_max

    rho_ns = 50.0 + (1.7 * ((Tpp-Tz) + 15.0)**ex)

    return {'pcs':pcs, 'rho_s':rho_ns}


def marks2017(Tpp,pp):

    """
    A new accumulated snow density model that accounts for compaction. The model
    builds upon :func:`~smrf.envphys.snow.piecewise_susong1999` by adding effects
    from compaction. Of four mechanisms for compaction, this model accounts for
    compaction by destructive metmorphosism and overburden. These two processes are
    accounted for by calculating a proportionalility using data from Kojima,
    Yosida and Mellor. The overburden is in part estimated using total storm deposition,
    where storms are defined in :func:`~smrf.envphys.storms.tracking_by_station`.
    Once this is determined the final snow density is applied through the entire storm
    only varying with hourly temperature.

    Snow Density:
        .. math::
            \\rho_{s} = \\rho_{ns} + (\\Delta \\rho_{c} + \\Delta \\rho_{m}) \\rho_{ns}

    Overburden Proportionality:
        .. math::
            \\Delta \\rho_{c} = 0.026 e^{-0.08 (T_{z} - T_{snow})}  SWE*  e^{-21.0 \\rho_{ns}}

    Metmorphosism Proportionality:
        .. math::
            \\Delta \\rho_{m} = 0.01 c11 e^{-0.04 (T_{z} - T_{snow})}

            c11 = c_min + (T_{z} - T_{snow}) C_{factor} + 1.0

    Constants:
                C_{factor} = 0.0013
                Tz = 0.0


    Args:
        Tpp: Numpy array of a single hour of temperature, use dew point if available [degrees C].

        pp: Numpy array representing the total amount of precip deposited during a storm
            in millimeters

    Returns:
        dictionary:

            - **rho_s** (*numpy.array*) - Density of the fresh snow in kg/m^3.

            - **swe** (*numpy.array*) - Snow water equivalent.

            - **pcs** (*numpy.array*) - Percent of the precipitation that is
             snow in values 0.0-1.0.

            - **rho_ns** (*numpy.array*) - Density of the uncompacted snow, which
             is equivalent to the output from
             :func:`~smrf.envphys.snow.piecewise_susong1999`.

            - **d_rho_c** (*numpy.array*) - Prportional coefficient for
            compaction from overburden.

            - **d_rho_m** (*numpy.array*) - Prportional coefficient for
             compaction from melt.

            - **rho_s** (*numpy.array*) - Final density of the snow [kg/m^3].

            - **rho** (*numpy.array*) - Density of the precipitation, which
             continuously ranges from low density snow to pure liquid water
             (50-1000 kg/m^3).

            - **zs** (*numpy.array*) - Snow height added from the precipitation.

    """

    ex_max = 1.75
    exr = 0.75
    ex_min = 1.0
    c1_min = 0.026
    c1_max = 0.069
    c1r = 0.043
    c_min = 0.0067
    cfac = 0.0013
    Tmin = -10.0
    Tmax = 0.0
    Tz = 0.0
    Tr0 = 0.5
    Pcr0 = 0.25
    Pc0 = 0.75

    water = 1000.0

    rho_ns = d_rho_m = d_rho_c = zs = rho_s = swe = pcs = np.zeros(Tpp.shape)
    rho = np.ones(Tpp.shape)

    if np.any(pp > 0):

        # check the precipitation temperature
        Tpp, tsnow = check_temperature(Tpp, Tmax=Tmax, Tmin=Tmin)

        # Calculate the percent snow and initial uncompacted density
        result = piecewise_susong1999(Tpp,pp,Tmax=Tmax, Tmin=Tmin)
        pcs = result['pcs']
        rho_orig = result['rho_s']

        swe = pp * pcs

        # Calculate the density only where there is swe
        swe_ind = swe > 0.0
        if np.any(swe_ind):

            s_pcs = pcs[swe_ind]
            s_pp = pp[swe_ind]
            s_swe = swe[swe_ind]
            s_tpp = Tpp[swe_ind] # transforms to a 1D array, will put back later
            s_tsnow = tsnow[swe_ind] # transforms to a 1D array, will put back later

            s_rho_ns = rho_orig[swe_ind] # transforms to a 1D array, will put back later

            #Convert to a percentage of water
            s_rho_ns /= water

            # proportional total storm mass compaction
            s_d_rho_c = (0.026 * np.exp(-0.08 * (Tz - s_tsnow)) * s_swe * np.exp(-21.0 * s_rho_ns))

            ind = s_rho_ns * water > 100.0
            c11 = np.ones(s_rho_ns.shape)

            c11[ind] = (c_min + ((Tz - s_tsnow[ind]) * cfac)) + 1.0

            s_d_rho_m = 0.01 * c11 * np.exp(-0.04 * (Tz - s_tsnow))

            # compute snow density, depth & combined liquid and snow density
            s_rho_s = s_rho_ns +((s_d_rho_c + s_d_rho_m) * s_rho_ns)

            s_zs = s_swe / s_rho_s

            # a mixture of snow and liquid
            s_rho = s_rho_s.copy()
            mix = (s_swe < s_pp) & (s_pcs > 0)
            if np.any(mix):
                s_rho[mix] = (s_pcs[mix] * s_rho_s[mix]) + (1.0 - s_pcs[mix])

            s_rho[s_rho > 1.0] = 1.0

            # put the values back into the larger array
            rho_ns[swe_ind] = s_rho_ns
            d_rho_m[swe_ind] = s_d_rho_m
            d_rho_c[swe_ind] = s_d_rho_c
            zs[swe_ind] = s_zs
            rho_s[swe_ind]= s_rho_s
            rho[swe_ind] = s_rho
            pcs[swe_ind] = s_pcs

    # convert densities from proportions, to kg/m^3 or mm/m^2
    rho_ns *= water
    rho_s *= water
    rho *= water

    return {'swe':swe, 'pcs':pcs,'rho_ns': rho_ns, 'd_rho_c' : d_rho_c, 'd_rho_m' : d_rho_m, 'rho_s' : rho_s, 'rho':rho, 'zs':zs}

if __name__ == '__main__':
    print("\nNothing implemented here.")
