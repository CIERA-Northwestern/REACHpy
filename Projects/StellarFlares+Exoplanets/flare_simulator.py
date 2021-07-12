
from astropy import table, constants as const, units as u
import numpy as np
import os
import mpmath

# Abbbreviations:
# eqd = equivalent duration
# ks = 1000 s (obvious perhaps :), but not a common unit)

#region defaults and constants
# some constants
h, c, k_B = const.h, const.c, const.k_B
default_flarespec_path = os.path.join(os.path.dirname(__file__), 'relative_energy_budget.ecsv')
default_flarespec = table.Table.read(default_flarespec_path, format='ascii.ecsv')
default_flarespec = default_flarespec.filled(0)
fuv = [912., 1700.] * u.AA
nuv = [1700., 3200.] * u.AA
version = '1.0'

# the default function for estimating flare peak flux
@u.quantity_input(eqd=u.s)
def boxcar_height_function_default(eqd):
    eqd_s = eqd.to('s').value
    return 0.3*eqd_s**0.6

# other flare defaults
flare_defaults = dict(eqd_min = 100.*u.s,
                      eqd_max = 1e6*u.s,
                      ks_rate = 8/u.d, # rate of ks flares for Si IV (Fig 6 of Loyd+ 2018)
                      cumulative_index = 0.75, # power law index of FUV flares for all stars (Table 5 of Loyd+ 2018)
                      boxcar_height_function = boxcar_height_function_default,
                      decay_boxcar_ratio = 1./2.,
                      BB_SiIV_Eratio=160,  # Hawley et al. 2003
                      T_BB = 9000*u.K,  # Hawley et al. 2003
                      clip_BB = True,
                      SiIV_quiescent=0.1*u.Unit('erg s-1 cm-2'), # for GJ 832 with bolometric flux equal to Earth
                      SiIV_normed_flare_spec=default_flarespec)
#endregion


#region boilerplate code
def _kw_or_default(kws, keys):
    """Boilerplate for pulling from the default dictionary if a desired key isn't present."""
    values = []
    for key in keys:
        if key not in kws or kws[key] is None:
            kws[key] = flare_defaults[key]
        values.append(kws[key])
    return values


def _check_unit(func, var, unit):
    """Boilerplate for checking units of a variable."""
    try:
        var.to(unit)
    except (AttributeError, u.UnitConversionError):
        raise ValueError('Variable {} supplied to the {} must be an '
                         'astropy.Units.Quantity object with units '
                         'convertable to {}'.format(var, func, unit))


def _integrate_spec_table(spec_table):
    """Integrate a spectrum defined in a table with 'w0', 'w1', and 'Edensity' columns."""
    return np.sum((spec_table['w1'] - spec_table['w0']) * spec_table['Edensity'])
#endregion code


#region documentation tools
# there is a lot of duplicated documetation here, so to make sure it is consistent I am going to define it in only one
# place and then insert it into the docstrings, at the cost of readability when actually looking at the source. Sorry
# about that. However, pulling up help on each function should work well, and, like I said, it's more consistent.
_fd = flare_defaults
_flare_params_doc = "flare_params : dictionary\n" \
                    "        Parameters of the flare model. If a parameter is not sepcified, \n" \
                    "        the default is taken from the flare_simulator.flare_defaults \n" \
                    "        dictionary. Parameters relevant to this function are:"
_param_doc_dic = dict(eqd_min = "eqd_min : astropy quantity, units of time\n"
                                "    Minimum flare equivalent duration to be considered.\n"
                                "    Default is {}."
                                "".format(_fd['eqd_min']),
                      eqd_max = "eqd_max : astropy quantity, units of time\n"
                                "    Maxium flare equivalent duration to be considered. \n"
                                "    Default is {}."
                                "".format(_fd['eqd_max']),
                      ks_rate = "ks_rate : astropy quantity, units of time-1\n"
                                "    Rate of Si IV flares with an equivalent duration of 1000 s. \n"
                                "    Default is {}."
                                "".format(_fd['ks_rate']),
                      cumulative_index= "cumulative_index : float\n"
                                        "    Cumulative index of a power-law relating the frequency of flares\n"
                                        "    greater than a given energy to that energy. Default is {}."
                                        "".format(_fd['cumulative_index']),
                      boxcar_height_function = "boxcar_height_function : function\n"
                                               "    Function relating the peak flare flux (height of the boxcar \n"
                                               "    portion of the boxcar-decay model) to the equivalent duration \n"
                                               "    of the flare. The function must accept an equivalent duration \n"
                                               "    as an astropy quantity with units of time as its only input. \n"
                                               "    Default is the function height = 0.3 * equivalent_duration**0.6",
                      decay_boxcar_ratio = "decay_boxcar_ratio : float\n"
                                           "    Ratio between the the amount of flare energy contained in \n"
                                           "    the boxcar portion of the boxcar-decay model and the decay \n"
                                           "    portion. This actually determines the time-constant of the \n"
                                           "    decay. I'm not sure if I actually like that... Default is {}."
                                           "".format(_fd['decay_boxcar_ratio']),
                      BB_SiIV_Eratio = "BB_SiIV_Eratio : float\n"
                                       "    Ratio of the blackbody energy to the Si IV energy of the flare.\n"
                                       "    Default is {}.".format(_fd['BB_SiIV_Eratio']),
                      T_BB = "T_BB : astropy quantity, units of temperature\n"
                             "    Temperature of the flare blackbody continuum. \n"
                             "    Default is {}.".format(_fd['T_BB']),
                      SiIV_quiescent = "SiIV_quiescent : astropy quantity, units of energy time-1 length-2\n"
                                       "    Quiescent flux of the star in the Si IV 1393,1402 AA lines. \n"
                                       "    Default is representative of an inactive M dwarf at the distance \n"
                                       "    where the bolometric irradiation equals that of Earth,\n"
                                       "     {}.".format(_fd['SiIV_quiescent']),
                      SiIV_normed_flare_spec = "SiIV_normed_flare_spec : astropy table\n"
                                               "    Spectral energy budget of the flare (excluding the blackbody) \n"
                                               "    normalized to the combined flux of the Si IV 1393,1402 AA lines. \n"
                                               "    The energy budget  should be an astropy table with columns of\n"
                                               "        'w0' : start of each spectral bin, units of length\n"
                                               "        'w1' : end of each spectral bin, units of length\n"
                                               "        'Edensity' : energy emitted by that flare in the spectral\n"
                                               "                     bin divided by the width of the bin, units of \n"
                                               "                     energy length-1\n"
                                               "    Default is loaded from the 'relative_energy_budget.ecsv' file.",
                      clip_BB = "clip_BB : True|False\n"
                                "    If True (default), do not include blackbody flux in the FUV range \n"
                                "    and shortward. This is done because BB flux is presumed to be \n"
                                "    included in the flare SED at EUV and FUV wavelengths assembled by \n"
                                "    Loyd+ 2018 that is the default here. However, should be changed to\n"
                                "    False if, e.g., a hotter or more energetic blackbody is adopted.")
_tbins_doc = 'tbins : astropy quantity array, units of time\n' \
             '        Edges of the lightcurve time bins.'
_wbins_doc = 'wbins : astropy quantity array, units of length\n' \
             '        Edges of the spectral bins.'
_t0_doc = 't0 : astropy quantity, units of time\n' \
          '        Start time of flare.'
_eqd_doc = 'eqd : astropy quantity, units of time\n' \
           '        Equivalent duration of flare in the Si IV 1393,1402 line \n' \
           '        (flare energy divided by star\'s quiescent luminosity\n' \
           '        in the same band).'

def add_indent(txt):
    return "    " + txt.replace('\n', '\n    ')
def _get_param_string(*keys):
    strings = [_param_doc_dic[key] for key in keys]
    strings = list(map(add_indent, strings))
    strings = list(map(add_indent, strings))
    return '\n'.join([_flare_params_doc] + strings)
def _format_doc(func, **kws):
    func.__doc__ = func.__doc__.format(**kws)
#endregion


#region fast planck function computations
_Li = mpmath.fp.polylog
def _P3(x):
    """Dang, I should have cited where I got this. Now it is lost."""
    e = np.exp(-x)
    return _Li(4, e) + x*_Li(3, e) + x**2/2*_Li(2, e) + x**3/6*_Li(1, e)
_P3 = np.vectorize(_P3)

@u.quantity_input(w=u.AA, T=u.K)
def _blackbody_partial_integral(w, T):
    """
    Integral of blackbody surface flux at wavelengths from 0 to w.

    Parameters
    ----------
    w : astropy quantity, units of length
        wavelength to which to integrate
    T : astropy quantity, units of temperature
        temperature of blackbody

    Returns
    -------
    I : astropy quantity
    """
    x = (h*c/w/k_B/T).to('').value
    I = 12 * np.pi * (k_B*T)**4 / c**2 / h**3 * _P3(x)
    return I.to('erg s-1 cm-2')


@u.quantity_input(wbins=u.AA, T=u.K)
def blackbody_binned(wbins, T, bolometric=None):
    """
    Quick computation of blackbody surface flux integrated within wbins.

    This is especially helpful if there are large wavelength bins where taking the value of the Planck function at the
    midpoint might give inaccurate results.

    Parameters
    ----------
    {wbins}
    T : astropy quantity, units of temperature
        temperature of blackbody
    bolometric : astropy quantity, units of energy time-1 length-2
        value of the bolometric blackbody flux by which to normalize the
        output.  A value of None gives the flux at the surface of the
        emitter.

    Returns
    -------
    flux_density : astropy quantity, units of energy time-1 length-3
        The flux spectral density of the blackbody in each wbin, generally in units of erg s-1 cm-2 AA-1.
    """

    # take difference of cumulative integral at each bin edge to get flux in each bin
    F = np.diff(_blackbody_partial_integral(wbins, T))

    # divide by bin widths to get flux density
    f = F / np.diff(wbins)

    # renormalize, if desired, and return
    if bolometric is None:
        return f.to('erg s-1 cm-2 AA-1')
    else:
        fbolo = const.sigma_sb*T**4
        fnorm = (f/fbolo).to(1/wbins.unit)
        return fnorm*bolometric
_format_doc(blackbody_binned, wbins=_wbins_doc)


@u.quantity_input(wbins=u.AA, T=u.K)
def blackbody_points(w, T, bolometric=None):
    """
    Compute the flux spectral density of the emission from a blackbody.

    Returns the value at each w, rather than the value averaged over wbins. For the latter, use blackbody_binned.

    Parameters
    ----------
    w : astropy quantity array, units of length
        Wavelengths at which to compute flux density.
    T : astropy quantity, units of temperature
        temperature of blackbody
    bolometric : astropy quantity, units of energy time-1 length-2
        value of the bolometric blackbody flux by which to normalize the
        output.  A value of None gives the flux at the surface of the
        emitter.

    Returns
    -------
    flux_density : astropy quantity, units of energy time-1 length-3
        The flux spectral density of the blackbody at each w, generally in units of erg s-1 cm-2 AA-1.
    """
    # compute flux density from Planck function (with that extra pi factor to get rid of per unit solid angle portion)
    f = np.pi * 2 * const.h * const.c ** 2 / w ** 5 / (np.exp(const.h * const.c / const.k_B / T / w) - 1)

    # return flux density, renormalized if desired
    if bolometric is None:
        return f.to('erg s-1 cm-2 AA-1')
    else:
        fbolo = const.sigma_sb*T**4
        fnorm = (f/fbolo).to(1/w.unit)
        return fnorm*bolometric
#endregion


#region utilities
def rebin(bins_new, bins_old, y):
    """
    Rebin some binned values.

    Parameters
    ----------
    bins_new : array
        New bin edges.
    bins_old : array
        Old bin edges.
    y : array
        Binned values (average of some function like a spectrum across
        each bin).

    Returns
    -------
    y_new : array
        Rebinned values.

    """
    # politely let user no that quantity input is not desired for this
    if any(isinstance(x, u.Quantity) for x in [bins_new, bins_old, y]):
        raise ValueError('No astropy Quantity input for this function, please.')
    if np.any(bins_old[1:] <= bins_old[:-1]) or np.any(bins_new[1:] <= bins_new[:-1]):
        raise ValueError('Old and new bin edges must be monotonically increasing.')

    # compute cumulative integral of binned data
    areas = y*np.diff(bins_old)
    I = np.cumsum(areas)
    I = np.insert(I, 0, 0)

    # compute average value in new bins
    Iedges = np.interp(bins_new, bins_old, I)
    y_new = np.diff(Iedges)/np.diff(bins_new)

    return y_new


def power_rv(min, max, cumulative_index, n):
    """
    Random values drawn from a power-law distribution.

    Parameters
    ----------
    min : float
        Minimum value of the distribution.
    max : float
        Maximum value of the distribution.
    cumulative_index : float
        Index of the cumulative distribution.
    n : integer
        Number of values to draw.

    Returns
    -------
    values : array
        Array of random values.
    """

    # politely let user know that, in this instance, astropy Quantities are not wanted
    if any(isinstance(x, u.Quantity) for x in [min, max, cumulative_index]):
        raise ValueError('No astropy Quantity input for this function, please.')

    # I found it easier to just make my own than figure out the numpy power, pareto, etc. random number generators
    a = cumulative_index
    norm = min**-a - max**-a
    # cdf = 1 - ((x**-a - max**-a)/norm)
    x_from_cdf = lambda c: ((1-c)*norm + max**-a)**(-1/a)
    x_uniform = np.random.uniform(size=n)
    return x_from_cdf(x_uniform)


def shot_times(rate, time_span):
    """
    Generate random times of events that when binned into even intervals would yield counts that are Poisson distributed.

    Parameters
    ----------
    rate : float
        Average rate of events.
    time_span : float
        Length of time over which to generate events.

    Returns
    -------
    times : array
        Times at which random events occurr.
    """

    # politely let user know that, in this instance, astropy Quantities are not wanted
    if any(isinstance(x, u.Quantity) for x in [rate, time_span]):
        raise ValueError('No astropy Quantity input for this function, please.')
    # generate wait times from exponential distribution (for poisson stats)
    # attempt drawing 10 std devs more "shots" than the number expected to fill time_span so chances are very low it
    # won't be filled
    avg_wait_time = 1. / rate
    navg = time_span / avg_wait_time
    ndraw = int(navg + 10*np.sqrt(navg))
    wait_times = np.random.exponential(avg_wait_time, size=ndraw)

    # cumulatively sum wait_times to get actual event times
    tshot = np.cumsum(wait_times)

    # if the last event occurs before user-specified length of time, try again. Else, return the times.
    if tshot[-1] < time_span:
        return shot_times(rate, time_span)
    else:
        return tshot[tshot < time_span]


def boxcar_decay(tbins, t0, area_box, height_box, area_decay):
    """
    Compute the lightcurve from one or more boxcar-decay functions.

    Parameters
    ----------
    tbins : array
        edges of the time bins used for the lightcurve
    t0 : float or array
        start times of the boxcar-decays
    area_box : float or array
        areas of the boxcar portion of the boxcar-decays
    height_box : float or array
        heights of the boxcar-decays
    area_decay : float or array
        areas of the decay portions of the boxcar-decays

    Returns
    -------
    y : array
        lightcurve values

    Notes
    -----
    This function is a bottleneck when creating a lightcurve from a long
    series of flares. If this code is to be adapted for quick simulation
    of years-long series of flares, this is where the speedup needs to
    happen.
    """

    # politely let user know that, in this instance, astropy Quantities are not wanted
    if any(isinstance(x, u.Quantity) for x in [tbins, t0, area_box, height_box, area_decay]):
        raise ValueError('No astropy Quantity input for this function, please.')

    # this is going to have to be ugly for it to be fast, I think

    # standardize t0, area_box, height_box, and area_decay for array input
    t0, area_box, height_box, area_decay = [np.reshape(a, [-1]) for a in [t0, area_box, height_box, area_decay]]

    # compute end of box, start of decay
    t1 = t0 + area_box/height_box

    # correct for portions hanging over ends of tbins
    t0 = np.copy(t0)
    t0[t0 < tbins[0]] = tbins[0]
    t1[t1 > tbins[-1]] = tbins[-1]

    # initialize y array
    y = np.zeros((len(t0), len(tbins)-1))
    i_rows = np.arange(y.shape[0])

    # add starting portion of box to first bin that is only partially covered by it
    i0 = np.searchsorted(tbins, t0, side='right')
    frac = (tbins[i0] - t0)/(tbins[i0] - tbins[i0-1])
    y[i_rows, i0-1] += frac*height_box

    # add box to bins fully covered by it
    inbox = (tbins[None, :-1] > t0[:, None]) & (tbins[None, 1:] < t1[:, None])
    y += height_box[:,None]*inbox

    # add ending fraction of box to last bin that is partially covered by it
    i1 = np.searchsorted(tbins, t1, side='left')
    frac = (t1 - tbins[i1-1])/(tbins[i1] - tbins[i1-1])
    y[i_rows, i1-1] += frac*height_box

    # deal with any cases where the box was entirely within a bin
    j = i0 == i1
    y[i_rows[j], i0[j]-1] = area_box[j]/(tbins[i0][j] - tbins[i0-1][j])

    # add decay
    # compute cumulative decay integral at all time points
    amp_decay = height_box
    tau_decay = area_decay / amp_decay
    with np.errstate(over='ignore', invalid='ignore'):
        Idecay = -amp_decay[:,None]*tau_decay[:,None]*np.exp(-(tbins[None,:] - t1[:,None])/tau_decay[:,None])
        ydecay = np.diff(Idecay, 1)/np.diff(tbins)
    keep = tbins[:-1] > t1[:, None]
    y[keep] += ydecay[keep]

    # add fractional piece of exponential
    i1 = np.searchsorted(tbins, t1, side='right')
    inrange = i1 < len(tbins)
    i_rows, i1 = i_rows[inrange], i1[inrange]
    Idecay1 = -amp_decay*tau_decay
    ydecay1 = (Idecay[i_rows, i1] - Idecay1[i_rows])/(tbins[i1] - tbins[i1-1])
    y[i_rows, i1-1] += ydecay1

    return np.sum(y, 0)
#endregion


#region front end functions
@u.quantity_input(eqd=u.AA, filter_response=u.Unit(''))
def filter_to_SiIV_energy(filter_wave, filter_response, energy, **flare_params):
    """
    Convenience function for converting the energy in a photometric filter to the Si IV energy of a flare.

    Parameters
    ----------
    filter_wave : astropy quantity array, units of length
        Wavelengths of filter response curve.
    filter_response : array, unitless
        Filter response at filter_wave.
    energy : float or astropy quantity, units of energy
        Energy of the flare in the specified filter.
    {flare_params}

    Returns
    -------
    energy_SiIV : float or astropy quantity
        Energy of the flare in the Si IV 1393,1402 AA line.
    """

    # get filter-convolved fraction of flare energy relative to Si IV
    w_mids = (filter_wave[1:] + filter_wave[:-1])/2.
    w_bins = np.insert(w_mids.value, [0,len(w_mids)],
                       filter_wave[[0,-1]].value)*filter_wave.unit
    flux = flare_spectrum(w_bins, 1.0, **flare_params)
    filter_fraction = np.sum(filter_response*flux*np.diff(w_bins))

    # then just invert to get the energy in Si IV given the filter energy
    return energy/filter_fraction
_format_doc(filter_to_SiIV_energy, flare_params=_get_param_string('BB_SiIV_Eratio', 'T_BB', 'SiIV_normed_flare_spec'))


@u.quantity_input(tbins=u.s, t0=u.s, eqd=u.s)
def flare_lightcurve(tbins, t0, eqd, **flare_params):
    """
    Return a lightcurve for a single flare normalized to quiescent flux.

    Parameters
    ----------
    {tbins}
    {t0}
    {eqd}
    {flare_params}

    Returns
    -------
    y : array
        Quiescent-normalized lightcurve of the flare.
    """

    # get relevant flare parameters
    values = _kw_or_default(flare_params, ['boxcar_height_function', 'decay_boxcar_ratio'])
    boxcar_height_function, decay_boxcar_ratio = values

    # compute boxcar parameters
    boxcar_height = boxcar_height_function(eqd)
    boxcar_area = eqd/(1 + decay_boxcar_ratio)
    decay_area = boxcar_area * decay_boxcar_ratio

    # make units uniform
    tunit = tbins.unit
    tbins, t0, eqd, boxcar_area, decay_area = [x.to(tunit).value for x in [tbins, t0, eqd, boxcar_area, decay_area]]
    y = boxcar_decay(tbins, t0, boxcar_area, boxcar_height, decay_area)

    return y
_format_doc(flare_lightcurve, flare_params=_get_param_string('boxcar_height_function', 'decay_boxcar_ratio'),
            tbins=_tbins_doc, eqd=_eqd_doc, t0=_t0_doc)


def flare_rate(**flare_params):
    """
    Rate of flares spanning the given energy range.

    Parameters
    ----------
    {flare_params}

    Returns
    -------
    rate : astropy quantity
    """
    # get relevant flare parameters and check units
    values = _kw_or_default(flare_params, ['eqd_min', 'eqd_max', 'ks_rate', 'cumulative_index'])
    eqd_min, eqd_max, ks_rate, cumulative_index = values
    _check_unit(flare_rate, ks_rate, 's-1')
    [_check_unit(flare_rate, v, 's') for v in [eqd_min, eqd_max]]

    # make sure no stupid input
    if eqd_min <= 0:
        raise ValueError('Flare rate diverges at eqd_min == 0. Only eqd_min > 0 makes sense.')

    # compute rate
    rate = ks_rate * ((eqd_min/u.ks).to('')**-cumulative_index - (eqd_max/u.ks).to('')**-cumulative_index)

    return rate.to('d-1')
_format_doc(flare_rate, flare_params=_get_param_string('eqd_min', 'eqd_max', 'ks_rate', 'cumulative_index'))


@u.quantity_input(time_span=u.s)
def flare_series(time_span, **flare_params):
    """
    Start times and equivalent durations for a randomly generated series of flares.

    Parameters
    ----------
    time_span : astropy quantity, units of time
    {flare_params}

    Returns
    -------
    t_flare : astropy quantity array, units of time
        Start times of the random flares.
    eqd : astropy quantity array, units of time
        Equivalent durations of the random flares.
    """
    values = _kw_or_default(flare_params, ['eqd_min', 'eqd_max', 'cumulative_index'])
    eqd_min, eqd_max, cumulative_index = values
    [_check_unit(flare_series, v, 's') for v in [eqd_min, eqd_max]]

    # get the expected flare rate
    rate = flare_rate(**flare_params)

    # draw flares at that rate
    tunit = time_span.unit
    rate = rate.to(tunit**-1).value
    time_span = time_span.value
    t_flare = shot_times(rate, time_span) * tunit
    n = len(t_flare)

    # draw energies for those flares
    eqd_min, eqd_max = [x.to(tunit).value for x in [eqd_min, eqd_max]]
    eqd = power_rv(eqd_min, eqd_max, cumulative_index, n) * tunit

    return t_flare, eqd
_format_doc(flare_series, flare_params=_get_param_string('eqd_min', 'eqd_max', 'cumulative_index'))


@u.quantity_input(tbins=u.s)
def flare_series_lightcurve(tbins, return_flares=False, **flare_params):
    """
    Generate a series of random flares and return their lightcurve.

    Parameters
    ----------
    {tbins}
    return_flares : bool
        If True, return the start times and equivalent durations of
        the flares.
    {flare_params}

    Returns
    -------
    y : array
        Quiescent-normalized ightcurve values in each tbin.
    tflares : astropy quantity array, units of time, optional
        Start time of random flares.
    eqd : astropy quantity array, units of time, optional
        Equivalent durations of the random flares.
    """
    # generate random flares
    time_span = tbins[-1] - tbins[0]
    tflares, eqds = flare_series(time_span, **flare_params)

    # make lightcurve from those flares
    y = flare_lightcurve(tbins, tflares, eqds, **flare_params)

    if return_flares:
        return y, (tflares, eqds)
    else:
        return y
_format_doc(flare_series_lightcurve, tbins=_tbins_doc,
            flare_params=_get_param_string('eqd_min', 'eqd_max', 'cumulative_index', 'boxcar_height_function',
                                           'decay_boxcar_ratio'))


@u.quantity_input(wbins=u.AA)
def flare_spectrum(wbins, SiIV, **flare_params):
    """
    Return the flare spectrum scaled to match the energy or equivalent duration specified by SiIV and binned according
    to wbins.

    Parameters
    ----------
    {wbins}
    SiIV : float or astropy quantity
        Equivalent duration or energy of the flare in the Si IV 1393,1402 AA
        line. This could also be peak flux or some other quantity, but note
        that you should probably specificy your own 'Edensity' column of  the
        SiIV_normed_flare_spec table to match if so.
    {flare_params}

    Returns
    -------
    spectrum : astropy quantity array, units variabile according to units of SiIV
        Energy spectral density or other spectral density of the flare spectrum
        in each wbin.
    """
    BBratio, T, flarespec, clip_BB  = _kw_or_default(flare_params, ['BB_SiIV_Eratio', 'T_BB', 'SiIV_normed_flare_spec',
                                                            'clip_BB'])

    # rebin energy density from specified flare SED (from MUSCLES data by default)
    fs_bins = np.append(flarespec['w0'], flarespec['w1'][-1]) * flarespec['w0'].unit
    fs_density = flarespec['Edensity'].quantity
    fs_bins = fs_bins.to(wbins.unit)
    FUV_and_lines = rebin(wbins.value, fs_bins.value, fs_density.value) * fs_density.unit * SiIV

    # get the blackbody (should not be included in SED) emission in each bin. Add to regions shortward of FUV as
    # desired by user
    BBbolo = BBratio * SiIV
    if clip_BB:
        red = (wbins[1:] > fuv[1])
        BBbins = np.insert(wbins[1:][red], 0, fuv[1])
        BB = blackbody_binned(BBbins, T, bolometric=BBbolo)

        # add SED and blackbody
        result = FUV_and_lines
        result[red] += BB
    else:
        BB = blackbody_binned(wbins, T, bolometric=BBbolo)
        result = FUV_and_lines + BB

    return result
_format_doc(flare_spectrum, wbins=_wbins_doc,
            flare_params=_get_param_string('BB_SiIV_Eratio', 'T_BB', 'SiIV_normed_flare_spec'))


@u.quantity_input(wbins=u.AA, tbins=u.s, t0=u.s, eqd=u.s)
def flare_spectra(wbins, tbins, t0, eqd, **flare_params):
    """
    Return a series of flare spectra averaged over each tbin for a flare starting at t0 with equivalent duration eqd
    in Si IV.

    Parameters
    ----------
    {wbins}
    {tbins}
    {t0}
    {eqd}
    {flare_params}

    Returns
    -------
    spectra : astropy quantity array, variable units
        Array of spectra in each tbin, where the array has dimensions
        (len(tbins)-1, len(wbins)-1). Units will match the product
        of the eqd and SiIV_quiescent units, divided by time and length.
    """
    # get quiescent Si IV flux
    SiIVq, = _kw_or_default(flare_params, ['SiIV_quiescent'])

    # get lightcurve of flare
    lightcurve = flare_lightcurve(tbins, t0, eqd, **flare_params)

    # get spectrum of flare
    spectrum = flare_spectrum(wbins, SiIVq, **flare_params)

    # multiply spectrum by (quiecent-normalized) lightcurve to get array of spectra in each tbin
    return np.outer(lightcurve, spectrum.value)*spectrum.unit
_format_doc(flare_spectra, wbins=_wbins_doc, tbins=_tbins_doc, t0=_t0_doc, eqd=_eqd_doc,
            flare_params=_get_param_string('SiIV_quiescent', 'boxcar_height_function', 'decay_boxcar_ratio',
                                           'BB_SiIV_Eratio', 'T_BB', 'SiIV_normed_flare_spec'))


@u.quantity_input(wbins=u.AA, tbins=u.s)
def flare_series_spectra(wbins, tbins, **flare_params):
    """
    Generate time-evolving spectra from a random series of flares.

    Parameters
    ----------
    {wbins}
    {tbins}
    {flare_params}

    Returns
    -------
    spectra : astropy quantity array
        Array of spectra in each tbin, where the array has dimensions
        (len(tbins)-1, len(wbins)-1). Units will match the product of
        the eqd and SiIV_quiescent units, divided by time and length.
    """
    # get quiescent Si IV flux
    SiIVq, = _kw_or_default(flare_params, ['SiIV_quiescent'])

    # get lightcurve of a series of random flares
    lightcurve = flare_series_lightcurve(tbins, **flare_params)

    # get spectrum of flares
    spectrum = flare_spectrum(wbins, SiIVq, **flare_params)

    # multiply spectrum by (quiecent-normalized) lightcurve to get array of spectra in each tbin
    return np.outer(lightcurve, spectrum.value)*spectrum.unit
_format_doc(flare_series_spectra, wbins=_wbins_doc, tbins=_tbins_doc,
            flare_params=_get_param_string('SiIV_quiescent', 'eqd_min', 'eqd_max', 'cumulative_index',
                                           'boxcar_height_function', 'decay_boxcar_ratio', 'BB_SiIV_Eratio', 'T_BB',
                                           'SiIV_normed_flare_spec'))
#endregion
