import os
import os.path as op

import numpy as num

from pyrocko.plot.dynamic_rupture import \
    clear_temp, make_colormap, RuptureMap as Map
from pyrocko import model, gmtpy, util, orthodrome as pod
from pyrocko import moment_tensor as pmt
from pyrocko.io import stationxml

from ewrica.scaling import length_blaser
from ewrica.dataset.faults import load_faults_from_database
from ewrica.dataset import tsumaps


plain = True

events = model.load_events('events.txt')

gmtpy.check_have_gmt()
factor_symbol_size = 4.0
beachball_symbol = 'd'


overview_kwargs = dict(
    lat=num.mean([ev.lat for ev in events]) + 2.,
    lon=num.mean([ev.lon for ev in events]) - 9.,
    radius=1500000.,
    width=20., height=10.,
    show_grid=False,
    show_topo=True,
    color_dry=(238, 236, 230),
    topo_cpt_wet='light_sea_uniform',
    topo_cpt_dry='light_land_uniform',
    illuminate=True,
    illuminate_factor_ocean=0.15,
    show_rivers=False,
    show_plates=False)

fault_list = load_faults_from_database(database='edsf')

# OVERVIEW MAP MEDITERANEAN
if False:
    m = Map(**overview_kwargs)

    for f in fault_list:
        m.gmt.psxy(
            in_columns=(f.lon, f.lat),
            W='1.0p,{}'.format(gmtpy.color('aluminium4')),
            *m.jxyr)

    rcmt = model.load_events(filename='rcmt_full.yaml')

    m.gmt.psxy(
        in_columns=(
            [e.lon for e in rcmt],
            [e.lat for e in rcmt],
            ['{}p'.format(e.magnitude / 2.) for e in rcmt]),
        G='darkorange3',
        S='c',
        *m.jxyr)

    m.save('figures/overview_map.png', psconvert=True, resolution=150.)

# OVERVIEW MAP TSUMAPS
if False:
    m = Map(**overview_kwargs)

    tsumap = tsumaps.load()
    latlon = tsumap['latlon']

    m.gmt.psxy(
        in_columns=(latlon[:, 1], latlon[:, 0]),
        S='c1.5p',
        G='darkslategrey',
        *m.jxyr)

    m.save('figures/tsumaps_grid.png', psconvert=True, resolution=150.)

# OVERVIEW MAP TSUMAPS CUM PROB
if True:
    tsumap = tsumaps.load()
    probs = num.sort(tsumap['prob'], axis=1)[:, ::-1]

    cmap = 'YlGnBu'
    make_colormap(0., 1., cmap=cmap)

    for n in [1, 7]:
        cum_probs = num.cumsum(probs, axis=1)[:, n-1]
        p_min = num.mean(cum_probs)

        m = Map(**overview_kwargs)
        m.gmt.psxy(
            in_columns=(
                tsumap['latlon'][cum_probs < p_min, 1],
                tsumap['latlon'][cum_probs < p_min, 0],
                cum_probs),
            S='c1.75p',
            G='darkslategray',
            *m.jxyr)
        m.gmt.psxy(
            in_columns=(
                tsumap['latlon'][:, 1], tsumap['latlon'][:, 0],
                cum_probs),
            S='c1.5p',
            C='my_{}.cpt'.format(cmap),
            *m.jxyr)

        m.gmt.psscale(
            C='my_{}.cpt'.format(cmap),
            D='jBL+w{:g}c/{:g}c+h+o{:g}c/{:g}c'.format(
                overview_kwargs['width'] / 3.,
                overview_kwargs['width'] / 30.,
                overview_kwargs['width'] / 16.,
                overview_kwargs['width'] / 150. + 3 * overview_kwargs['width'] / 30.),
            # H='{}'.format(m._fontsize/2.),
            B='af+l{}'.format('cum. prob. for {} scenarios'.format(n)),
            # F='',
            *m.jxyr)

        m.save(
            op.join('figures', 'map_grid_cumprob_n_{}.png'.format(n)),
            psconvert=True,
            resolution=150.)

    clear_temp(cpts=[cmap])

# OVERVIEW MAP SAMOS + WAVE GAUGE
if False:
    ev = [ev for ev in events if ev.name.startswith('Samos')][0]

    m = Map(
        lat=ev.lat,
        lon=ev.lon - 1.,
        radius=250000.,
        width=20., height=20.,
        show_grid=False,
        show_topo=True,
        color_dry=(238, 236, 230),
        topo_cpt_wet='light_sea_uniform',
        topo_cpt_dry='light_land_uniform',
        illuminate=True,
        illuminate_factor_ocean=0.15,
        show_rivers=False,
        show_plates=False)

    mag = ev.magnitude

    devi = ev.moment_tensor.deviatoric()
    beachball_size = mag*factor_symbol_size
    mt = devi.m_up_south_east()
    mt = mt / ev.moment_tensor.scalar_moment() \
        * pmt.magnitude_to_moment(5.0)
    m6 = pmt.to6(mt)
    data = (ev.lon, ev.lat, 10) + tuple(m6) + (1, 0, 0)

    if m.gmt.is_gmt5():
        kwargs = dict(
            M=True,
            S='%s%g' % (beachball_symbol[0], (beachball_size) / gmtpy.cm))
    else:
        kwargs = dict(
            S='%s%g' % (beachball_symbol[0],
                        (beachball_size)*2 / gmtpy.cm))

    m.gmt.psmeca(
        in_rows=[data],
        G=gmtpy.color('scarletred1'),
        E='white',
        W='1p,%s' % gmtpy.color('scarletred3'),
        *m.jxyr,
        **kwargs)

    m.add_label(
        ev.lat,
        ev.lon,
        '{year} Mw {magnitude:.1f}'.format(
            year=util.tts(ev.time, format='%Y-%m-%d'),
            magnitude=mag),
        color=gmtpy.color('dimgrey'))

    m.gmt.psxy(
        in_columns=([24.9411], [37.438]),
        G=gmtpy.color('chocolate3'),
        S='c15p',
        *m.jxyr)

    m.add_label(
        37.438,
        24.9411,
        'Wave Gauge Syros',
        color=gmtpy.color('dimgrey'))

    m.save('figures/samos_map.png', psconvert=True, resolution=150.)

# SAMOS STATION MAP
if False:
    ev = [ev for ev in events if ev.name.startswith('Samos')][0]
    m = Map(
        lat=ev.effective_lat,
        lon=ev.effective_lon,
        radius=200000.,
        width=20., height=20.,
        show_grid=False,
        show_topo=True,
        color_dry=(238, 236, 230),
        topo_cpt_wet='light_sea_uniform',
        topo_cpt_dry='light_land_uniform',
        illuminate=True,
        illuminate_factor_ocean=0.15,
        show_rivers=False,
        show_plates=False)

    faults_used = load_faults_from_database(
        database='GEM',
        lat=ev.effective_lat,
        lon=ev.effective_lon,
        radius=length_blaser(magnitude=7.03, rake=-115., use_errors=False))

    faults = load_faults_from_database(
        database='GEM',
        lat=ev.effective_lat,
        lon=ev.effective_lon,
        radius=350000.)

    for f in faults:
        m.gmt.psxy(
            in_columns=(f.lon, f.lat),
            W='1.0p,{}'.format(gmtpy.color('aluminium4')),
            *m.jxyr)

    for f in faults_used:
        m.gmt.psxy(
            in_columns=(f.lon, f.lat),
            W='2.0p,{}'.format(gmtpy.color('aluminium5')),
            *m.jxyr)

    sxs = []
    for fn in [f for f in os.listdir('.')
               if f.endswith('.xml') and f.startswith('stations')]:
        sxs.append(stationxml.load_xml(filename=fn))

    sx = stationxml.primitive_merge(sxs)

    stations = sx.get_pyrocko_stations()

    # gnss_camp = model.gnss.GNSSCampaign.load(
    #     filename='../../data/events/Samos_20201030_115127/gnss/processed_ingv_static/GPS_static_Samos_2020_Betaversion_3d.yaml')

    lats_gnss, lons_gnss, labels_gnss = [], [], []
    lats_acc, lons_acc, labels_acc = [], [], []
    lats_bb, lons_bb, labels_bb = [], [], []

    dists = pod.distance_accurate50m_numpy(
        ev.effective_lat,
        ev.effective_lon,
        [s.lat for s in stations],
        [s.lon for s in stations])

    for s, dist in zip(stations, dists):
        if s.station in '094A SAMU IKAR'.split():
            lats_gnss.append(s.lat)
            lons_gnss.append(s.lon)
            labels_gnss.append(s.station)

        elif s.station in '0911 0918 0920 3511 3512 3514 3517 3523 3524 3528 3533 3538 GMLD KUSD'.split() and s.station not in labels_acc:
            lats_acc.append(s.lat)
            lons_acc.append(s.lon)
            labels_acc.append(s.station)

        elif dist >= 1.e5 and dist <= 2.e5:
            if s.station not in ['BUHA', 'AKHS', 'KIRA', 'ESEN', 'MULA']:
                continue

            resp = None
            for _, _, channel in sx.iter_network_station_channels(
                    *s.nsl(), s.channels[0].name, time=None, timespan=None):
                resp = channel.response

            if resp is None:
                continue

            if resp.instrument_sensitivity.input_units.name != 'M/S':
                continue

            lats_bb.append(s.lat)
            lons_bb.append(s.lon)
            labels_bb.append(s.station)

        else:
            continue

    # coordinates = gnss_camp.coordinates
    # m.gmt.psxy(
    #     in_columns=(coordinates[:, 1], coordinates[:, 0]),
    #     S='t20p',
    #     G=gmtpy.color('orange1'),
    #     W='1.0p,%s' % gmtpy.color('orange3'),
    #     *m.jxyr)

    m.gmt.psxy(
        in_columns=(lons_gnss, lats_gnss),
        S='t20p',
        G=gmtpy.color('skyblue1'),
        W='1.0p,%s' % gmtpy.color('skyblue3'),
        *m.jxyr)

    m.gmt.psxy(
        in_columns=(lons_acc, lats_acc),
        S='t20p',
        G=gmtpy.color('chocolate1'),
        W='1.0p,%s' % gmtpy.color('chocolate3'),
        *m.jxyr)

    m.gmt.psxy(
        in_columns=(lons_bb, lats_bb),
        S='t20p',
        G=gmtpy.color('butter1'),
        W='1.0p,%s' % gmtpy.color('butter3'),
        *m.jxyr)

    # for s in gnss_camp.stations:
    #     m.add_label(s.lat, s.lon, s.code)

    if not plain:
        for i in range(len(labels_gnss)):
            m.add_label(lats_gnss[i], lons_gnss[i], labels_gnss[i])

        for i in range(len(labels_acc)):
            m.add_label(lats_acc[i], lons_acc[i], labels_acc[i])

        for i in range(len(labels_bb)):
            m.add_label(lats_bb[i], lons_bb[i], labels_bb[i])

    mag = ev.magnitude

    devi = ev.moment_tensor.deviatoric()
    beachball_size = mag*factor_symbol_size
    mt = devi.m_up_south_east()
    mt = mt / ev.moment_tensor.scalar_moment() \
        * pmt.magnitude_to_moment(5.0)
    m6 = pmt.to6(mt)
    data = (ev.lon, ev.lat, 10) + tuple(m6) + (1, 0, 0)

    if m.gmt.is_gmt5():
        kwargs = dict(
            M=True,
            S='%s%g' % (beachball_symbol[0], (beachball_size) / gmtpy.cm))
    else:
        kwargs = dict(
            S='%s%g' % (beachball_symbol[0],
                        (beachball_size)*2 / gmtpy.cm))

    m.gmt.psmeca(
        in_rows=[data],
        G=gmtpy.color('scarletred1'),
        E='white',
        W='1p,%s' % gmtpy.color('scarletred3'),
        *m.jxyr,
        **kwargs)

    legend_kwargs = dict(
        D='g24.75/36.2+w200p+o50p/60p',
        F='+gwhite+pthicker',
        in_rows=[
            [],
            # ['S 10p T 10p {} darkslategrey 1i Static GNSS'.format(
            #     gmtpy.color('orange1'))],
            ['S 10p T 10p {} darkslategrey 1i HR GNSS'.format(
                gmtpy.color('skyblue1'))],
            ['S 10p T 10p {} darkslategrey 1i Broad Band'.format(
                gmtpy.color('butter1'))],
            ['S 10p T 10p {} darkslategrey 1i Strong Motion'.format(
                gmtpy.color('chocolate1'))]],
        J='',
        R='')

    legend_kwargs['in_rows'][0] = ['H - - stations']
    m.gmt.pslegend(**legend_kwargs)

    m.save('figures/stationmap{}.png'.format('_plain' if plain else ''),
           psconvert=True, resolution=150.)
