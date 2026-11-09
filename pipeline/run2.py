import numpy as np, pandas as pd, json, sys, os
from multiprocessing import Pool
from model2 import *

SPOTS = [  # id, name, lat, lon (Surfline), facing hint, kind, level, note
    ('sola', 'Solastranden', 58.8774, 5.5898, 255, 'Beach break', 'All levels', 'Long sandy beach by Sola airport; usually smaller than the beaches further south.'),
    ('hellesto', 'Hellestø', 58.8407, 5.5577, 290, 'Beach break', 'All levels', 'Consistent peaks both ways; some shelter from north winds.'),
    ('sele', 'Sele', 58.8121, 5.5448, 280, 'Beach break', 'All levels', 'Quieter beach between Hellestø and Bore.'),
    ('bore', 'Bore', 58.7993, 5.5452, 270, 'Beach break', 'All levels', 'Usually the biggest of the northern beaches; works on all tides.'),
    ('reve', 'Reve havn', 58.7726, 5.5166, 240, 'Beach by harbour wall', 'Intermediate', 'Beside the harbour wall near Revtangen; smaller than the open beaches.'),
    ('point_perfect', 'Point Perfect (Orre)', 58.6969, 5.5382, 230, 'Reef / point', 'Intermediate and up', 'Point break at Orre; best on south-west to south swell, north-west swell breaks poorly.'),
    ('refsnes', 'Refsnes', 58.6915, 5.5486, 240, 'Beach break', 'All levels', 'Open beach just south of Point Perfect.'),
    ('brusand', 'Brusand', 58.5300, 5.7632, 225, 'Beach break', 'All levels', 'Long southern beach; north-west wind blows more from the side here.'),
]
# Swell window per spot: swell directions (from, degrees, clockwise lo -> hi) that get full marks; the score
# fades over 30 degrees outside it to 30%. Revtangen splits the coast: the spots north of it are open to W-NW
# and sheltered from S; those south-east of it are open to SW-S, and NW swell runs along their coast at a
# grazing angle. Point Perfect is set from surfing experience (NW swell breaks poorly there, S is better); the
# others from where the mainland coastline (Kartverket) blocks open sea.
SWELL_WINDOW = {'sola': (260, 330), 'hellesto': (230, 320), 'sele': (230, 315), 'bore': (235, 325),
                'reve': (240, 340), 'point_perfect': (190, 255), 'refsnes': (200, 270), 'brusand': (180, 260)}
grids = {k: dict(v, elev=np.array(v['elev'], dtype=float)) for k, v in json.load(open('grids_kartverket.json')).items()}
bw = pd.read_csv('bw_all.csv'); bw['time'] = pd.to_datetime(bw.time)
yr = pd.read_csv('yr_all.csv'); yr['time'] = pd.to_datetime(yr.time)
pts = bw.drop_duplicates(['lat', 'lon'])[['lat', 'lon']].values
om = pd.read_csv('om_all.csv') if os.path.exists('om_all.csv') else pd.DataFrame(columns=['time', 'lat', 'lon', 'hs', 'tm', 'dir_from'])
om['time'] = pd.to_datetime(om.time)
EXT_HOURS = 144  # extended forecast: 6 days from the first BarentsWatch time
ISO = lambda t: t.strftime('%Y-%m-%dT%H:%M:%SZ')


def extend(off, ola, olo):
    """Continue the BarentsWatch series at (ola, olo) with MFWAM from the nearest Open-Meteo cell, out to EXT_HOURS.
    MFWAM gives mean period and its own height/direction bias, so scale it to BarentsWatch over the overlap."""
    off = off[['hs', 'tp', 'dir_to']].assign(wave_source='barentswatch')
    if om.empty: return off
    cells = om.drop_duplicates(['lat', 'lon'])
    d = np.hypot(cells.lat - ola, (cells.lon - olo)*np.cos(np.radians(ola)))
    c = cells.iloc[int(np.argmin(d.values))]
    m = om[(om.lat == c.lat) & (om.lon == c.lon)].drop_duplicates('time').set_index('time').sort_index()
    t = off.index.intersection(m.index)
    if len(t) >= 4:
        b, f = off.loc[t], m.loc[t]
        kh = float(np.clip(np.median(b.hs.values/np.maximum(f.hs.values, 0.05)), 0.6, 1.6))
        kt = float(np.clip(np.median(b.tp.values/np.maximum(f.tm.values, 1.0)), 1.0, 1.6))  # peak / mean period
        dd = np.radians(((b.dir_to.values + 180) - f.dir_from.values + 540) % 360 - 180)
        ddir = float(np.degrees(np.arctan2(np.mean(np.sin(dd)), np.mean(np.cos(dd)))))
    else:
        kh, kt, ddir = 1.0, 1.2, 0.0
    end = off.index.min() + pd.Timedelta(hours=EXT_HOURS)
    tail = m[(m.index > off.index.max()) & (m.index <= end)]
    ext = pd.DataFrame({'hs': (tail.hs*kh).round(2), 'tp': (tail.tm*kt).round(1),
                        'dir_to': ((tail.dir_from + ddir + 180) % 360).round(), 'wave_source': 'meteofrance'}, index=tail.index)
    return pd.concat([off, ext])


def run(spec):
    sid, name, lat, lon, hint, kind, level, note = spec
    sp = Spot(sid, lat, lon, grids[sid], hint)
    th = np.radians(sp.seaward)
    # BarentsWatch boundary point: closest to 5 km seaward, in >= 12 m of water
    tx, ty = 5000*np.sin(th), 5000*np.cos(th)
    best = None
    for la, lo in pts:
        x, y = sp.to_xy(la, lo)
        if abs(x) > 12000 or abs(y) > 12000: continue
        h = float(sp.h_at(np.array([x]), np.array([y]))[0])
        d = np.hypot(x-tx, y-ty) + (1e5 if h < 12 else 0)
        if best is None or d < best[0]: best = (d, la, lo, h)
    _, ola, olo, oh = best
    near = None
    for la, lo in pts:
        x, y = sp.to_xy(la, lo); d = np.hypot(x, y)
        h = float(sp.h_at(np.array([x]), np.array([y]))[0])
        if h > 1 and d < 2500 and (near is None or d < near[0]): near = (d, la, lo)
    off = extend(bw[(bw.lat == ola) & (bw.lon == olo)].drop_duplicates('time').set_index('time').sort_index(), ola, olo)
    nr = bw[(bw.lat == near[1]) & (bw.lon == near[2])].drop_duplicates('time').set_index('time') if near else None
    wy = yr[yr.spot_id == sid]
    PD = np.arange(0, 3001, 25); PZ = sp.h_at(PD*np.sin(th), PD*np.cos(th))
    rows, spatial = [], []
    for t, r in off.iterrows():
        w = wy.iloc[(wy.time - t).abs().argmin()]
        U, gust, wdir = float(w.wind_speed), float(w.wind_gust), float(w.wind_dir)
        Hs, Tp, dfrom = float(r.hs), float(r.tp), (float(r.dir_to)+180) % 360
        Hws, Tws = sp.windsea(U, wdir); Hws = min(Hws, 0.95*Hs)
        Hsw = float(np.sqrt(max(Hs**2-Hws**2, 0)))
        Tsw = Tp if Hsw >= Hws else max(Tp, 1.15*Tws)
        res, paths = sp.trace(Hsw, Tsw, dfrom)
        ok = np.isfinite(res['bx']); inz = ok & (np.hypot(res['bx'], res['by']) < ZONE_R)
        hb_ws = 0.0
        if Hws > 0.05:
            r2, _ = sp.trace(Hws, max(Tws, 2.5), wdir, n_rays=300)
            z2 = np.isfinite(r2['bx']) & (np.hypot(r2['bx'], r2['by']) < ZONE_R)
            if z2.any(): hb_ws = float(np.nanmedian(r2['hb'][z2]))
        if inz.any():
            hb_sw = float(np.nanmedian(res['hb'][inz])); db = float(np.nanmedian(res['db'][inz]))
            bx = float(np.nanmedian(res['bx'][inz])); by = float(np.nanmedian(res['by'][inz])); peel = sp.peel(res, inz)
        else:
            hb_sw, db, bx, by, peel = 0.0, np.nan, np.nan, np.nan, np.nan
        hb = float(np.hypot(hb_sw, hb_ws))
        slope = sp.slope_at(bx if np.isfinite(bx) else 0, by if np.isfinite(by) else 0)
        irb = slope/np.sqrt(max(hb, 0.05)/(G*Tsw**2/(2*np.pi)))
        btype = 'spilling' if irb < 0.4 else ('plunging' if irb <= 2 else 'surging')
        clean = Hsw**2/max(Hs**2, 1e-6)
        rel = sp.rel_wind(wdir)
        fit = swell_fit(dfrom, SWELL_WINDOW[sid])
        sc = score_parts(hb, Tsw, U, rel, clean, irb, peel, fit)
        wrel = 'offshore' if rel < 60 else ('cross-shore' if rel < 120 else 'onshore')
        bdist = float(PD[np.argmax(PZ >= db)]) if np.isfinite(db) and (PZ >= db).any() else None
        T = ISO(t)
        rows.append(dict(spot_id=sid, spot=name, time=T, hs_offshore=Hs, tp=Tp, wave_from=round(dfrom),
                         swell_h=round(Hsw, 2), swell_t=round(Tsw, 1), windsea_h=round(Hws, 2), windsea_t=round(Tws, 1),
                         wind_ms=U, gust_ms=gust, wind_from=round(wdir), wind_rel=wrel,
                         breaker_h=round(hb, 2), break_depth=round(db, 1) if np.isfinite(db) else None,
                         break_x=round(bx) if np.isfinite(bx) else None, break_y=round(by) if np.isfinite(by) else None,
                         break_dist=bdist, peel_angle=round(peel) if np.isfinite(peel) else None,
                         iribarren=round(irb, 2), breaker_type=btype, swell_fit=round(fit, 2),
                         bw_nearshore_hs=float(nr.loc[t].hs) if nr is not None and t in nr.index else None,
                         score=sc['score'], rating=label(sc['score'], hb), wave_source=r.wave_source))
        for i in np.where(ok)[0][::4]:
            spatial.append(dict(kind='break', time=T, ray=None, seq=None, x=round(res['bx'][i]), y=round(res['by'][i]), depth=None,
                                hb=round(float(res['hb'][i]), 2), in_zone=int(inz[i])))
        for i in np.where(ok)[0][::20]:
            P = paths[i][::2] + [paths[i][-1]]
            for j, (px, py) in enumerate(P):
                if abs(px) > 6000 or abs(py) > 6000: continue
                spatial.append(dict(kind='ray', time=T, ray=int(i), seq=j, x=round(px), y=round(py), depth=None, hb=None, in_zone=None))
    for j in range(len(sp.Y)):
        for i in range(len(sp.X)):
            if abs(sp.X[i]) <= 4000 and abs(sp.Y[j]) <= 3500 and np.isfinite(sp.elev[j, i]):
                spatial.append(dict(kind='cell', time=None, ray=None, seq=None, x=round(sp.X[i]), y=round(sp.Y[j]),
                                    depth=round(float(-sp.elev[j, i]), 1), hb=None, in_zone=None))
    prof = [dict(spot_id=sid, dist=int(d), depth=round(float(z), 1)) for d, z in zip(PD, PZ)]
    wsel = wy[(wy.time >= off.index.min()) & (wy.time <= off.index.max() + pd.Timedelta(hours=12))]
    wind = [dict(spot_id=sid, time=ISO(r.time), wind_speed=r.wind_speed, wind_gust=r.wind_gust, wind_dir=r.wind_dir,
                 wind_rel=('offshore' if sp.rel_wind(r.wind_dir) < 60 else ('cross-shore' if sp.rel_wind(r.wind_dir) < 120 else 'onshore')))
            for r in wsel.itertuples()]
    meta = dict(id=sid, name=name, lat=lat, lon=lon, kind=kind, level=level, note=note,
                swell_window=list(SWELL_WINDOW[sid]), facing=round(sp.seaward), bathy_facing=round(sp.bathy_bearing), offshore_wind_from=round(sp.offshore_wind_from),
                bw_lat=float(ola), bw_lon=float(olo), bw_depth=round(oh, 1),
                bw_near_lat=float(near[1]) if near else None, bw_near_lon=float(near[2]) if near else None)
    print(sid, meta['facing'], meta['bathy_facing'], hint, (ola, olo, round(oh)), 'scores', [r['score'] for r in rows], flush=True)
    return rows, spatial, prof, wind, meta


if __name__ == '__main__':
    with Pool(min(4, os.cpu_count() or 2)) as p: out = p.map(run, SPOTS)
    pd.DataFrame([r for o in out for r in o[0]]).to_csv('forecast.csv', index=False)
    pd.DataFrame([r for o in out for r in o[2]]).to_csv('profile.csv', index=False)
    pd.DataFrame([r for o in out for r in o[3]]).to_csv('wind.csv', index=False)
    json.dump([o[4] for o in out], open('spots.json', 'w'), ensure_ascii=False, indent=1)
    for o in out:
        pd.DataFrame(o[1]).to_csv(f"spatial_{o[4]['id']}.csv", index=False)
