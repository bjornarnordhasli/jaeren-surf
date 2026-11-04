"""Fetch fresh inputs for the Jæren surf model.

Writes bw_all.csv (BarentsWatch Bølgevarsel wave points, 3-hourly), yr_all.csv (MET Norway
locationforecast wind per spot) and om_all.csv (Météo-France MFWAM waves via Open-Meteo, 7 days, used to
extend the forecast past BarentsWatch's ~60 h). Needs HTTPS access to www.barentswatch.no, api.met.no and
marine-api.open-meteo.com.
"""
import math, os, struct, sys, time, datetime as dt, urllib.request, json, csv, concurrent.futures as cf

# MET Norway asks for an identifying User-Agent with contact info; the GitHub workflow sets SURF_UA to the repo URL
UA = os.environ.get('SURF_UA') or 'jaeren-surf-dashboard/1.0 (personal surf forecast)'
SPOTS = {'sola': (58.8774, 5.5898), 'hellesto': (58.8407, 5.5577), 'sele': (58.8121, 5.5448), 'bore': (58.7993, 5.5452),
         'reve': (58.7726, 5.5166), 'point_perfect': (58.6969, 5.5382), 'refsnes': (58.6915, 5.5486), 'brusand': (58.5300, 5.7632)}
BOX = (5.30, 58.48, 5.85, 58.92)


def get(url, binary=False):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    with urllib.request.urlopen(req, timeout=40) as r:
        b = r.read()
        return b if binary else b.decode()


# ---- minimal Mapbox Vector Tile (protobuf) decoder: point features only ----
def _varint(b, p):
    r = s = 0
    while True:
        x = b[p]; p += 1; r |= (x & 0x7f) << s; s += 7
        if not x & 0x80: return r, p

def _fields(b, p, end):
    while p < end:
        k, p = _varint(b, p); f, w = k >> 3, k & 7
        if w == 0: v, p = _varint(b, p); yield f, w, v
        elif w == 1: yield f, w, b[p:p+8]; p += 8
        elif w == 2: l, p = _varint(b, p); yield f, w, b[p:p+l]; p += l
        elif w == 5: yield f, w, b[p:p+4]; p += 4
        else: raise ValueError('wire type %d' % w)

def decode_points(buf, z, tx, ty):
    n = 2**z; out = []
    zz = lambda v: -((v+1) >> 1) if v & 1 else v >> 1
    for f, w, lay in _fields(buf, 0, len(buf)):
        if f != 3: continue
        keys, vals, feats, extent = [], [], [], 4096
        for f2, w2, v in _fields(lay, 0, len(lay)):
            if f2 == 3: keys.append(v.decode())
            elif f2 == 4:
                val = None
                for f3, w3, x in _fields(v, 0, len(v)):
                    if f3 == 1: val = x.decode()
                    elif f3 == 2: val = struct.unpack('<f', x)[0]
                    elif f3 == 3: val = struct.unpack('<d', x)[0]
                    elif f3 in (4, 5): val = x
                    elif f3 == 6: val = zz(x)
                    elif f3 == 7: val = bool(x)
                vals.append(val)
            elif f2 == 2:
                tags, geom = [], []
                for f3, w3, x in _fields(v, 0, len(v)):
                    if f3 in (2, 4):
                        arr, p = [], 0
                        while p < len(x): a, p = _varint(x, p); arr.append(a)
                        if f3 == 2: tags = arr
                        else: geom = arr
                feats.append((tags, geom))
            elif f2 == 5: extent = v
        for tags, geom in feats:
            if len(geom) < 3: continue
            o = {keys[tags[i]]: vals[tags[i+1]] for i in range(0, len(tags), 2)}
            px = (tx + zz(geom[1])/extent)/n; py = (ty + zz(geom[2])/extent)/n
            o['lon'] = round(px*360-180, 4); o['lat'] = round(math.degrees(math.atan(math.sinh(math.pi*(1-2*py)))), 4)
            out.append(o)
    return out

def tile_xy(lon, lat, z):
    n = 2**z; lr = math.radians(lat)
    return int((lon+180)/360*n), int((1-math.log(math.tan(lr)+1/math.cos(lr))/math.pi)/2*n)


def fetch_bw(hours=72):
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
    start = now - dt.timedelta(hours=now.hour % 3)
    times = [start + dt.timedelta(hours=3*i) for i in range(hours//3 + 1)]
    tiles = []
    for z in (10, 11):
        x0, y0 = tile_xy(BOX[0], BOX[3], z); x1, y1 = tile_xy(BOX[2], BOX[1], z)
        tiles += [(z, x, y) for x in range(x0, x1+1) for y in range(y0, y1+1)]
    jobs = [(t, z, x, y) for t in times for (z, x, y) in tiles]
    def one(j):
        t, z, x, y = j
        url = f"https://www.barentswatch.no/bwapi/v1/tiles/waveforecastpoints/{t.strftime('%Y-%m-%dT%H:%M:%S.000Z')}/{z}/{x}/{y}.pbf"
        try: return t, z, decode_points(get(url, True), z, x, y)
        except Exception as e: return t, z, []
    rows, seen = [], set()
    with cf.ThreadPoolExecutor(8) as ex:
        for t, z, pts in ex.map(one, jobs):
            for o in pts:
                if not (BOX[1] <= o['lat'] <= BOX[3] and BOX[0] <= o['lon'] <= BOX[2]): continue
                if o.get('total_significant_wave_height') is None: continue
                k = (t, o['lat'], o['lon'])
                if k in seen: continue
                seen.add(k)
                rows.append([t.strftime('%Y-%m-%dT%H:%M:%SZ'), o['lat'], o['lon'], o['total_significant_wave_height'],
                             o['total_peak_period'], o['total_mean_wave_direction'], z])
    rows.sort()
    with open('bw_all.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time', 'lat', 'lon', 'hs', 'tp', 'dir_to', 'zoom']); w.writerows(rows)
    return len(rows), len({r[0] for r in rows})


def fetch_yr():
    rows = []
    for sid, (la, lo) in SPOTS.items():
        # api.met.no sometimes serves a truncated (cached) body or a 5xx for one exact URL; retry, then
        # ask for the point rounded to 3 decimals (~50 m away), which is a separate cache entry
        tries = [(la, lo)]*2 + [(round(la, 3), round(lo, 3))]*3 + [(round(la, 2), round(lo, 2))]
        for k, (qa, qo) in enumerate(tries):
            try:
                j = json.loads(get(f'https://api.met.no/weatherapi/locationforecast/2.0/complete?lat={qa}&lon={qo}'))
                break
            except (ValueError, OSError) as e:
                print(f'yr {sid} ({qa}, {qo}) failed: {str(e)[:80]}', flush=True)
                if k == len(tries) - 1: raise
                time.sleep(10)
        up = j['properties']['meta']['updated_at']
        for s in j['properties']['timeseries'][:90]:
            d = s['data']['instant']['details']
            pr = (s['data'].get('next_1_hours') or s['data'].get('next_6_hours') or {}).get('details', {}).get('precipitation_amount', '')
            rows.append([sid, s['time'], d['wind_speed'], d.get('wind_speed_of_gust', ''), d['wind_from_direction'], d.get('air_temperature', ''), pr, up])
    with open('yr_all.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['spot_id', 'time', 'wind_speed', 'wind_gust', 'wind_dir', 'air_temp', 'precip', 'updated']); w.writerows(rows)
    return len(rows)


def fetch_om(days=7, step=0.05):
    """MFWAM total sea on a lat/lon lattice over BOX; Open-Meteo snaps each to its nearest model cell."""
    lats = [round(BOX[1] + step*i, 3) for i in range(int((BOX[3]-BOX[1])/step) + 2)]
    lons = [round(BOX[0] + step*i, 3) for i in range(int((BOX[2]-BOX[0])/step) + 2)]
    pts = [(la, lo) for la in lats for lo in lons]
    rows, seen = [], set()
    for k in range(0, len(pts), 50):
        chunk = pts[k:k+50]
        url = ('https://marine-api.open-meteo.com/v1/marine?models=meteofrance_wave&timezone=GMT'
               f'&forecast_days={days}&hourly=wave_height,wave_period,wave_direction'
               f"&latitude={','.join(str(p[0]) for p in chunk)}&longitude={','.join(str(p[1]) for p in chunk)}")
        j = json.loads(get(url))
        for o in (j if isinstance(j, list) else [j]):
            cell = (round(o['latitude'], 4), round(o['longitude'], 4))
            if cell in seen: continue
            seen.add(cell)
            h = o['hourly']
            for t, hs, tm, d in zip(h['time'], h['wave_height'], h['wave_period'], h['wave_direction']):
                if hs is None or tm is None or d is None or int(t[11:13]) % 3: continue
                rows.append([t + ':00Z', cell[0], cell[1], hs, tm, d])
    with open('om_all.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time', 'lat', 'lon', 'hs', 'tm', 'dir_from']); w.writerows(rows)
    return len(rows), len({(r[1], r[2]) for r in rows})


if __name__ == '__main__':
    print('yr rows', fetch_yr())
    n, nt = fetch_bw()
    print('barentswatch rows', n, 'times', nt)
    if nt < 4: sys.exit('BarentsWatch returned too few forecast times')
    try:
        n, nc = fetch_om()
        print('open-meteo MFWAM rows', n, 'cells', nc)
    except Exception as e:  # the extended forecast is optional; the core forecast still runs without it
        print('open-meteo MFWAM fetch failed, no extended forecast:', e)
        open('om_all.csv', 'w').write('time,lat,lon,hs,tm,dir_from\n')
