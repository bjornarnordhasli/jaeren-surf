"""Build grids_kartverket.json (the bathymetry run2.py uses) from Kartverket's open nautical-chart depth
data ("Dybdedata", CC BY 4.0), on exactly the same cells as the EMODnet crops in grids_emodnet.json.

Input: the Dybdedata GML files (EPSG:25832) for the municipalities around the spots (Sola 1124,
Klepp 1120, Hå 1119). Usage (from this folder):
    py build_grids_kartverket.py --download <dir>   # order them from Geonorge into <dir>, then build
    py build_grids_kartverket.py <dir>              # build from GML files already in <dir>
Only needed again when a spot is added or Kartverket updates its charts; the daily refresh doesn't run it.

Depth points: soundings (Dybdepunkt), shoals (Grunne), contour vertices (Dybdekurve, every 20 m),
coastline (Kystkontur, 0 m) and rocks (Skjær, 0 m), linearly interpolated (Delaunay). Cells inside
Kartverket's land polygons stay land (EMODnet height, at least +1 m). Cells more than MAX_GAP from any
Kartverket depth keep the EMODnet value. Depths are below chart datum (lowest tide), as charted.
"""
import glob, io, json, math, os, sys, time, urllib.request, zipfile, xml.etree.ElementTree as ET
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

MAX_GAP = 600.0   # m; beyond this from any chart depth, keep EMODnet
STEP = 20.0       # m; densify contour and coastline vertices
DYBDEDATA = '2751aacf-5472-4850-a208-3532a51c529a'  # Geonorge metadata id of Kartverket "Dybdedata"
MUNICIPALITIES = [('1124', 'Sola'), ('1120', 'Klepp'), ('1119', 'Hå')]


def download(out):
    """Order Dybdedata GML (EPSG:25832) per municipality from the Geonorge download API and unzip into out."""
    hdr = {'User-Agent': 'jaeren-surf-dashboard/1.0 (bathymetry)', 'Content-Type': 'application/json'}
    order = {'email': '', 'orderLines': [{'metadataUuid': DYBDEDATA,
             'areas': [{'code': c, 'type': 'kommune', 'name': n} for c, n in MUNICIPALITIES],
             'projections': [{'code': '25832'}], 'formats': [{'name': 'GML'}]}]}
    get = lambda u, d=None: json.loads(urllib.request.urlopen(urllib.request.Request(u, data=d, headers=hdr), timeout=120).read())
    r = get('https://nedlasting.geonorge.no/api/order', json.dumps(order).encode())
    for _ in range(40):
        if r.get('files') and all(f.get('status') == 'ReadyForDownload' for f in r['files']): break
        time.sleep(15); r = get(f"https://nedlasting.geonorge.no/api/order/{r['referenceNumber']}")
    os.makedirs(out, exist_ok=True)
    for f in r['files']:
        b = urllib.request.urlopen(urllib.request.Request(f['downloadUrl'], headers=hdr), timeout=600).read()
        zipfile.ZipFile(io.BytesIO(b)).extractall(out)
        print('downloaded', f['name'], flush=True)


def utm32(lat, lon):
    """WGS84/ETRS89 lat/lon (deg, arrays) -> UTM zone 32N easting/northing (m), Krüger series."""
    a, f, k0, lon0 = 6378137.0, 1/298.257222101, 0.9996, math.radians(9.0)
    n = f/(2-f); A = a/(1+n)*(1 + n**2/4 + n**4/64)
    al = [n/2 - 2*n**2/3 + 5*n**3/16, 13*n**2/48 - 3*n**3/5, 61*n**3/240]
    e = math.sqrt(f*(2-f))
    phi, lam = np.radians(lat), np.radians(lon) - lon0
    t = np.sinh(np.arctanh(np.sin(phi)) - e*np.arctanh(e*np.sin(phi)))
    xi, eta = np.arctan2(t, np.cos(lam)), np.arctanh(np.sin(lam)/np.sqrt(1+t**2))
    E = eta + sum(al[j]*np.cos(2*(j+1)*xi)*np.sinh(2*(j+1)*eta) for j in range(3))
    N = xi + sum(al[j]*np.sin(2*(j+1)*xi)*np.cosh(2*(j+1)*eta) for j in range(3))
    return 500000 + k0*A*E, k0*A*N


def coords(el):
    txt = ' '.join((x.text or '') for x in el.iter() if x.tag.split('}')[-1] in ('pos', 'posList'))
    v = np.array(txt.split(), dtype=float)
    return v.reshape(-1, 2)


def densify(xy, step=STEP):
    out = [xy[:1]]
    for p, q in zip(xy[:-1], xy[1:]):
        n = max(1, int(np.hypot(*(q-p))//step))
        out.append(p + (q-p)*np.arange(1, n+1)[:, None]/n)
    return np.vstack(out)


def read(paths):
    pts, land = [], []
    for path in paths:
        for _, el in ET.iterparse(path, events=('end',)):
            if el.tag.split('}')[-1] not in ('featureMember', 'member'): continue
            for f in el:
                ft = f.tag.split('}')[-1]
                kids = {c.tag.split('}')[-1]: c for c in f}
                dyb = kids.get('dybde')
                if ft in ('Dybdepunkt', 'Grunne') and dyb is not None:
                    xy = coords(kids['posisjon']); pts.append(np.column_stack([xy, np.full(len(xy), float(dyb.text))]))
                elif ft == 'Dybdekurve' and dyb is not None:
                    xy = densify(coords(kids['grense'])); pts.append(np.column_stack([xy, np.full(len(xy), float(dyb.text))]))
                elif ft == 'Kystkontur':
                    xy = densify(coords(kids['grense'])); pts.append(np.column_stack([xy, np.zeros(len(xy))]))
                elif ft == 'Skjær':
                    xy = coords(kids['posisjon']); pts.append(np.column_stack([xy, np.zeros(len(xy))]))
                elif ft == 'Landareal':
                    for ring in kids['område'].iter():
                        if ring.tag.split('}')[-1] in ('exterior', 'interior'):
                            land.append((ring.tag.endswith('exterior'), coords(ring)))
            el.clear()
    p = np.vstack(pts)
    p = p[np.unique(np.round(p[:, :2], 1), axis=0, return_index=True)[1]]
    return p, land


def inside(px, py, ring):
    """Even-odd ray casting of many points against one polygon ring."""
    x, y = ring[:, 0], ring[:, 1]
    res = np.zeros(px.shape, bool)
    x0, y0, x1, y1 = x[:-1], y[:-1], x[1:], y[1:]
    m = (px.max() >= min(x.min(), x.max())) and (px.min() <= x.max()) and (py.max() >= y.min()) and (py.min() <= y.max())
    if not m: return res
    for a, b, c, d in zip(x0, y0, x1, y1):
        cond = (b > py) != (d > py)
        if cond.any():
            xc = a + (py - b)*(c - a)/((d - b) if d != b else 1e-12)
            res ^= cond & (px < xc)
    return res


def build(grids, gml_paths):
    pts, land = read(gml_paths)
    print(f'{len(pts)} depth points, {len(land)} land rings', flush=True)
    tree = cKDTree(pts[:, :2])
    itp = LinearNDInterpolator(pts[:, :2], pts[:, 2])
    out = {}
    for sid, g in grids.items():
        elev = np.array(g['elev'], dtype=float); ny, nx = elev.shape
        lons = g['x0'] + g['dx']*(np.arange(nx)+.5); lats = g['y0'] + g['dy']*(np.arange(ny)+.5)
        LON, LAT = np.meshgrid(lons, lats)
        E, N = utm32(LAT, LON)
        gap, _ = tree.query(np.column_stack([E.ravel(), N.ravel()]))
        gap = gap.reshape(E.shape)
        depth = itp(E, N)
        island = np.zeros(E.shape, bool)
        bb = (E.min(), E.max(), N.min(), N.max())
        for ext, ring in land:
            if ring[:, 0].max() < bb[0] or ring[:, 0].min() > bb[1] or ring[:, 1].max() < bb[2] or ring[:, 1].min() > bb[3]: continue
            island ^= inside(E, N, ring)   # interiors (lakes) toggle back
        new = elev.copy()
        use = ~island & (gap <= MAX_GAP) & np.isfinite(depth)
        new[use] = -np.maximum(depth[use], 0.1)
        new[island] = np.maximum(np.nan_to_num(elev[island], nan=1.0), 1.0)
        sea_e = elev < 0
        print(f'{sid}: {use.sum()} sea cells from Kartverket, {island.sum()} land, '
              f'{(~use & ~island).sum()} kept from EMODnet; mean change on sea cells '
              f'{np.nanmean((new - elev)[use & sea_e]):+.1f} m', flush=True)
        out[sid] = dict(g, elev=np.round(new, 1).tolist())
    return out


if __name__ == '__main__':
    args = sys.argv[1:]
    if args[0] == '--download':
        args = args[1:]; download(args[0])
    src = args[0]
    here = os.path.dirname(os.path.abspath(__file__))
    dst = args[1] if len(args) > 1 else os.path.join(here, 'grids_kartverket.json')
    grids = json.load(open(os.path.join(here, 'grids_emodnet.json')))
    out = build(grids, sorted(glob.glob(os.path.join(src, '*Dybdedata*.gml'))))
    json.dump(out, open(dst, 'w'))
    print('wrote', dst)
