"""Tailored surf forecast for Jæren spots (generalised from the Point Perfect model).

Per spot: bathymetry crop (Kartverket chart depths over EMODnet) -> ray tracing (refraction + shoaling) of BarentsWatch swell
and a yr-wind-derived wind sea -> depth-limited breaking -> breaker height/depth/peel -> score.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

G = 9.81; GAMMA = 0.78; ZONE_R = 450.0

def wavenum(T, h):
    w = 2*np.pi/T; h = np.maximum(h, 0.05)
    k = w**2/G / np.sqrt(np.tanh(w**2/G*h))
    for _ in range(4):
        f = G*k*np.tanh(k*h) - w**2
        df = G*np.tanh(k*h) + G*k*h/np.cosh(k*h)**2
        k = k - f/df
    return k

def c_cg(T, h):
    k = wavenum(T, h); w = 2*np.pi/T; c = w/k
    n = 0.5*(1 + 2*k*h/np.sinh(np.minimum(2*k*h, 50)))
    return c, c*n


class Spot:
    def __init__(self, sid, lat, lon, grid, facing_hint):
        self.id, self.lat, self.lon = sid, lat, lon
        elev = grid['elev']; ny, nx = elev.shape
        lons = grid['x0'] + grid['dx']*(np.arange(nx)+.5)
        lats = grid['y0'] + grid['dy']*(np.arange(ny)+.5)
        self.mlat = 111320.0; self.mlon = 111320.0*np.cos(np.radians(lat))
        self.X = (lons-lon)*self.mlon; self.Y = (lats-lat)*self.mlat
        self.lons, self.lats, self.elev = lons, lats, elev
        depth = -np.nan_to_num(elev, nan=5.0)
        self.dsm = gaussian_filter(depth, 0.8)
        self.Yi = self.Y[::-1]
        self.itp = RegularGridInterpolator((self.Yi, self.X), self.dsm[::-1], bounds_error=False, fill_value=None)
        self.gy, self.gx = np.gradient(self.dsm, self.Y, self.X)
        # shore orientation from a broad (~1 km) smoothed depth gradient around the spot
        big = gaussian_filter(depth, 4.0)
        gy2, gx2 = np.gradient(big, self.Y, self.X)
        iy = np.argmin(abs(self.Y)); ix = np.argmin(abs(self.X))
        sx, sy = gx2[iy-3:iy+4, ix-4:ix+5].mean(), gy2[iy-3:iy+4, ix-4:ix+5].mean()
        b = float(np.degrees(np.arctan2(sx, sy)) % 360)
        diff = abs((b - facing_hint + 180) % 360 - 180)
        self.bathy_bearing = b
        self.seaward = b if diff <= 60 else float(facing_hint)
        self.offshore_wind_from = (self.seaward + 180) % 360

    def h_at(self, x, y):
        xc = np.clip(x, self.X[0], self.X[-1]); yc = np.clip(y, self.Yi[0], self.Yi[-1])
        return self.itp(np.column_stack([yc, xc]))

    def to_xy(self, lat, lon):
        return (lon-self.lon)*self.mlon, (lat-self.lat)*self.mlat

    def rel_wind(self, wdir):
        """0 = straight offshore, 180 = straight onshore."""
        return abs(((wdir - self.offshore_wind_from) + 180) % 360 - 180)

    def windsea(self, U, wdir):
        a = 180 - self.rel_wind(wdir)                       # 0 = straight onshore
        F = 200e3 if a < 70 else (40e3 if a < 100 else 2e3)
        U = max(U, 0.5); Fd = G*F/U**2
        Hs = min(0.0016*np.sqrt(Fd), 0.2433)*U**2/G
        Tp = min(0.286*Fd**(1/3), 8.13)*U/G
        return Hs, Tp

    def trace(self, H0, T, dir_from, x0=0.0, y0=0.0, n_rays=600, span=8000.0, R=8000.0, dt=1.0, nmax=3500):
        th = np.radians((dir_from + 180) % 360)
        ux, uy = np.sin(th), np.cos(th); px, py = -uy, ux
        s = np.linspace(-span, span, n_rays)
        x = x0 - R*ux + s*px; y = y0 - R*uy + s*py
        a = np.full(n_rays, np.arctan2(uy, ux))
        _, cg0 = c_cg(T, self.h_at(x, y))
        alive = np.ones(n_rays, bool)
        res = {k: np.full(n_rays, np.nan) for k in ('bx', 'by', 'hb', 'db', 'ab')}
        paths = [[(x[i], y[i])] for i in range(n_rays)]
        b0 = np.gradient(s); eps = 15.0
        for step in range(nmax):
            if not alive.any(): break
            h = self.h_at(x, y); c, cg = c_cg(T, h)
            dcx = (c_cg(T, self.h_at(x+eps, y))[0] - c_cg(T, self.h_at(x-eps, y))[0])/(2*eps)
            dcy = (c_cg(T, self.h_at(x, y+eps))[0] - c_cg(T, self.h_at(x, y-eps))[0])/(2*eps)
            b = np.hypot(np.gradient(x), np.gradient(y))
            kr = np.sqrt(np.clip(b0/np.maximum(b, 1e-3), 0.05, 4.0))
            H = H0*np.sqrt(cg0/np.maximum(cg, 1e-3))*np.minimum(kr, 1.6)
            brk = alive & ((H > GAMMA*h) | (h < 0.4))
            if brk.any():
                res['bx'][brk] = x[brk]; res['by'][brk] = y[brk]; res['hb'][brk] = H[brk]
                res['db'][brk] = h[brk]; res['ab'][brk] = a[brk]
                alive &= ~brk
            a = np.where(alive, a + (dcx*np.sin(a) - dcy*np.cos(a))*dt, a)
            x = np.where(alive, x + c*np.cos(a)*dt, x)
            y = np.where(alive, y + c*np.sin(a)*dt, y)
            if step % 25 == 0:
                for i in np.where(alive)[0]: paths[i].append((x[i], y[i]))
        for i in range(n_rays):
            if np.isfinite(res['bx'][i]): paths[i].append((res['bx'][i], res['by'][i]))
        return res, paths

    def peel(self, res, inz):
        angs = []
        for i in np.where(inz)[0]:
            jx = np.argmin(abs(self.X-res['bx'][i])); jy = np.argmin(abs(self.Y-res['by'][i]))
            g = np.array([self.gx[jy, jx], self.gy[jy, jx]])
            if np.hypot(*g) < 1e-4: continue
            ray = np.array([np.cos(res['ab'][i]), np.sin(res['ab'][i])])
            angs.append(np.degrees(np.arccos(np.clip(ray @ (-g/np.hypot(*g)), -1, 1))))
        return float(np.median(angs)) if angs else np.nan

    def slope_at(self, x, y):
        jx = np.argmin(abs(self.X-x)); jy = np.argmin(abs(self.Y-y))
        return float(max(np.hypot(self.gx[jy, jx], self.gy[jy, jx]), 0.005))


def swell_fit(dfrom, window, fade=30.0):
    """1 when the swell comes from inside the spot's window (lo, hi in degrees from, clockwise), fading
    linearly to 0 over `fade` degrees outside it."""
    lo, hi = window
    width = (hi - lo) % 360
    off = (dfrom - lo) % 360
    if off <= width: return 1.0
    out = min(off - width, 360 - off)   # degrees past hi, or before lo
    return float(max(0.0, 1 - out/fade))


def score_parts(hb, T, wind, rel, clean, irb, peel, fit=1.0):
    # breaker size: target 1.5-2.0 m (full marks), tapering to 0 at 0.3 m and to ~0.05 by 4 m
    HB_LO, HB_HI = 1.5, 2.0
    if hb < 0.3: s_size = 0.0
    elif hb < HB_LO: s_size = ((hb-0.3)/(HB_LO-0.3))**1.3
    elif hb <= HB_HI: s_size = 1.0
    else: s_size = max(0.05, 1-(hb-HB_HI)/2.0)
    s_per = float(np.clip((T-5)/7, 0, 1))
    if wind < 2: s_wind = 1.0
    elif rel < 60: s_wind = 1.0 if wind < 9 else max(0.5, 1-(wind-9)/12)
    elif rel < 120: s_wind = max(0.15, 0.9-wind/14)
    else: s_wind = max(0.05, 1-wind/9)
    s_brk = 1.0 if 0.4 <= irb <= 2.0 else 0.7
    if not np.isfinite(peel): s_peel = 0.5
    elif peel < 8: s_peel = 0.45
    elif peel <= 50: s_peel = 1.0
    else: s_peel = max(0.4, 1-(peel-50)/60)
    s_dir = 0.3 + 0.7*fit   # swell from outside the spot's window: at most 30% of the score
    total = 10*s_size*s_wind*(0.4+0.6*s_per)*(0.6+0.4*clean)*s_brk*(0.5+0.5*s_peel)*s_dir
    return dict(s_size=s_size, s_period=s_per, s_wind=s_wind, s_clean=clean, s_breaker=s_brk, s_peel=s_peel, s_dir=s_dir, score=round(float(total), 1))


def label(sc, hb):
    if hb < 0.3: return 'Flat / sheltered'
    return ['Very poor', 'Poor', 'Poor to fair', 'Fair', 'Fair to good', 'Good', 'Very good', 'Epic'][int(np.digitize(sc, [0.8, 2, 3.2, 4.4, 5.6, 7, 8.5]))]
