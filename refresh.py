"""Refresh the local Jæren surf dashboard.

1. pipeline/fetch.py  - downloads BarentsWatch wave forecast tiles and yr (MET Norway) wind for each spot
2. pipeline/run2.py   - runs the surf model (ray tracing over bathymetry, breaking, rating) for all spots
3. copies the results into data/, which the dashboard page reads

Usage (from this folder):  python refresh.py
Then serve the folder:     python -m http.server 8000   ->  http://localhost:8000/dashboard/
"""
import os, shutil, subprocess, sys, glob, time, json, datetime as dt

ROOT = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(ROOT, 'pipeline')
DATA = os.path.join(ROOT, 'data')
OUTPUTS = ['forecast.csv', 'wind.csv', 'profile.csv', 'spots.json']
ENV = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')  # keep æ/ø UTF-8 on Windows


def step(name, script):
    t = time.time()
    print(f'--> {name} ...', flush=True)
    r = subprocess.run([sys.executable, script], cwd=PIPE, env=ENV)
    if r.returncode != 0:
        sys.exit(f'{name} failed (exit code {r.returncode}); data/ was left unchanged.')
    print(f'    done in {time.time() - t:.0f} s', flush=True)


if __name__ == '__main__':
    step('Fetching BarentsWatch waves and yr wind', 'fetch.py')
    step('Running the surf model for all spots (several minutes)', 'run2.py')
    files = OUTPUTS + [os.path.basename(p) for p in glob.glob(os.path.join(PIPE, 'spatial_*.csv'))]
    for f in files:
        shutil.copy2(os.path.join(PIPE, f), os.path.join(DATA, f))
    with open(os.path.join(DATA, 'meta.json'), 'w') as fh:
        json.dump({'refreshed_at': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}, fh)
    print(f'Updated {len(files)} files in data/. Reload the dashboard page to see them.')
