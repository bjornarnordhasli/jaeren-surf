# Jæren surf forecast dashboard

Tailored surf forecast for eight Jæren spots (Solastranden, Hellestø, Sele, Bore, Reve havn,
Point Perfect (Orre), Refsnes, Brusand). Exported from the claude.ai dashboard of the same name.

## Run it

```
pip install -r requirements.txt
python refresh.py                 # fetch new forecasts + run the model (takes several minutes)
python -m http.server 8000        # then open http://localhost:8000/dashboard/
```

The page must be served over http (opening `index.html` directly from disk blocks the data loads).

## On GitHub (phone access)

`.github/workflows/refresh.yml` runs `refresh.py` on GitHub's servers every morning (04:30 UTC), on
every push to `main`, and on demand, then publishes `dashboard/` and `data/` to GitHub Pages at
`https://<user>.github.io/<repo>/`. To refresh from a phone: GitHub app → the repo → Actions →
"Refresh forecast" → Run workflow; the page shows "Forecast updated …" when the new data is live
(about 10 minutes). If a run fails it waits 15 minutes and tries once more. One-time setup: repo Settings → Pages → Source: "GitHub Actions".
`data/` already contains the forecast from 28 Oct 2026, so the page works before the first refresh.

## Layout

| Path | What it is |
|---|---|
| `dashboard/index.html` | The dashboard page (HTML, CSS and D3 drawing code in one file) |
| `dashboard/dash-shim.js` | Small stand-in for the claude.ai dashboard runtime: `dash.data(id)` loads `../data/<id>.csv`/`.json`, plus `onData`, `params`/`setParams`, `colors` |
| `dashboard/vendor/d3.min.js` | D3 v7.9.0 (ISC licence) |
| `data/` | Model outputs the page reads: `forecast.csv`, `wind.csv`, `profile.csv`, `region.csv`, `spots.json`, `spatial_<spot>.csv` |
| `pipeline/fetch.py` | Downloads BarentsWatch Bølgevarsel wave tiles (public vector tiles), MET Norway locationforecast (yr) wind and Météo-France MFWAM waves from Open-Meteo (for the 6-day view) |
| `pipeline/model2.py` | The surf model: wave rays over bathymetry, shoaling, breaking (H > 0.78 h), rating |
| `pipeline/run2.py` | Runs the model for every spot and writes the output files; spot list is at the top |
| `pipeline/grids_kartverket.json` | Bathymetry the model uses: Kartverket chart depths gridded per spot (built by `build_grids_kartverket.py`) |
| `pipeline/grids_emodnet.json` | EMODnet bathymetry crops around each spot: the cell layout, land heights and fallback depths |
| `pipeline/build_grids_kartverket.py` | Downloads Kartverket "Dybdedata" (Sola, Klepp, Hå) from Geonorge and builds `grids_kartverket.json`; run it again only when adding a spot |
| `pipeline/browser_fetch.js`, `split_fallback.py` | Fallback fetch through a browser tab, for machines that can't reach the APIs directly |
| `refresh.py` | Runs fetch + model and copies results into `data/` |

## The model, briefly

- **Waves:** BarentsWatch Jæren wave model, about 5 km offshore of each spot (Hs, Tp, direction "towards").
  It runs about 60 h ahead. The page's "6 days" option continues each spot with Météo-France MFWAM
  (Open-Meteo marine API, nearest model cell), scaled to BarentsWatch over the hours both cover; those
  rows have `wave_source=meteofrance` in `forecast.csv`. If Open-Meteo can't be reached, the refresh
  still runs and the 6-day option stays off.
- **Wind:** MET Norway locationforecast 2.0 at each spot (the forecast behind yr). If MET fails for a spot
  after retries, that spot uses the nearest spot's wind; the refresh only stops if every spot fails.
- **Swell vs wind sea:** wind sea estimated from yr wind (JONSWAP fetch-limited growth); the rest is swell.
- **Bathymetry:** Kartverket nautical-chart depths ("Dybdedata", © Kartverket, CC BY 4.0: soundings,
  shoals, contours, coastline) interpolated onto the EMODnet cells (~60 × 115 m); EMODnet fills cells
  more than 600 m from any charted depth. Depths are below chart datum. Compared with EMODnet alone the
  chart data is 1–2.5 m deeper within ~500 m of the beaches, so waves break closer in and slightly bigger.
- **Transformation:** linear-wave ray tracing over that bathymetry, refraction,
  shoaling and ray-spacing; a ray breaks where H > 0.78 × depth. Breaks within 450 m of the spot count.
- **Rating 0–10:** breaker size (full marks for 1.5–2 m), period, wind vs the spot's shoreline,
  swell share, breaker type (Iribarren) and peel angle.
- Not included: tide; reef-scale bathymetry (Kartverket's modern multibeam surveys off Sola and
  Hellestø are classified, and elsewhere the inshore soundings are from 1979–85).
