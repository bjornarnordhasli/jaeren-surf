// Minimal stand-in for the claude.ai Dashboard runtime, so the page runs as a plain local web page.
// Provides window.dash: data(id), onData(fn), params(), setParams(), resetParams(), colors.
// Datasets are files in ../data/ (CSV parsed with d3.csvParse, JSON as-is).
(function () {
  const FILES = {
    forecast: '../data/forecast.csv',
    wind: '../data/wind.csv',
    profile: '../data/profile.csv',
    region: '../data/region.csv',
    spots: '../data/spots.json',
    meta: '../data/meta.json',
  };
  const urlFor = (id) => FILES[id] || (id.startsWith('spatial_') ? `../data/${id}.csv` : null);
  const DEFAULTS = { spot: 'point_perfect' };
  const store = {};
  const subs = [];
  let params = { ...DEFAULTS };
  try { Object.assign(params, JSON.parse(localStorage.getItem('jaeren-surf-params') || '{}')); } catch (e) { /* no storage */ }

  let queued = false;
  function notify() {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      for (const fn of subs) { try { fn(); } catch (e) { console.error(e); } }
    });
  }

  function data(id) {
    if (store[id]) return store[id];
    const url = urlFor(id);
    if (!url) return (store[id] = { status: 'missing', data: [], message: 'Unknown dataset ' + id });
    store[id] = { status: 'loading', data: [] };
    fetch(url + '?v=' + Date.now())
      .then((r) => { if (!r.ok) throw new Error(url + ' returned ' + r.status); return url.endsWith('.json') ? r.json() : r.text().then((t) => d3.csvParse(t)); })
      .then((d) => { store[id] = { status: 'ok', data: d, columns: d.columns || Object.keys(d[0] || {}) }; notify(); })
      .catch((e) => { console.error(e); store[id] = { status: 'error', data: [], message: String(e) }; notify(); });
    return store[id];
  }

  window.dash = {
    data,
    onData(fn) { subs.push(fn); try { fn(); } catch (e) { console.error(e); } },
    params: () => ({ ...params }),
    setParams(p) { params = { ...params, ...p }; try { localStorage.setItem('jaeren-surf-params', JSON.stringify(params)); } catch (e) { /* ignore */ } notify(); },
    resetParams() { params = { ...DEFAULTS }; notify(); },
    refresh(id) { delete store[id]; data(id); },
    setLink() {}, calc() {}, loader() {},
    colors: ['#2f6db5', '#e0742c', '#7b5cc6', '#1f9a86', '#c8487a', '#8c8f22', '#4f9fd8', '#b0573a'],
  };
  window.addEventListener('resize', notify);
})();
