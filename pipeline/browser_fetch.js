// Fallback fetch when the cloud sandbox can't reach api.met.no / barentswatch.no.
// Run in the built-in browser on a tab at https://api.met.no/ (the BarentsWatch tiles allow cross-origin reads).
// Returns "@@yr\n<csv>\n@@bw\n<csv>\n#...." (padded so the tool saves it to a file).
(async () => {
  const SP = { sola: [58.8774, 5.5898], hellesto: [58.8407, 5.5577], sele: [58.8121, 5.5448], bore: [58.7993, 5.5452], reve: [58.7726, 5.5166], point_perfect: [58.6969, 5.5382], refsnes: [58.6915, 5.5486], brusand: [58.53, 5.7632] };
  const yr = ['spot_id,time,wind_speed,wind_gust,wind_dir,air_temp,precip,updated'];
  for (const [id, [la, lo]] of Object.entries(SP)) {
    const j = await fetch(`/weatherapi/locationforecast/2.0/complete?lat=${la}&lon=${lo}`).then((r) => r.json());
    for (const s of j.properties.timeseries.slice(0, 90)) { const d = s.data.instant.details; yr.push([id, s.time, d.wind_speed, d.wind_speed_of_gust ?? '', d.wind_from_direction, d.air_temperature ?? '', s.data.next_1_hours?.details?.precipitation_amount ?? s.data.next_6_hours?.details?.precipitation_amount ?? '', j.properties.meta.updated_at].join(',')); }
  }
  function rd(buf) { let p = 0; const u = new Uint8Array(buf), dv = new DataView(buf);
    const vi = () => { let r = 0, s = 0, b; do { b = u[p++]; r += (b & 0x7f) * 2 ** s; s += 7; } while (b & 0x80); return r; };
    const msg = (end, cb) => { while (p < end) { const k = vi(); cb(Math.floor(k / 8), k & 7); } };
    const skip = (w) => { if (w == 0) vi(); else if (w == 1) p += 8; else if (w == 2) { const n = vi(); p += n; } else if (w == 5) p += 4; };
    const td = new TextDecoder(), L = [];
    msg(u.length, (f, w) => { if (f == 3 && w == 2) { const ln = vi(), end = p + ln; const l = { keys: [], vals: [], feats: [], extent: 4096 };
      msg(end, (f2, w2) => { if (f2 == 3) { const n = vi(); l.keys.push(td.decode(u.slice(p, p + n))); p += n; }
        else if (f2 == 4) { const l2 = vi(), e2 = p + l2; let v = null; msg(e2, (f3, w3) => { if (f3 == 1) { const n = vi(); v = td.decode(u.slice(p, p + n)); p += n; } else if (f3 == 2) { v = dv.getFloat32(p, true); p += 4; } else if (f3 == 3) { v = dv.getFloat64(p, true); p += 8; } else if (f3 == 4 || f3 == 5) v = vi(); else if (f3 == 6) { const z = vi(); v = z % 2 ? -(z + 1) / 2 : z / 2; } else if (f3 == 7) v = !!vi(); else skip(w3); }); l.vals.push(v); }
        else if (f2 == 2) { const l2 = vi(), e2 = p + l2; const F = { tags: [], geom: [] }; msg(e2, (f3, w3) => { if (f3 == 2 || f3 == 4) { const l3 = vi(), e3 = p + l3; const a = []; while (p < e3) a.push(vi()); if (f3 == 2) F.tags = a; else F.geom = a; } else skip(w3); }); l.feats.push(F); }
        else if (f2 == 5) l.extent = vi(); else skip(w2); }); L.push(l); } else skip(w); });
    return L; }
  const tx = (lon, z) => Math.floor((lon + 180) / 360 * 2 ** z), ty = (lat, z) => { const r = lat * Math.PI / 180; return Math.floor((1 - Math.log(Math.tan(r) + 1 / Math.cos(r)) / Math.PI) / 2 * 2 ** z); };
  const now = new Date(); now.setUTCMinutes(0, 0, 0); now.setUTCHours(now.getUTCHours() - (now.getUTCHours() % 3));
  const bw = ['time,lat,lon,hs,tp,dir_to,zoom'], seen = new Set();
  for (let i = 0; i <= 24; i++) {
    const t = new Date(now.getTime() + i * 3 * 3600e3), iso = t.toISOString(), key = iso.slice(0, 19) + 'Z';
    for (const z of [10, 11]) for (let x = tx(5.30, z); x <= tx(5.85, z); x++) for (let y = ty(58.92, z); y <= ty(58.48, z); y++) {
      const r = await fetch(`https://www.barentswatch.no/bwapi/v1/tiles/waveforecastpoints/${iso}/${z}/${x}/${y}.pbf`, { credentials: 'omit' });
      if (!r.ok) continue; const n = 2 ** z, zz = (v) => (v % 2 ? -(v + 1) / 2 : v / 2);
      for (const l of rd(await r.arrayBuffer())) for (const F of l.feats) { if (F.geom.length < 3) continue; const o = {}; for (let k = 0; k < F.tags.length; k += 2) o[l.keys[F.tags[k]]] = l.vals[F.tags[k + 1]];
        const lon = +((x + zz(F.geom[1]) / l.extent) / n * 360 - 180).toFixed(4), lat = +(Math.atan(Math.sinh(Math.PI * (1 - 2 * (y + zz(F.geom[2]) / l.extent) / n))) * 180 / Math.PI).toFixed(4);
        if (lat < 58.48 || lat > 58.92 || lon < 5.30 || lon > 5.85 || o.total_significant_wave_height == null) continue;
        const s = key + lat + ',' + lon; if (seen.has(s)) continue; seen.add(s);
        bw.push([key, lat, lon, o.total_significant_wave_height, o.total_peak_period, o.total_mean_wave_direction, z].join(',')); }
    }
  }
  return '@@yr\n' + yr.join('\n') + '\n@@bw\n' + bw.join('\n') + '\n#' + '.'.repeat(52000);
})()
