'use strict';
// ───── Análisis de la base: BPM, dónde empieza el compás, bombo y caja en cada compás ─────
// Todo en el móvil. Entra: audio mono (Float32Array) y su frecuencia de muestreo.
// Sale: { bpm, off (s hasta el primer compás), intro (compases antes de la batería),
//         bars: [{k:[16], s:[16]}] golpes por paso (0..1), pat: patrón típico, changes: [compás], conf }

(function (G) {
  function fft(re, im) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const a = -2 * Math.PI / len, wr = Math.cos(a), wi = Math.sin(a), h = len >> 1;
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < h; k++) {
          const p = i + k, q = p + h, xr = re[q] * cr - im[q] * ci, xi = re[q] * ci + im[q] * cr;
          re[q] = re[p] - xr; im[q] = im[p] - xi; re[p] += xr; im[p] += xi;
          const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
  }
  const pause = () => new Promise(r => setTimeout(r, 0));
  const pct = (arr, p) => { const a = Array.from(arr).sort((x, y) => x - y); return a.length ? a[Math.min(a.length - 1, Math.floor(p * a.length))] : 0; };

  async function analyzeBeat(input, sr0, onProgress) {
    // a ~22 kHz para ir rápido
    const dec = sr0 > 32000 ? 2 : 1, sr = sr0 / dec, len = Math.floor(input.length / dec);
    const x = new Float32Array(len);
    for (let i = 0; i < len; i++) { let s = 0; for (let j = 0; j < dec; j++) s += input[i * dec + j]; x[i] = s / dec; }
    const N = 1024, H = 256, fr = sr / H, nF = Math.max(0, Math.floor((len - N) / H));
    const win = new Float32Array(N); for (let i = 0; i < N; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / N);
    const bin = f => Math.round(f * N / sr);
    const B = { k: [bin(40), bin(130)], s: [bin(1500), bin(5000)], b: [bin(150), bin(400)], h: [bin(7000), Math.min(N / 2 - 1, bin(11000))] };
    const K = new Float32Array(nF), S = new Float32Array(nF), Hh = new Float32Array(nF), E = new Float32Array(nF), LowE = new Float32Array(nF);
    const re = new Float32Array(N), im = new Float32Array(N), L = new Float32Array(N / 2), Lp = new Float32Array(N / 2);
    // para distinguir partes: timbre (12 bandas) y armonía (12 notas) de cada frame
    const NFEAT = 24, FT = new Float32Array(nF * NFEAT);
    const bandEdges = []; for (let i = 0; i <= 12; i++) bandEdges.push(Math.max(1, bin(60 * Math.pow(10000 / 60, i / 12))));
    const pcOf = new Int8Array(N / 2).fill(-1);
    for (let k = bin(100); k <= bin(2000); k++) pcOf[k] = ((Math.round(12 * Math.log2((k * sr / N) / 440)) % 12) + 12 + 9) % 12;
    for (let f = 0; f < nF; f++) {
      const o = f * H; let e = 0;
      for (let i = 0; i < N; i++) { const v = x[o + i]; re[i] = v * win[i]; im[i] = 0; e += v * v; }
      fft(re, im);
      for (let k = 0; k < N / 2; k++) L[k] = Math.log1p(100 * Math.hypot(re[k], im[k]));
      const fl = (a, b) => { let s = 0; for (let k = a; k <= b; k++) { const d = L[k] - Lp[k]; if (d > 0) s += d; } return s / (b - a + 1); };
      if (f > 0) { K[f] = fl(B.k[0], B.k[1]); S[f] = fl(B.s[0], B.s[1]) + 0.4 * fl(B.b[0], B.b[1]); Hh[f] = fl(B.h[0], B.h[1]); }
      let le = 0; for (let k = B.k[0]; k <= B.k[1]; k++) le += L[k]; LowE[f] = le;
      E[f] = Math.sqrt(e / N);
      { const o2 = f * NFEAT;
        for (let b = 0; b < 12; b++) { let s2 = 0; const a = bandEdges[b], z = Math.max(a, bandEdges[b + 1] - 1); for (let k = a; k <= z; k++) s2 += L[k]; FT[o2 + b] = s2 / (z - a + 1); }
        for (let k = bin(100); k <= bin(2000); k++) FT[o2 + 12 + pcOf[k]] += L[k]; }
      Lp.set(L);
      if (f % 3000 === 0) { onProgress && onProgress(0.7 * f / nF); await pause(); }
    }
    // curva de golpes sin la media local
    const O = new Float32Array(nF), w = Math.round(fr * 0.4);
    { let acc = 0; const raw = new Float32Array(nF); for (let f = 0; f < nF; f++) raw[f] = K[f] + S[f] + 0.5 * Hh[f];
      for (let f = 0; f < nF; f++) { acc += raw[f]; if (f > 2 * w) acc -= raw[f - 2 * w - 1]; const m = acc / Math.min(f + 1, 2 * w + 1); O[f] = Math.max(0, raw[f] - m); } }
    // tempo: autocorrelación con preferencia por ~92 BPM (rap)
    const maxLag = Math.ceil(60 * fr / 55) * 4 + 4, ac = new Float32Array(maxLag);
    for (let l = 1; l < maxLag; l++) { let s = 0; for (let f = l; f < nF; f++) s += O[f] * O[f - l]; ac[l] = s / (nF - l); }
    const acAt = l => { const i = Math.floor(l), t = l - i; return i + 1 < maxLag ? ac[i] * (1 - t) + ac[i + 1] * t : 0; };
    let best = 90, bestS = -1;
    for (let b = 60; b <= 180; b += 0.25) {
      const lag = 60 * fr / b; let s = 0;
      for (let k = 1; k <= 4; k++) s += acAt(lag * k) / k;
      const prior = Math.exp(-0.5 * Math.pow(Math.log2(b / 92) / 0.45, 2));
      if (s * prior > bestS) { bestS = s * prior; best = b; }
    }
    onProgress && onProgress(0.8); await pause();
    // fase y BPM fino: caja en 2 y 4, bombo en el 1
    const NB = 128; let top = { score: -1e9 };
    const norm = a => { const m = pct(a, 0.98) || 1; return Float32Array.from(a, v => Math.min(1.5, v / m)); };
    const Kn = norm(K), Sn = norm(S);
    for (let b = best - 0.75; b <= best + 0.75; b += 0.05) {
      const barF = 16 * 60 * fr / b / 4, fK = new Float32Array(NB), fS = new Float32Array(NB);
      for (let f = 0; f < nF; f++) { const p = Math.floor(((f % barF) / barF) * NB) % NB; fK[p] += Kn[f]; fS[p] += Sn[f]; }
      const g = (a, p) => a[(p + NB) % NB] + 0.5 * (a[(p + 1) % NB] + a[(p + NB - 1) % NB]);
      for (let p = 0; p < NB; p++) {
        const sc = g(fS, p + 32) + g(fS, p + 96) + 0.6 * g(fK, p) - 0.5 * (g(fS, p) + g(fS, p + 64));
        const peak = (sc / (nF / barF));   // por compás, para comparar tempos
        if (peak > top.score) top = { score: peak, bpm: b, p };
      }
    }
    const bpm = Math.round(top.bpm * 10) / 10, sf = 60 * fr / top.bpm / 4, barF = 16 * sf;
    let phi = top.p / NB * barF;
    // ajuste fino de la fase al golpe real (±1/2 paso)
    { let bs = -1, bo = 0; for (let d = -Math.round(sf / 2); d <= Math.round(sf / 2); d++) { let s = 0; for (let t = phi + d; t < nF; t += barF) { const i = Math.round(t); if (i >= 0 && i < nF) s += Kn[i] + 0.5 * Sn[Math.min(nF - 1, Math.round(i + 4 * sf))]; } if (s > bs) { bs = s; bo = d; } } phi = ((phi + bo) % barF + barF) % barF; }
    onProgress && onProgress(0.9); await pause();
    // golpes compás a compás
    const nBars = Math.floor((nF - phi) / barF);
    const mx = (a, t) => { let m = 0; for (let i = Math.round(t) - 2; i <= Math.round(t) + 2; i++) if (i >= 0 && i < nF && a[i] > m) m = a[i]; return m; };
    const kAll = [], sAll = [], raw = [];
    for (let n = 0; n < nBars; n++) {
      const kk = [], ss = [];
      for (let s = 0; s < 16; s++) { const t = phi + (n * 16 + s) * sf; kk.push(mx(Kn, t)); ss.push(mx(Sn, t)); }
      const dk = []; for (let s = 0; s < 16; s++) { const t = Math.round(phi + (n * 16 + s) * sf); let pk = -1e9, pi = t; for (let i = t - 2; i <= t + 3; i++) if (i >= 0 && i < nF && LowE[i] > pk) { pk = LowE[i]; pi = i; } const j = Math.min(nF - 1, pi + Math.round(fr * 0.12)); dk.push(pk - LowE[j]); }
      raw.push({ kk, ss, dk }); kk.forEach(v => kAll.push(v)); ss.forEach(v => sAll.push(v));
    }
    const kRef = pct(kAll, 0.95) || 1, sRef = pct(raw.map(r => Math.max(r.ss[4], r.ss[12])), 0.75) || 1;
    // bombo: golpe grave fuerte; caja: golpe medio/agudo fuerte que no es solo bombo
    // el bajo también suena en graves: el bombo se distingue porque cae rápido
    const strong = []; raw.forEach(r => r.kk.forEach((v, i) => { if (v / kRef >= 0.7) strong.push(r.dk[i]); }));
    const dkRef = Math.max(1, pct(strong, 0.7));
    const bars = raw.map(r => {
      const k = [], s = [];
      for (let i = 0; i < 16; i++) {
        const fall = 1;
        const kq = Math.min(1, r.kk[i] / kRef) * fall, sq = Math.min(1, r.ss[i] / sRef), both = kq >= 0.95 && sq >= 0.95;
        k.push(kq >= 0.4 && (kq > sq * 0.95 || both) ? Math.round(kq * 100) / 100 : 0);
        s.push(sq >= 0.55 && (sq > kq * 1.05 || both) ? Math.round(sq * 100) / 100 : 0);
      }
      return { k, s };
    });
    // cuándo entra la batería
    const be = bars.map(b => b.k.reduce((a, v) => a + v, 0) + b.s.reduce((a, v) => a + v, 0));
    const ref = pct(be, 0.75) || 1; let intro = be.findIndex(v => v > 0.35 * ref); if (intro < 0) intro = 0;
    const drum = bars.slice(intro);
    const pat = { k: [], s: [] };
    for (let s = 0; s < 16; s++) { pat.k.push(drum.filter(b => b.k[s]).length / Math.max(1, drum.length)); pat.s.push(drum.filter(b => b.s[s]).length / Math.max(1, drum.length)); }
    // cambios de parte: energía de 4 en 4 compases desde que entra la batería
    const barE = []; for (let n = 0; n < nBars; n++) { let s = 0, c = 0; for (let f = Math.round(phi + n * barF); f < Math.round(phi + (n + 1) * barF) && f < nF; f++) { s += E[f]; c++; } barE.push(c ? s / c : 0); }
    const changes = [], eMax = Math.max(...barE, 1e-9);
    for (let n = intro + 4; n + 4 <= nBars; n += 4) {
      const a = barE.slice(n - 4, n).reduce((x, y) => x + y, 0) / 4, b = barE.slice(n, n + 4).reduce((x, y) => x + y, 0) / 4;
      const pa = be.slice(n - 4, n).reduce((x, y) => x + y, 0) / 4, pb = be.slice(n, n + 4).reduce((x, y) => x + y, 0) / 4;
      if (Math.abs(b - a) / eMax > 0.18 || Math.abs(pb - pa) / (ref || 1) > 0.45) changes.push({ bar: n, up: b > a });
    }
    const parts = segmentParts(FT, NFEAT, nBars, phi, barF, nF, barE, be, intro);
    // confianza: qué claro está la caja en 2 y 4
    const conf = Math.max(0, Math.min(1, (pat.s[4] + pat.s[12]) / 2));
    onProgress && onProgress(1);
    let off = phi / fr; const barS = barF / fr;
    if (off > barS - 0.08) { off = Math.max(0, off - barS); bars.unshift({ k: new Array(16).fill(0), s: new Array(16).fill(0) }); intro += 1; changes.forEach(c => c.bar++); }
    parts.forEach(p => { if (off !== phi / fr) p.start += 1; }); if (off !== phi / fr && parts.length) { parts[0].start = 0; parts[0].len += 1; }
    return { bpm, off, parts, intro, bars, pat, changes, conf, dur: len / sr, v: 2 };
  }

  // ── Partes del beat: bloques de 4 compases, cambios de timbre/armonía/energía y repeticiones ──
  function segmentParts(FT, D, nBars, phi, barF, nF, barE, be, intro) {
    if (nBars < 8) return [{ start: 0, len: nBars, label: 'Beat', g: 0 }];
    // vector por compás
    const V = [];
    for (let n = 0; n < nBars; n++) {
      const v = new Float32Array(D); let c = 0;
      for (let f = Math.round(phi + n * barF); f < Math.round(phi + (n + 1) * barF) && f < nF; f++, c++) for (let d = 0; d < D; d++) v[d] += FT[f * D + d];
      if (c) for (let d = 0; d < D; d++) v[d] /= c;
      // armonía normalizada (que no dependa del volumen)
      let cs = 0; for (let d = 12; d < 24; d++) cs += v[d]; if (cs > 0) for (let d = 12; d < 24; d++) v[d] = v[d] / cs * 12;
      V.push(v);
    }
    // z-score por dimensión
    for (let d = 0; d < D; d++) { let m = 0, q = 0; V.forEach(v => { m += v[d]; }); m /= nBars; V.forEach(v => { q += (v[d] - m) ** 2; }); const sd = Math.sqrt(q / nBars) || 1; V.forEach(v => { v[d] = (v[d] - m) / sd; }); }
    const eMax = Math.max(...barE, 1e-9), dMax = Math.max(...be, 1e-9);
    const cos = (a, b) => { let x = 0, y = 0, z = 0; for (let i = 0; i < a.length; i++) { x += a[i] * b[i]; y += a[i] * a[i]; z += b[i] * b[i]; } return x / (Math.sqrt(y * z) || 1); };
    const blockVec = (a, b) => { const v = new Float32Array(D + 2); for (let n = a; n < b; n++) for (let d = 0; d < D; d++) v[d] += V[n][d] / (b - a); let e = 0, k = 0; for (let n = a; n < b; n++) { e += barE[n] / eMax; k += be[n] / dMax; } v[D] = 3 * e / (b - a); v[D + 1] = 3 * k / (b - a); return v; };
    // bloques de 4 compases alineados con la entrada de la batería
    const first = intro % 4, starts = [];
    if (first) starts.push(0);
    for (let n = first; n < nBars; n += 4) starts.push(n);
    const blocks = starts.map((a, i) => ({ a, b: Math.min(nBars, starts[i + 1] ?? nBars) })).filter(x => x.b > x.a);
    blocks.forEach(x => { x.v = blockVec(x.a, x.b); x.e = x.v[D] / 3; x.k = x.v[D + 1] / 3; });
    // novedad: cuánto cambia el sonido de un bloque al siguiente (timbre + armonía + batería)
    const nov = blocks.map((x, i) => (i ? 1 - cos(blocks[i - 1].v.slice(0, D), x.v.slice(0, D)) + 1.2 * Math.abs(x.k - blocks[i - 1].k) + 0.8 * Math.abs(x.e - blocks[i - 1].e) : 0));
    const inner = nov.slice(1).sort((p, q) => p - q), med = inner[Math.floor(inner.length / 2)] || 0;
    const thr = Math.max(0.25, Math.min(3 * med, med + 0.5));
    const cuts = [0];
    for (let i = 1; i < blocks.length; i++) {
      const drumIn = blocks[i].a === intro && intro > 0;
      const peak = nov[i] >= (nov[i - 1] || 0) * 0.9 && nov[i] >= (nov[i + 1] || 0) * 0.9;
      if (drumIn || (peak && nov[i] > thr)) {
        const last = cuts[cuts.length - 1];
        if (i - last >= 2 || i === blocks.length - 1 || last === 0 || drumIn) cuts.push(i);
        else if (nov[i] > nov[last]) cuts[cuts.length - 1] = i;
      }
    }
    const segs = cuts.map((c, j) => { const bEnd = j + 1 < cuts.length ? cuts[j + 1] : blocks.length; return { start: blocks[c].a, len: blocks[bEnd - 1].b - blocks[c].a, bi: [c, bEnd] }; });
    segs.forEach(sg => { sg.v = blockVec(sg.start, sg.start + sg.len); sg.e = sg.v[D] / 3; sg.k = sg.v[D + 1] / 3; });
    // agrupar partes que se parecen (repeticiones)
    let g = 0;
    segs.forEach((sg, i) => {
      for (let j = 0; j < i; j++) if (cos(sg.v.slice(0, D), segs[j].v.slice(0, D)) > 0.6) { sg.g = segs[j].g; return; }
      sg.g = g++;
    });
    // nombres: la parte repetida con más energía es el estribillo; la otra repetida, el verso
    const groups = {}; segs.forEach(sg => { (groups[sg.g] = groups[sg.g] || { g: sg.g, n: 0, e: 0, k: 0, len: 0 }); const G = groups[sg.g]; G.n++; G.e += sg.e; G.k += sg.k; G.len += sg.len; });
    const GL = Object.values(groups).map(G => ({ ...G, e: G.e / G.n, k: G.k / G.n, len: G.len / G.n, score: G.e / G.n + 0.5 * G.k / G.n }));
    const body = segs.filter(sg => !(sg.start === 0 && sg.k < 0.35 && intro > 0));
    const rep = GL.filter(G => G.n >= 2).sort((a, b) => b.score - a.score);
    let chorus = null, verse = null;
    if (rep.length >= 2) { chorus = rep[0].len <= rep[1].len ? rep[0].g : rep[1].g; const other = rep.find(G => G.g !== chorus); verse = other.g; if (rep[0].score - rep[1].score > 0.12) chorus = rep[0].g, verse = rep[1].g; }
    else if (rep.length === 1) { const G = rep[0]; const others = GL.filter(x => x.g !== G.g && x.len >= 8); if (others.length && others.some(o => o.score > G.score)) { verse = G.g; chorus = others.sort((a, b) => b.score - a.score)[0].g; } else verse = G.g; }
    segs.forEach((sg, i) => {
      if (i === 0 && intro > 0 && sg.start + sg.len <= intro + 1) sg.label = 'Intro';
      else if (sg.g === chorus) sg.label = 'Estribillo';
      else if (sg.g === verse) sg.label = 'Verso';
      else if (i === segs.length - 1 && sg.len <= 8) sg.label = 'Outro';
      else if (i === 0 && sg.len <= 8) sg.label = 'Intro';
      else if (sg.len <= 8 && sg.k < 0.5) sg.label = 'Puente';
      else if (!verse) { verse = sg.g; sg.label = 'Verso'; }
      else sg.label = sg.len <= 8 ? 'Puente' : 'Verso';
    });
    // si todo es lo mismo, es un loop
    if (segs.filter(sg => sg.label !== 'Intro' && sg.label !== 'Outro').every(sg => sg.g === segs.find(x => x.label !== 'Intro')?.g) && segs.length <= 3) segs.forEach(sg => { if (sg.label !== 'Intro' && sg.label !== 'Outro') sg.label = 'Loop'; });
    return segs.map(sg => ({ start: sg.start, len: sg.len, label: sg.label, g: sg.g, e: Math.round(sg.e * 100) / 100 }));
  }
  G.analyzeBeat = analyzeBeat;
})(typeof window !== 'undefined' ? window : globalThis);
