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
    for (let f = 0; f < nF; f++) {
      const o = f * H; let e = 0;
      for (let i = 0; i < N; i++) { const v = x[o + i]; re[i] = v * win[i]; im[i] = 0; e += v * v; }
      fft(re, im);
      for (let k = 0; k < N / 2; k++) L[k] = Math.log1p(100 * Math.hypot(re[k], im[k]));
      const fl = (a, b) => { let s = 0; for (let k = a; k <= b; k++) { const d = L[k] - Lp[k]; if (d > 0) s += d; } return s / (b - a + 1); };
      if (f > 0) { K[f] = fl(B.k[0], B.k[1]); S[f] = fl(B.s[0], B.s[1]) + 0.4 * fl(B.b[0], B.b[1]); Hh[f] = fl(B.h[0], B.h[1]); }
      let le = 0; for (let k = B.k[0]; k <= B.k[1]; k++) le += L[k]; LowE[f] = le;
      E[f] = Math.sqrt(e / N);
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
    // confianza: qué claro está la caja en 2 y 4
    const conf = Math.max(0, Math.min(1, (pat.s[4] + pat.s[12]) / 2));
    onProgress && onProgress(1);
    let off = phi / fr; const barS = barF / fr;
    if (off > barS - 0.08) { off = Math.max(0, off - barS); bars.unshift({ k: new Array(16).fill(0), s: new Array(16).fill(0) }); intro += 1; changes.forEach(c => c.bar++); }
    return { bpm, off, intro, bars, pat, changes, conf, dur: len / sr, v: 1 };
  }
  G.analyzeBeat = analyzeBeat;
})(typeof window !== 'undefined' ? window : globalThis);
