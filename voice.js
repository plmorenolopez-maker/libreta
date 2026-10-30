'use strict';
// ───── Golpes de voz: dónde empieza cada sílaba en una toma ─────
// Entra audio mono y su frecuencia de muestreo. Sale una lista de { t (s), v (fuerza 0..1) }.
// Busca subidas bruscas de energía en la banda de la voz (250–3500 Hz): cada sílaba rapeada da una.
(function (G) {
  function fftR(re, im) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
    for (let len = 2; len <= n; len <<= 1) {
      const a = -2 * Math.PI / len, wr = Math.cos(a), wi = Math.sin(a), h = len >> 1;
      for (let i = 0; i < n; i += len) { let cr = 1, ci = 0; for (let k = 0; k < h; k++) { const p = i + k, q = p + h, xr = re[q] * cr - im[q] * ci, xi = re[q] * ci + im[q] * cr; re[q] = re[p] - xr; im[q] = im[p] - xi; re[p] += xr; im[p] += xi; const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t; } }
    }
  }
  function detectOnsets(x, sr, opts = {}) {
    const H = Math.round(sr * 0.005), N = 1024, nF = Math.max(0, Math.floor((x.length - N) / H));
    const win = new Float32Array(N); for (let i = 0; i < N; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / N);
    const k0 = Math.round(250 * N / sr), k1 = Math.round(3500 * N / sr), kc = Math.round(1500 * N / sr);
    const E = new Float32Array(nF), F = new Float32Array(nF), re = new Float32Array(N), im = new Float32Array(N), prev = new Float32Array(N / 2), cur = new Float32Array(N / 2);
    for (let f = 0; f < nF; f++) {
      const o = f * H;
      for (let i = 0; i < N; i++) { re[i] = x[o + i] * win[i]; im[i] = 0; }
      fftR(re, im);
      let e = 0, fl = 0;
      for (let k = k0; k <= k1; k++) { const m = re[k] * re[k] + im[k] * im[k]; e += m; cur[k] = Math.log1p(1000 * Math.sqrt(m)); const d = cur[k] - prev[k]; if (d > 0) fl += d * (k > kc ? 1.3 : 1); }
      E[f] = 10 * Math.log10(e + 1e-12); F[f] = f ? fl / (k1 - k0 + 1) : 0;
      prev.set(cur);
    }
    // suavizado ligero
    const Fs = new Float32Array(nF); for (let f = 0; f < nF; f++) { let s = 0, c = 0; for (let j = -2; j <= 2; j++) if (f + j >= 0 && f + j < nF) { s += F[f + j] * (j === 0 ? 2 : 1); c += j === 0 ? 2 : 1; } Fs[f] = s / c; }
    const sorted = a => Array.from(a).sort((p, q) => p - q);
    const eS = sorted(E), floor = eS[Math.floor(eS.length * 0.15)] || -120, loud = eS[Math.floor(eS.length * 0.97)] || 0;
    const fS = sorted(Fs), fRef = fS[Math.floor(fS.length * 0.98)] || 1;
    const gate = floor + Math.max(10, (loud - floor) * 0.35), minGap = Math.round((opts.minGap || 0.085) / 0.005), look = Math.round(0.04 / 0.005);
    const out = []; let last = -1e9;
    for (let f = 1; f < nF - 1; f++) {
      if (Fs[f] < 0.22 * fRef) continue;
      let peak = true; for (let j = Math.max(0, f - look); j <= Math.min(nF - 1, f + look); j++) if (Fs[j] > Fs[f]) { peak = false; break; }
      if (!peak) continue;
      // tiene que sonar de verdad poco después (no es ruido)
      let mx = -1e9; for (let j = f; j < Math.min(nF, f + 16); j++) if (E[j] > mx) mx = E[j];
      if (mx < gate) continue;
      if (f - last < minGap) { if (out.length && Fs[f] > out[out.length - 1].raw) { out[out.length - 1] = { t: f * H / sr, raw: Fs[f] }; last = f; } continue; }
      // el inicio está un poco antes del pico de la subida
      let s = f; while (s > 0 && Fs[s - 1] > 0.5 * Fs[f] && f - s < 6) s--;
      out.push({ t: (s * H + N / 2) / sr + 0.008, raw: Fs[f] }); last = f;
    }
    const mxv = Math.max(...out.map(o => o.raw), 1e-9);
    return out.map(o => ({ t: Math.round(o.t * 1000) / 1000, v: Math.round(o.raw / mxv * 100) / 100 }));
  }
  G.detectOnsets = detectOnsets;
})(typeof window !== 'undefined' ? window : globalThis);
