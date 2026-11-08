export const f = (x, d = 0) => x.toFixed(d);

// Temperature -> colour ramp (cold blue -> hot magenta), piecewise linear
const STOPS = [[0, [47, 107, 255]], [0.25, [34, 211, 197]], [0.5, [255, 210, 63]], [0.75, [255, 122, 47]], [1, [255, 45, 111]]];
export const makeColor = (lo, hi) => (v) => {
  const x = Math.min(1, Math.max(0, (v - lo) / (hi - lo)));
  let i = 0;
  while (i < 3 && x > STOPS[i + 1][0]) i++;
  const [x0, c0] = STOPS[i], [x1, c1] = STOPS[i + 1], u = (x - x0) / (x1 - x0);
  return `rgb(${c0.map((c, k) => Math.round(c + (c1[k] - c) * u)).join()})`;
};
export const gradientCss = () => `linear-gradient(90deg,${STOPS.map(([p, c]) => `rgb(${c}) ${p * 100}%`).join()})`;

// Series `a` sampled on grid T, linearly interpolated at time t (clamped at both ends)
export const at = (T, a, t) => {
  const n = T.length - 1;
  if (t <= 0) return a[0];
  if (t >= T[n]) return a[n];
  let l = 0, h = n;
  while (h - l > 1) { const m = (l + h) >> 1; if (T[m] > t) h = m; else l = m; }
  return a[l] + ((a[h] - a[l]) * (t - T[l])) / (T[h] - T[l]);
};
