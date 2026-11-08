'use client';
import { memo, useEffect, useRef, useState } from 'react';
import D from '../data/ria.json';
import { at, f, gradientCss, makeColor } from '../lib/helpers';

const { t: T, meta: M } = D;
const TE = T[T.length - 1];                       // [s] end time
const OUT = D.tc[D.tc.length - 1];                // [C] core outlet cell
const NH = 14, Y0 = 108, Y1 = 452, DY = (Y1 - Y0) / 10, WH = 820 / NH, DUR = 26000; // geometry [svg units], playback [ms]
const LO = Math.floor(M.Tin) - 1, HI = Math.ceil(M.Tout_peak) + 2;
const col = makeColor(LO, HI), COLD = col(M.Tin);
const LOOP = M.tau_core + M.tau_hl + M.tau_hx + M.tau_cold;

const PH = [
  [0, 'Reactivity ramp', `External reactivity rises to +${f(M.rho_ext)} pcm in ${f(M.t_ramp, 1)} s, i.e. ${f(M.rho_ext / M.beta, 2)} β: power takes off within a few tenths of a second.`],
  [M.t_ramp, 'Prompt peak', `Power reaches ×${f(M.n_peak)} at ${f(M.t_peak, 2)} s. Fuel goes from ${f(M.Tf0)} to ${f(M.Tf_peak)} °C and Doppler feedback curbs the rise. The water has not reacted yet: heat must first cross pellet and cladding.`],
  [0.5, 'Heat reaches the water', `The hot wave climbs the core (transit ≈ ${f(M.tau_core, 1)} s), then travels through hot leg and heat exchanger. Core outlet peaks at ${f(M.Tout_peak)} °C, up from ${f(M.Tout0)} °C.`],
  [3, 'Towards a new equilibrium', `At ${f(M.t_end)} s: power ×${f(M.n_end, 2)}, net reactivity ${f(M.rho_end)} pcm. The external reactivity stays inserted, so power does not return to ×1.`],
];
const LIM = [
  'Only the core is solved: one 1-D channel, 10 axial cells, constant flow and inlet temperature, single-phase water (no boiling, no momentum equation).',
  `The loop is not closed: hot leg and exchanger are delayed plug flow of the core outlet (τ = ${f(M.tau_hl, 1)} s and ${f(M.tau_hx, 1)} s, assumed), the cold leg is held at T_in. After one loop (≈ ${f(LOOP)} s) the return of hot water is missing, so late-time results are weak.`,
  `Point kinetics with 6 delayed groups, no spatial effects. β = ${f(M.beta)} pcm is low for UO₂ and the feedback coefficients are illustrative: the peak is very sensitive to both.`,
  'Lumped fuel (2 nodes per cell), constant properties, no gap or burnup effects, no failure criteria: qualitative, not a safety analysis.',
  'Playback uses a square-root time scale so the prompt peak stays visible.',
];

// ---- charts (time axis in sqrt scale, same as the slider) ----
const C = { W: 500, H: 220, L: 40, R: 10, T: 28, B: 26 };
const X = (t) => C.L + (C.W - C.L - C.R) * Math.sqrt(t / TE);
const PMAX = 10 ** Math.ceil(Math.log10(M.n_peak));
const TMIN = Math.floor(M.Tin / 20) * 20, TMAX = Math.ceil(M.Tf_peak / 20) * 20;
const range = (a, b, g) => Array.from({ length: Math.floor((b - a) / g) + 1 }, (_, i) => a + i * g);
const S1 = [{ a: D.n, c: '#ff4fa3' }];
const S2 = [{ a: D.tf, c: '#ffb454', n: 'fuel' }, { a: OUT, c: '#4dd6e0', n: 'core outlet' }];
const K1 = range(0, Math.log10(PMAX), 1).map((e) => 10 ** e);
const K2 = range(Math.ceil(TMIN / 100) * 100, TMAX, 100);

const Static = memo(function Static({ title, series, y0, y1, log, ticks }) {
  const Y = (v) => {
    const u = log ? Math.log10(v / y0) / Math.log10(y1 / y0) : (v - y0) / (y1 - y0);
    return C.H - C.B - (C.H - C.B - C.T) * u;
  };
  return (
    <>
      <text x={C.L} y="16" className="ct">{title}</text>
      {ticks.map((v) => (
        <g key={v}>
          <line x1={C.L} x2={C.W - C.R} y1={Y(v)} y2={Y(v)} className="gr" />
          <text x={C.L - 6} y={Y(v) + 4} className="ax" textAnchor="end">{v}</text>
        </g>
      ))}
      {[0, 0.1, 0.5, 2, 5, TE].map((t) => (
        <text key={t} x={X(t)} y={C.H - 8} className="ax" textAnchor="middle">{t}{t === TE ? ' s' : ''}</text>
      ))}
      {series.map((q, i) => (
        <g key={i}>
          <path d={'M' + T.map((t, k) => `${X(t).toFixed(1)} ${Y(q.a[k]).toFixed(1)}`).join('L')} fill="none" stroke={q.c} strokeWidth="2" />
          {q.n && <text x={150 + i * 110} y="16" className="ct" style={{ fill: q.c }}>{q.n}</text>}
        </g>
      ))}
    </>
  );
});
const Chart = ({ t, ...p }) => (
  <svg viewBox={`0 0 ${C.W} ${C.H}`}>
    <Static {...p} />
    <line x1={X(t)} x2={X(t)} y1={C.T} y2={C.H - C.B} className="cur" />
  </svg>
);

export default function RiaViz() {
  const [s, setS] = useState(0);                  // slider position 0..1000 (t = TE (s/1000)^2)
  const [playing, setPlaying] = useState(false);
  const sRef = useRef(0);
  const t = TE * (s / 1000) ** 2;

  useEffect(() => {
    if (!matchMedia('(prefers-reduced-motion: reduce)').matches) setPlaying(true);
  }, []);
  useEffect(() => {
    if (!playing) return;
    let id, last = performance.now();
    const tick = (now) => {
      sRef.current = Math.min(1000, sRef.current + ((now - last) * 1000) / DUR);
      last = now;
      setS(sRef.current);
      if (sRef.current >= 1000) setPlaying(false); else id = requestAnimationFrame(tick);
    };
    id = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(id);
  }, [playing]);

  const toggle = () => {
    if (!playing && sRef.current >= 1000) { sRef.current = 0; setS(0); }
    setPlaying((p) => !p);
  };
  const seek = (e) => { sRef.current = +e.target.value; setS(sRef.current); };

  const n = at(T, D.n, t);
  const ph = PH.filter((q) => t >= q[0]).pop();

  return (
    <main>
      <h1>A RIA seen from the primary circuit</h1>
      <p className="sub">Reactivity-initiated accident in a point-kinetics + 1-D thermal-hydraulics model. Colour is the moderator water temperature; the pink halo around the core is power.</p>
      <div className="stage">
        <svg viewBox="60 40 920 480" role="img" aria-label="2D primary circuit: core, hot leg, heat exchanger and cold leg coloured by moderator water temperature">
          <defs><linearGradient id="sec" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stopColor="#7f8bff" /><stop offset="1" stopColor="#dfe9ff" /></linearGradient></defs>
          <rect x="84" y="102" width="132" height="356" rx="6" fill="none" stroke="#ff4fa3" strokeWidth="5"
            style={{ filter: 'drop-shadow(0 0 12px #ff4fa3)', opacity: Math.min(1, Math.log10(Math.max(n, 1)) / Math.log10(M.n_peak)) }} />
          {D.tc.map((a, j) => <rect key={'c' + j} x="90" y={Y1 - (j + 1) * DY} width="120" height={DY + 0.5} fill={col(at(T, a, t))} />)}
          {D.tc.map((_, i) => {
            const x = (i + 0.5) / 10;
            return <rect key={'x' + i} x="790" y={Y0 + i * DY} width="120" height={DY + 0.5}
              fill={col(M.Tin + (at(T, OUT, t - M.tau_hl - M.tau_hx * x) - M.Tin) * (1 - x))} />;
          })}
          {Array.from({ length: NH }, (_, m) => (
            <g key={m}>
              <rect x={90 + m * WH} y="72" width={WH + 0.5} height="36" fill={col(at(T, OUT, t - (M.tau_hl * (m + 0.5)) / NH))} />
              <rect x={90 + m * WH} y="452" width={WH + 0.5} height="36" fill={COLD} />
            </g>
          ))}
          {[0, 1, 2].map((k) => <rect key={k} x={121 + k * 24} y={Y0} width="5" height={Y1 - Y0} rx="2" fill="rgba(14,11,36,.5)" />)}
          <rect x="916" y="108" width="12" height="344" rx="3" fill="url(#sec)" />
          <text className="lb" x="954" y="280" transform="rotate(90 954 280)" textAnchor="middle">secondary side</text>
          <path className="flow" d="M150 470V90H850V470H150Z" />
          <text className="lb" x="500" y="140" textAnchor="middle">Hot leg</text>
          <text className="lb" x="500" y="436" textAnchor="middle">Cold leg</text>
          <text className="lb" x="226" y="290">Core</text>
          <text className="lb" x="774" y="290" textAnchor="end">Heat exchanger</text>
          <text className="lb" x="500" y="222" textAnchor="middle">relative power n/n₀</text>
          <text className="big" x="500" y="306" textAnchor="middle">×{f(n, n < 10 ? 2 : 1)}</text>
          <text className="lb" x="500" y="344" textAnchor="middle">t = {f(t, 3)} s</text>
        </svg>
      </div>
      <div className="ctl">
        <button onClick={toggle}>{playing ? 'Pause' : s >= 1000 ? 'Replay' : 'Play'}</button>
        <input type="range" min="0" max="1000" value={s} onChange={seek} aria-label="Transient time" />
        <output>{f(t, 2)} s</output>
      </div>
      <div className="row">
        <div className="cap"><h2>{ph[1]}</h2><p>{ph[2]}</p></div>
        <div className="ro">
          <div><b>{f(at(T, D.rho, t))}</b><span>net reactivity [pcm]</span></div>
          <div><b>{f(at(T, D.tf, t))}</b><span>mean fuel T [°C]</span></div>
          <div><b>{f(at(T, OUT, t))}</b><span>core outlet T [°C]</span></div>
        </div>
      </div>
      <div className="lg"><span>{LO}</span><i style={{ background: gradientCss() }} /><span>{HI}</span><span>water temperature [°C]</span></div>
      <div className="charts">
        <Chart t={t} title="Power n/n₀" series={S1} y0={1} y1={PMAX} log ticks={K1} />
        <Chart t={t} title="Temperatures [°C]" series={S2} y0={TMIN} y1={TMAX} ticks={K2} />
      </div>
      <section className="lim"><h2>Limitations</h2><ul>{LIM.map((l) => <li key={l}>{l}</li>)}</ul></section>
    </main>
  );
}
