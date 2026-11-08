# RIA transient: primary circuit visualization

Next.js (App Router) front end for a reactivity-initiated-accident transient in a water cooled nuclear reactor computed by a Python
point-kinetics + 1-D thermal-hydraulics solver. Colour = moderator water temperature.

```bash
npm install && npm run dev                                              # http://localhost:3000
python scripts/export_data.py --rho-pcm 600 --t-ramp 0.1 --t-end 20     # regenerate data/ria.json (numpy, scipy)
```

Layout: `solver/` physics, `scripts/export_data.py` solver to JSON, `data/ria.json`, `components/RiaViz.jsx`, `lib/helpers.js`.

## Limitations
Only the core is solved (1-D channel, constant flow and inlet T, single phase). Hot leg and exchanger are delayed
plug flow of the core outlet (assumed delays) and the cold leg is fixed at T_in: the loop is not closed, so results
after about one loop transit (~8 s) are weak. Point kinetics, low beta, illustrative feedback coefficients, lumped
fuel: qualitative only, not a safety analysis.

Deploy: import the GitHub repo on Vercel (zero config).
