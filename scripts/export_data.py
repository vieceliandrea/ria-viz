"""Run the RIA solver (./solver) and export a compact JSON for the web app -> data/ria.json
Usage: python scripts/export_data.py --rho-pcm 600 --t-ramp 0.1 --t-end 20"""
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "solver"))
from main_ria import run_case  # noqa: E402

# Loop transport delays [s]: ASSUMED, visualization only (the solver has a fixed inlet temperature)
TAU_HL, TAU_HX, TAU_COLD = 1.5, 3.0, 2.5

ap = argparse.ArgumentParser()
ap.add_argument("--rho-pcm", type=float, default=600.0)  # external reactivity [pcm]
ap.add_argument("--t-ramp", type=float, default=0.1)     # ramp duration [s]
ap.add_argument("--t-end", type=float, default=20.0)     # [s]
a = ap.parse_args()

r = run_case(rho_max_pcm=a.rho_pcm, t_ramp=a.t_ramp, t_end=a.t_end, n_eval=int(a.t_end * 100) + 1)
th, kin, t = r["sysm"].th, r["sysm"].kin, r["t"]
idx = [i for i in range(len(t)) if i <= 100 or i % 5 == 0 or i == len(t) - 1]  # dense early, coarse late
rd = lambda x, d: [round(float(v), d) for v in np.asarray(x)[idx]]
Tc = r["Tc"] - 273.15                                    # [C], shape (N cells, time)
v = th.W / (th.c.rho_c * th.g.flow_area)                 # [m/s] coolant velocity

out = dict(
    t=rd(t, 3), n=rd(r["n_rel"], 4), rho=rd(r["rho_tot"] * 1e5, 2), tf=rd(r["Tf_avg"] - 273.15, 2),
    tc=[rd(Tc[j], 2) for j in range(th.N)],
    meta=dict(
        Tin=th.Tin - 273.15, Tout0=float(Tc[-1, 0]), Tout_peak=float(Tc[-1].max()),
        Tf0=float(r["Tf_avg"][0] - 273.15), Tf_peak=float(r["Tf_avg"].max() - 273.15),
        n_peak=float(r["n_rel"].max()), t_peak=float(t[r["n_rel"].argmax()]), n_end=float(r["n_rel"][-1]),
        rho_end=float(r["rho_tot"][-1] * 1e5), rho_ext=a.rho_pcm, beta=float(kin.p.beta * 1e5),
        t_ramp=a.t_ramp, t_end=a.t_end, tau_core=float(th.g.L / v),
        tau_hl=TAU_HL, tau_hx=TAU_HX, tau_cold=TAU_COLD,
    ),
)
path = ROOT / "data" / "ria.json"
path.write_text(json.dumps(out, separators=(",", ":")))
print(f"wrote {path} ({path.stat().st_size // 1024} kB), n_peak={out['meta']['n_peak']:.1f}")
