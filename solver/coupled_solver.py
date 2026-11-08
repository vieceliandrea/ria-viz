"""
coupled_solver.py
Assembla il vettore di stato globale y = [n, C_1..6, Tf_c(0..N-1), Tf_s(0..N-1), Tc(0..N-1)]
e integra nel tempo con solver stiff implicito (scipy BDF, stessa famiglia
algoritmica di DASSL/Petzold usato in TOPAZ: backward differentiation formula,
ordine variabile, passo adattivo, Jacobiano numerico).

Il problema e' un classico caso "stiff" multi-scala (TOPAZ VI, p.46-47):
- kinetics: scale temporali da Lambda/beta (~ms) a 1/lambda_i (~0.1-80 s)
- TH fuel: costante di tempo conduzione gap (~s)
- TH coolant: tempo di transito canale (~0.1-1 s)
Un solver esplicito richiederebbe dt << Lambda/beta ~ 1e-2 s anche per
seguire transitori da decine di secondi: improponibile. BDF implicito e'
la scelta corretta, esattamente come la scelta DASSL in TOPAZ.
"""

import numpy as np
from scipy.integrate import solve_ivp

from kinetics import PointKinetics, KineticsParams
from thermal_hydraulics import ChannelTH, ChannelGeom, FuelProps, CoolantProps
from reactivity_feedback import ReactivityModel, FeedbackParams, ExternalReactivity


class CoupledRIASystem:
    def __init__(self, kin: PointKinetics, th: ChannelTH, reac: ReactivityModel):
        self.kin = kin
        self.th = th
        self.reac = reac
        self.n_kin = 1 + kin.ngroups          # 7
        self.n_th = 3 * th.N                  # 3N
        self.n_tot = self.n_kin + self.n_th

    def pack_ic(self, n0: float = 1.0) -> np.ndarray:
        y_kin0 = self.kin.steady_state(n0)
        y_th0 = self.th.steady_state(n0)
        return np.concatenate((y_kin0, y_th0))

    def unpack(self, y: np.ndarray):
        y_kin = y[:self.n_kin]
        y_th = y[self.n_kin:]
        return y_kin, y_th

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        y_kin, y_th = self.unpack(y)
        n_rel = y_kin[0]

        Tf_avg = self.th.avg_fuel_T(y_th)
        Tc_avg = self.th.avg_coolant_T(y_th)
        rho = self.reac.total_rho(t, Tf_avg, Tc_avg)

        dydt_kin = self.kin.rhs(t, y_kin, rho)
        dydt_th = self.th.rhs(t, y_th, n_rel)
        return np.concatenate((dydt_kin, dydt_th))

    def solve(self, t_span, n0: float = 1.0, max_step=np.inf, rtol=1e-7, atol=1e-9):
        y0 = self.pack_ic(n0)

        # check di consistenza iniziale: dy/dt deve essere ~0 a t=0 se rho_ext(0)=0
        dydt0 = self.rhs(t_span[0], y0)
        res0 = np.max(np.abs(dydt0))

        sol = solve_ivp(
            fun=self.rhs,
            t_span=t_span,
            y0=y0,
            method="BDF",          # implicito, stiff -- analogo a DASSL (Petzold)
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=True,
        )
        if not sol.success:
            raise RuntimeError(f"Integrazione fallita: {sol.message}")
        return sol, res0


def build_default_system(rho_max=0.006, t_ramp=0.1, W_flow=0.30,
                          P0_total=5.0e4, Tin=559.15):
    """
    Factory per il caso di default. Tin in Kelvin.
    W_flow, P0_total dimensionati per canale singolo equivalente con
    salita di temperatura realistica (~30-40 K) a piena potenza.
    """
    kin = PointKinetics(KineticsParams())

    geom = ChannelGeom(N=10, L=3.6)
    fuel = FuelProps()
    cool = CoolantProps()
    th = ChannelTH(geom, fuel, cool, W_flow=W_flow, Tin=Tin, P0_total=P0_total)

    yss_th = th.steady_state(1.0)
    Tf_ref = th.avg_fuel_T(yss_th)
    Tc_ref = th.avg_coolant_T(yss_th)

    fb = FeedbackParams(alpha_D=-3.0e-5, alpha_m=-2.0e-4, Tf_ref=Tf_ref, Tc_ref=Tc_ref)
    ext = ExternalReactivity(kind="ramp", rho_max=rho_max, t_start=0.0, t_ramp=t_ramp)
    reac = ReactivityModel(fb, ext)

    return CoupledRIASystem(kin, th, reac)


if __name__ == "__main__":
    sysm = build_default_system()
    sol, res0 = sysm.solve(t_span=(0.0, 5.0), n0=1.0)
    print(f"Residuo |dy/dt| a t=0 (deve essere ~0): {res0:.3e}")
    print(f"n_relativo finale (t={sol.t[-1]:.2f}s): {sol.y[0,-1]:.4f}")
    print(f"n_relativo massimo: {np.max(sol.y[0,:]):.4f} a t={sol.t[np.argmax(sol.y[0,:])]:.3f}s")
