"""
kinetics.py
Point Reactor Kinetics Equations (PRKE), 6 gruppi di precursori.

dn/dt   = (rho - beta)/Lambda * n + sum_i lambda_i * C_i
dC_i/dt = (beta_i/Lambda) * n - lambda_i * C_i

n          : densita' neutronica relativa (n=1 a potenza nominale)
C_i        : concentrazione precursori gruppo i, normalizzata come n
rho        : reattivita' totale [-] (delta-k/k), funzione di stato e tempo
Lambda     : prompt neutron generation time [s]
beta_i     : delayed neutron fraction gruppo i [-]
lambda_i   : decay constant gruppo i [1/s]

Dati di default: U-235 termico (Keepin), valori standard di letteratura.
"""

import numpy as np
from dataclasses import dataclass, field


@dataclass
class KineticsParams:
    Lambda: float = 2.0e-5          # [s] prompt generation time (PWR tipico ~1e-5..1e-4 s)
    beta_i: np.ndarray = field(default_factory=lambda: np.array(
        [2.18e-4, 1.024e-3, 6.05e-4, 1.31e-3, 2.20e-4, 6.0e-5]))   # U-235 (Keepin 1965)
    lambda_i: np.ndarray = field(default_factory=lambda: np.array(
        [1.249e-2, 3.013e-2, 1.116e-1, 3.013e-1, 1.140e0, 3.014e0]))  # [1/s]

    @property
    def beta(self) -> float:
        return float(np.sum(self.beta_i))


class PointKinetics:
    """
    Stato: y = [n, C_1, ..., C_6]   (lunghezza 7)
    Condizione iniziale stazionaria: dC_i/dt = 0 => C_i0 = (beta_i/(lambda_i*Lambda)) * n0
    """

    def __init__(self, params: KineticsParams = None):
        self.p = params if params is not None else KineticsParams()
        self.ngroups = len(self.p.beta_i)

    def steady_state(self, n0: float = 1.0) -> np.ndarray:
        C0 = (self.p.beta_i / (self.p.lambda_i * self.p.Lambda)) * n0
        return np.concatenate(([n0], C0))

    def rhs(self, t: float, y: np.ndarray, rho: float) -> np.ndarray:
        """
        Residuo dy/dt dato rho(t) GIA' VALUTATO (incluso feedback) dal chiamante.
        Tenuto separato da reactivity feedback per modularita' (coupling in main).
        """
        n = y[0]
        C = y[1:]
        beta = self.p.beta
        dndt = (rho - beta) / self.p.Lambda * n + np.sum(self.p.lambda_i * C)
        dCdt = (self.p.beta_i / self.p.Lambda) * n - self.p.lambda_i * C
        return np.concatenate(([dndt], dCdt))

    # --- sanity check di modulo, eseguito solo se invocato direttamente ---
    def _check_steady(self):
        y0 = self.steady_state(1.0)
        dydt = self.rhs(0.0, y0, rho=0.0)
        assert np.allclose(dydt, 0.0, atol=1e-12), "Stato stazionario non verificato!"
        print("Check OK: rho=0, y=y_ss => dy/dt = 0 (entro 1e-12)")


if __name__ == "__main__":
    k = PointKinetics()
    print(f"beta totale = {k.p.beta:.6e}")
    k._check_steady()
