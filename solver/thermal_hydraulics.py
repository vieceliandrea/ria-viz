"""
thermal_hydraulics.py
Canale 1-D, mesh control-volume assiale (filosofia TOPAZ: control volumes per
continuity/energy, stream tubes per momentum -- qui semplificato a portata
imposta, niente equazione di momento risolta esplicitamente: W(t) e' un
input, valido per transitori RIA a portata di pompa costante o con legge nota).

Per ogni nodo assiale k = 0..N-1:
  - Fuel: 2 nodi radiali lumped (centro Tf_c, superficie Tf_s) con
    conduttanza equivalente UA_fuel (gap+pellet semplificato a singola R).
  - Coolant: 1 nodo per cella, energia con upwind differencing (identico
    allo schema upwind descritto in TOPAZ III-7..III-9), continuita'
    incompressibile (portata W costante lungo il canale, no aree variabili).

Variabili di stato per nodo: [Tf_c, Tf_s, Tc]  -> 3*N stati totali.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class ChannelGeom:
    N: int = 10                  # numero di celle assiali
    L: float = 3.6                # [m] lunghezza attiva
    rod_pitch: float = 1.26e-2    # [m]
    rod_diam: float = 9.5e-3      # [m] diametro pastiglia+clad esterno
    clad_thick: float = 0.57e-3   # [m]
    n_rods_per_channel: int = 1   # canale equivalente per 1 barretta (scalabile)

    @property
    def dz(self) -> float:
        return self.L / self.N

    @property
    def flow_area(self) -> float:
        # area di flusso equivalente per cella unitaria a passo quadro - area barretta
        return self.rod_pitch**2 - np.pi * (self.rod_diam / 2) ** 2

    @property
    def heated_perimeter(self) -> float:
        return np.pi * self.rod_diam

    @property
    def Dh(self) -> float:
        # diametro idraulico equivalente (definizione 4A/P)
        return 4.0 * self.flow_area / self.heated_perimeter


@dataclass
class FuelProps:
    rho_f: float = 10970.0    # [kg/m3] UO2
    cp_f: float = 300.0       # [J/kg/K] media operativa
    k_gap: float = 5680.0     # [W/m2/K] conduttanza gap+pellet equivalente (lumped)
    rho_c_clad: float = 6560.0   # [kg/m3] zircaloy
    cp_clad: float = 330.0       # [J/kg/K]


@dataclass
class CoolantProps:
    rho_c: float = 720.0      # [kg/m3] acqua sottoraffreddata PWR tipica (~300 C, 15.5 MPa)
    cp_c: float = 5500.0      # [J/kg/K]
    k_c: float = 0.56         # [W/m/K]
    mu_c: float = 9.0e-5      # [Pa s]

    @property
    def Pr_c(self) -> float:
        return self.mu_c * self.cp_c / self.k_c


class ChannelTH:
    """
    Stato y_th = [Tf_c_0..N-1, Tf_s_0..N-1, Tc_0..N-1]  (3N variabili)
    """

    def __init__(self, geom: ChannelGeom, fuel: FuelProps, coolant: CoolantProps,
                 W_flow: float, Tin: float, P0_total: float):
        self.g = geom
        self.f = fuel
        self.c = coolant
        self.W = W_flow            # [kg/s] portata nel canale (costante o funzione di t esterna)
        self.Tin = Tin             # [K] temperatura ingresso
        self.N = geom.N

        # potenza nominale per cella: distribuzione assiale tipo coseno troncato
        z = (np.arange(self.N) + 0.5) * geom.dz - geom.L / 2
        shape = np.cos(np.pi * z / (1.3 * geom.L))   # forma normalizzata
        shape = np.clip(shape, 0, None)
        shape /= np.sum(shape) * geom.dz
        self.q_shape = shape           # [1/m], integra a 1 su L
        self.P0_total = P0_total       # [W] potenza nominale canale

        # masse/capacita' per nodo (per unita' di lunghezza * dz)
        r_pellet = (geom.rod_diam - 2 * geom.clad_thick) / 2
        A_pellet = np.pi * r_pellet**2
        A_clad = np.pi * (geom.rod_diam / 2)**2 - A_pellet

        self.M_fc = fuel.rho_f * A_pellet * geom.dz   # massa "centro" pellet (lumped 1 nodo equivalente)
        self.M_fs = fuel.rho_c_clad * A_clad * geom.dz  # massa clad (nodo superficie)
        self.M_c_cell = coolant.rho_c * geom.flow_area * geom.dz

        self.UA_gap = fuel.k_gap * geom.heated_perimeter * geom.dz  # [W/K] conduttanza pellet->clad lumped

    # ---------------- correlazioni costitutive (quasi-steady, stile TOPAZ III-22/23) ----------------
    def htc_dittus_boelter(self, Tc: np.ndarray) -> np.ndarray:
        """h = 0.023 Re^0.8 Pr^0.4 k/Dh  -- riscaldamento (Dittus-Boelter)."""
        v = self.W / (self.c.rho_c * self.g.flow_area)
        Re = self.c.rho_c * v * self.g.Dh / self.c.mu_c
        Pr = self.c.Pr_c
        Nu = 0.023 * Re**0.8 * Pr**0.4
        h = Nu * self.c.k_c / self.g.Dh
        return np.full_like(Tc, h)   # proprieta' costanti -> h uniforme (estendibile a f(T))

    # ---------------- stato stazionario ----------------
    def steady_state(self, n_rel: float = 1.0) -> np.ndarray:
        """Risolve lo stazionario per condizione di potenza relativa n_rel (=1 nominale)."""
        N = self.N
        q_cell = n_rel * self.P0_total * self.q_shape * self.g.dz   # [W] per cella

        Tc = np.zeros(N)
        Tf_s = np.zeros(N)
        Tf_c = np.zeros(N)
        h = self.htc_dittus_boelter(np.zeros(N))[0]
        Tprev = self.Tin
        for k in range(N):
            # energia coolant: W*cp*(Tc_k - Tprev) = q_cell_k (stazionario, upwind)
            Tc[k] = Tprev + q_cell[k] / (self.W * self.c.cp_c)
            Tprev = Tc[k]
            As = self.g.heated_perimeter * self.g.dz
            Tf_s[k] = Tc[k] + q_cell[k] / (h * As)
            Tf_c[k] = Tf_s[k] + q_cell[k] / self.UA_gap
        return np.concatenate((Tf_c, Tf_s, Tc))

    # ---------------- derivate per integrazione transitoria ----------------
    def rhs(self, t: float, y_th: np.ndarray, n_rel: float) -> np.ndarray:
        N = self.N
        Tf_c = y_th[0:N]
        Tf_s = y_th[N:2 * N]
        Tc = y_th[2 * N:3 * N]

        q_cell = n_rel * self.P0_total * self.q_shape * self.g.dz   # [W] potenza istantanea per cella

        h = self.htc_dittus_boelter(Tc)
        As = self.g.heated_perimeter * self.g.dz

        q_gap = self.UA_gap * (Tf_c - Tf_s)              # pellet -> clad
        q_conv = h * As * (Tf_s - Tc)                      # clad -> coolant

        dTfc_dt = (q_cell - q_gap) / (self.M_fc * self.f.cp_f)
        dTfs_dt = (q_gap - q_conv) / (self.M_fs * self.f.cp_clad)

        Tc_upwind_in = np.empty(N)
        Tc_upwind_in[0] = self.Tin
        Tc_upwind_in[1:] = Tc[:-1]
        dTc_dt = (q_conv + self.W * self.c.cp_c * (Tc_upwind_in - Tc)) / (self.M_c_cell * self.c.cp_c)

        return np.concatenate((dTfc_dt, dTfs_dt, dTc_dt))

    # ---------------- quantita' derivate utili al coupling ----------------
    def avg_fuel_T(self, y_th: np.ndarray) -> float:
        N = self.N
        Tf_c = y_th[0:N]
        Tf_s = y_th[N:2 * N]
        w = self.q_shape  # pesatura sulla potenza locale
        w = w / np.sum(w)
        Tf_avg_node = 0.6 * Tf_c + 0.4 * Tf_s   # peso maggiore al centro (gradiente radiale tipico)
        return float(np.sum(w * Tf_avg_node))

    def avg_coolant_T(self, y_th: np.ndarray) -> float:
        N = self.N
        Tc = y_th[2 * N:3 * N]
        w = self.q_shape / np.sum(self.q_shape)
        return float(np.sum(w * Tc))


if __name__ == "__main__":
    geom = ChannelGeom()
    fuel = FuelProps()
    cool = CoolantProps()
    print(f"Pr_c (check) = {cool.Pr_c:.3f}  (atteso ~0.85-0.9 per PWR a piena T)")
    print(f"Dh = {geom.Dh*1e3:.3f} mm, flow_area = {geom.flow_area*1e6:.2f} mm2")

    ch = ChannelTH(geom, fuel, cool, W_flow=0.30, Tin=559.0, P0_total=2.0e5)
    yss = ch.steady_state(1.0)
    dydt = ch.rhs(0.0, yss, n_rel=1.0)
    print("max |dy/dt| allo stazionario:", np.max(np.abs(dydt)))
    print("T uscita coolant stazionaria:", yss[2*geom.N+geom.N-1] - 273.15, "C")
