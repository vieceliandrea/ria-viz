"""
main_ria.py
Caso applicativo: RIA (Reactivity Initiated Accident) - inserimento rampa di
reattivita' (es. espulsione/eiezione parziale di una barra di controllo),
con coupling point-kinetics <-> termoidraulica 1-D canale caldo.

Esegue:
  1. Costruzione sistema (kinetics + TH + feedback)
  2. Integrazione stiff (BDF, vedi coupled_solver.py)
  3. Sanity checks fisici (Sezione "Check" obbligatoria)
  4. Plot: potenza, reattivita' nette, temperature fuel/coolant
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from coupled_solver import build_default_system


def run_case(rho_max_pcm=600.0, t_ramp=0.1, t_end=10.0, n_eval=2000):
    """rho_max_pcm in pcm (1 pcm = 1e-5 delta-k/k)."""
    rho_max = rho_max_pcm * 1e-5
    sysm = build_default_system(rho_max=rho_max, t_ramp=t_ramp)

    sol, res0 = sysm.solve(t_span=(0.0, t_end), n0=1.0, rtol=1e-8, atol=1e-10)

    t_eval = np.linspace(0.0, t_end, n_eval)
    Y = sol.sol(t_eval)

    n_kin = sysm.n_kin
    n_rel = Y[0, :]
    C = Y[1:n_kin, :]

    N = sysm.th.N
    Tf_c = Y[n_kin:n_kin + N, :]
    Tf_s = Y[n_kin + N:n_kin + 2 * N, :]
    Tc = Y[n_kin + 2 * N:n_kin + 3 * N, :]

    Tf_avg = np.array([sysm.th.avg_fuel_T(Y[n_kin:, i]) for i in range(len(t_eval))])
    Tc_avg = np.array([sysm.th.avg_coolant_T(Y[n_kin:, i]) for i in range(len(t_eval))])
    rho_tot = np.array([sysm.reac.total_rho(t_eval[i], Tf_avg[i], Tc_avg[i]) for i in range(len(t_eval))])
    rho_ext = np.array([sysm.reac.ext(ti) for ti in t_eval])

    return dict(sysm=sysm, t=t_eval, n_rel=n_rel, C=C, Tf_c=Tf_c, Tf_s=Tf_s, Tc=Tc,
                Tf_avg=Tf_avg, Tc_avg=Tc_avg, rho_tot=rho_tot, rho_ext=rho_ext, res0=res0)


def run_checks(res):
    """[Check] obbligatorio: consistenza iniziale, limiti fisici, ordini di grandezza."""
    print("=" * 70)
    print("[CHECK] Sanity checks fisici")
    print("=" * 70)

    print(f"1) Residuo |dy/dt| a t=0 (consistenza IC): {res['res0']:.3e}  -> deve essere ~1e-10..1e-12")
    assert res['res0'] < 1e-6, "IC non consistente!"

    n_max = np.max(res['n_rel'])
    t_peak = res['t'][np.argmax(res['n_rel'])]
    print(f"2) Picco di potenza relativa: {n_max:.3f} x nominale, a t={t_peak:.3f} s")
    print("   -> per RIA con feedback Doppler negativo: picco prompt poi turn-over (atteso)")

    n_final = res['n_rel'][-1]
    rho_final = res['rho_tot'][-1]
    print(f"3) n_rel a t_end: {n_final:.4f}, rho_tot a t_end: {rho_final*1e5:.2f} pcm")
    print("   -> NOTA: canale single-pass senza pozzo termico esterno (no scambiatore).")
    print("      Il moderatore si scalda monotonamente finche' rho_tot->0 e n_rel si stabilizza")
    print("      su un nuovo livello (non torna a 1): costante di tempo = tempo di transito")
    print("      canale, molto piu' lenta del picco prompt. Check corretto: rho_tot deve")
    print("      DECRESCERE in modulo rispetto al picco subito dopo l'inserimento (no divergenza).")
    rho_after_ramp = res['rho_tot'][np.searchsorted(res['t'], 0.15)]
    assert abs(rho_final) <= abs(rho_after_ramp) * 1.05, \
        "rho_tot non decresce nel tempo: possibile instabilita' del feedback!"
    assert np.all(res['n_rel'] < 200), "Potenza diverge oltre limiti fisici plausibili: instabilita'!"

    dTf = res['Tf_avg'][-1] - res['Tf_avg'][0]
    dTc = res['Tc_avg'][-1] - res['Tc_avg'][0]
    print(f"4) Delta T fuel medio: {dTf:.2f} K, Delta T coolant medio: {dTc:.2f} K")
    print("   -> ordine di grandezza atteso per RIA moderato: decine-centinaia K fuel, pochi-decine K coolant")
    assert dTf > 0, "Il combustibile dovrebbe scaldarsi durante un inserimento positivo di reattivita'"

    print("5) Limite rho_ext->0: testato separatamente in coupled_solver.py (res0~1e-12)")
    print("=" * 70)


def make_plots(res, outpath="/mnt/user-data/outputs/ria_transient.png"):
    t = res['t']
    fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axs[0].plot(t, res['n_rel'], color='crimson')
    axs[0].set_ylabel("Potenza relativa n(t)/n0")
    axs[0].set_yscale("log")
    axs[0].grid(True, which='both', alpha=0.3)
    axs[0].set_title("RIA: rampa di reattivita' esterna + feedback Doppler/moderatore")

    axs[1].plot(t, res['rho_ext'] * 1e5, label="rho esterna", color='black', ls='--')
    axs[1].plot(t, res['rho_tot'] * 1e5, label="rho netta (tot)", color='blue')
    axs[1].set_ylabel("Reattivita' [pcm]")
    axs[1].legend()
    axs[1].grid(True, alpha=0.3)

    axs[2].plot(t, res['Tf_avg'] - 273.15, label="T fuel media", color='darkorange')
    axs[2].plot(t, res['Tc_avg'] - 273.15, label="T coolant media", color='teal')
    axs[2].set_ylabel("Temperatura [C]")
    axs[2].set_xlabel("Tempo [s]")
    axs[2].legend()
    axs[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=130)
    print(f"Plot salvato in: {outpath}")


if __name__ == "__main__":
    res = run_case(rho_max_pcm=600.0, t_ramp=0.1, t_end=20.0)
    run_checks(res)
    make_plots(res)
