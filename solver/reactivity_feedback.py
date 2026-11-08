"""
reactivity_feedback.py
Modello di reattivita': esterna (rampa/step, input RIA) + feedback Doppler
(combustibile) + feedback moderatore (densita'/temperatura coolant).

rho_fb = alpha_D * (Tf_avg - Tf_avg_0) + alpha_m * (Tc_avg - Tc_avg_0)

Convenzione: alpha_D, alpha_m in [1/K] (reattivita' per Kelvin), entrambi
tipicamente negativi in un PWR a combustibile basso-arricchito (feedback
stabilizzante).
"""

from dataclasses import dataclass


@dataclass
class FeedbackParams:
    alpha_D: float = -3.0e-5     # [1/K] coefficiente Doppler tipico PWR: -1..-5 pcm/K
    alpha_m: float = -2.0e-4     # [1/K] coefficiente moderatore (puo' essere piu' negativo a fine vita)
    Tf_ref: float = 0.0          # [K] settati a runtime = valori stazionari pre-transitorio
    Tc_ref: float = 0.0


class ExternalReactivity:
    """Inserimento di reattivita' esterna: step, rampa, o combinazione (RIA tipico)."""

    def __init__(self, kind: str = "step", rho_max: float = 0.0,
                 t_start: float = 0.0, t_ramp: float = 0.0):
        self.kind = kind
        self.rho_max = rho_max
        self.t_start = t_start
        self.t_ramp = t_ramp

    def __call__(self, t: float) -> float:
        if t < self.t_start:
            return 0.0
        if self.kind == "step":
            return self.rho_max
        elif self.kind == "ramp":
            if self.t_ramp <= 0:
                return self.rho_max
            frac = min((t - self.t_start) / self.t_ramp, 1.0)
            return self.rho_max * frac
        else:
            raise ValueError(f"Tipo reattivita' esterna non riconosciuto: {self.kind}")


class ReactivityModel:
    def __init__(self, fb: FeedbackParams, ext: ExternalReactivity):
        self.fb = fb
        self.ext = ext

    def total_rho(self, t: float, Tf_avg: float, Tc_avg: float) -> float:
        rho_ext = self.ext(t)
        rho_doppler = self.fb.alpha_D * (Tf_avg - self.fb.Tf_ref)
        rho_moder = self.fb.alpha_m * (Tc_avg - self.fb.Tc_ref)
        return rho_ext + rho_doppler + rho_moder

    # --- check di modulo ---
    def _check_zero_feedback_at_ref(self):
        rho = self.total_rho(t=-1.0, Tf_avg=self.fb.Tf_ref, Tc_avg=self.fb.Tc_ref)
        assert abs(rho) < 1e-15, "Feedback non nullo a T=T_ref con t<t_start!"
        print("Check OK: rho_feedback(T=T_ref) = 0, rho_ext(t<t_start)=0 => rho_tot=0")


if __name__ == "__main__":
    fb = FeedbackParams(Tf_ref=900.0, Tc_ref=580.0)
    ext = ExternalReactivity(kind="ramp", rho_max=0.005, t_start=0.0, t_ramp=0.1)
    model = ReactivityModel(fb, ext)
    model._check_zero_feedback_at_ref()
    print("rho_ext(t=0.05) =", ext(0.05), " (atteso 0.0025, meta' rampa)")
    print("rho_ext(t=0.2)  =", ext(0.2), " (atteso 0.005, saturato)")
