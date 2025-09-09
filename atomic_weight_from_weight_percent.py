from typing import Sequence

def atomic_weight_from_weight_percent(
    masses_u: Sequence[float],
    weight_percents: Sequence[float],
) -> float:
    """
    Compute the atomic weight (average atomic mass) from isotopic masses and
    their composition given in weight % (or mass fractions).

    Parameters
    ----------
    masses_u : sequence of float
        Isotopic masses (in unified atomic mass units, u).
    weight_percents : sequence of float
        Composition by *mass*. Can be weight percent (summing ~100)
        or mass fractions (summing ~1). Order must match `masses_u`.

    Returns
    -------
    float
        Atomic weight (u), i.e., the mole-fraction-weighted mean mass.

    Theory
    ------
    Let m_i be isotopic masses, and w_i mass fractions (∑w_i=1). Then
        n_i ∝ w_i / m_i  (moles of isotope i),
        x_i = n_i / ∑n_i = (w_i / m_i) / ∑(w_j / m_j),
        Atomic weight = ∑ x_i m_i = 1 / ∑(w_i / m_i).
    If inputs are in percent (W_i, ∑W_i≈100), the equivalent form is
        Atomic weight = (∑W_i) / ∑(W_i / m_i).
    """
    if len(masses_u) != len(weight_percents):
        raise ValueError("masses_u and weight_percents must have the same length.")
    if any(m <= 0 for m in masses_u):
        raise ValueError("All masses must be positive.")
    if any(w < 0 for w in weight_percents):
        raise ValueError("All weight entries must be nonnegative.")

    total = float(sum(weight_percents))
    if total == 0:
        raise ValueError("Sum of weight percents/fractions is zero.")

    # Auto-detect: if they look like percents (sum ~100), use percent formula; else normalize to fractions
    if abs(total - 100.0) < 1.0:
        # Use percent form directly: AW = (sum W_i) / sum(W_i / m_i)
        denom = sum((w / m) for w, m in zip(weight_percents, masses_u))
        if denom == 0:
            raise ZeroDivisionError("Denominator became zero; check inputs.")
        return total / denom
    else:
        # Treat as arbitrary mass weights → normalize to fractions
        w_frac = [w / total for w in weight_percents]
        denom = sum((w / m) for w, m in zip(w_frac, masses_u))
        if denom == 0:
            raise ZeroDivisionError("Denominator became zero; check inputs.")
        return 1.0 / denom


aw = atomic_weight_from_weight_percent([238.0496,239.0522,240.0538,241.0568,242.0587], [15,35,15,20,15])
print(f"{aw:.6f} u")  

# Same numbers as mass fractions (0.30, 0.70) also work:
aw2 = atomic_weight_from_weight_percent([238.0496,239.0522,240.0538,241.0568,242.0587], [.15, .35, .15, .20, .15])
print(f"{aw2:.6f} u")  # -> 10.777778 u
