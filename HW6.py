import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Geometry (cm)
# ------------------------------------------------------------
a_core = 250.0   # half-thickness of core (full core 500 cm)
t_refl = 30.0    # reflector thickness on each side

# Given current at core/reflector interface (from problem)
J_interface = 1.0e10  # n / cm^2-s

# ------------------------------------------------------------

# ------------------------------------------------------------
materials = {
    # These four will be used as CORE materials
    "H2O": {
        "Sigma_s": 3.443,       # <-- PUT Lamarsh Σ_s
        "Sigma_a": 0.0222,     # example
        "A": 18.0              # water molecule effective A (approx)
    },
    "D2O": {
        "Sigma_s": 0.4519,       # <-- PUT Lamarsh Σ_s
        "Sigma_a": 4.42e-5,     # example
        "A": 20.0              # heavy water molecule effective A (approx)
    },
    "Be": {
        "Sigma_s": 0.7589,       # <-- PUT Lamarsh Σ_s
        "Sigma_a": 0.001137,    # example
        "A": 9.0               # Be-9
    },
    "Graphite": {
        "Sigma_s": 0.3811,       # <-- PUT Lamarsh Σ_s
        "Sigma_a": 2.728e-4,     # example
        "A": 12.0              # C-12
    },

    # These two will be used as REFLECTOR materials
    "Fe": {
        "Sigma_s": 0.9251,        # <-- PUT Σ_s from notes / Lamarsh
        "Sigma_a": 0.2164,      # example
        "A": 56.0              # Fe-56
    },
    "Zr": {
        "Sigma_s": 0.2746,        # <-- PUT Σ_s from notes / Lamarsh
        "Sigma_a": 0.007938,      # example
        "A": 91.0              # ~Zr-91
    }
}

# Compute μ_bar, D and L from Σ_s, Σ_a, A
for name, props in materials.items():
    Ss = props["Sigma_s"]
    Sa = props["Sigma_a"]
    A = props["A"]
    mu_bar = 2.0 / (3.0 * A)  # μ̄ = 2/(3A)
    props["mu_bar"] = mu_bar
    props["D"] = 1.0 / (3.0 * Ss * (1.0 - mu_bar))  # D = 1/[3 Σ_s (1-μ̄)]
    props["L"] = np.sqrt(props["D"] / Sa)


# ------------------------------------------------------------
# Root-finding (bisection) for core geometric buckling B1
# ------------------------------------------------------------
def find_B1(D_core, D_ref, L_ref, a, t):
    """
    Solve for geometric buckling B1 using the interface / vacuum matching:

      cos(B a) / (D_core * B * sin(B a))
        = (L_ref / D_ref) * tanh((t + z0) / L_ref),

    where z0 = 2.13 * D_ref is the extrapolation length.
    """
    z0 = 2.13 * D_ref                   # extrapolation distance (cm)
    A_dimless = (t + z0) / L_ref        # dimensionless
    K = (L_ref / D_ref) * np.tanh(A_dimless)

    def f(B):
        # Avoid singularities at B=0, B a = n*pi
        return np.cos(B * a) / (D_core * B * np.sin(B * a)) - K

    # Bracket the first root between ~0 and pi/a
    B_low = 1e-6
    B_high = 0.99 * np.pi / a

    f_low = f(B_low)
    f_high = f(B_high)
    if f_low * f_high > 0:
        raise RuntimeError("Could not bracket root for B1; adjust bounds or parameters.")

    # Bisection method
    for _ in range(100):
        B_mid = 0.5 * (B_low + B_high)
        f_mid = f(B_mid)
        if f_low * f_mid <= 0:
            B_high, f_high = B_mid, f_mid
        else:
            B_low, f_low = B_mid, f_mid
        if abs(B_high - B_low) < 1e-8:
            break

    return 0.5 * (B_low + B_high)


# ------------------------------------------------------------
# Compute flux and current for core + 1-region reflector (RIGHT HALF)
# ------------------------------------------------------------
def core_plus_reflector_profile(D_core, D_ref, L_ref,
                                a, t, J0,
                                n_core=300, n_refl=200):
    """
    Compute flux and current on the RIGHT HALF (x >= 0)
    for a core (0..a) + reflector (a..a+t) using one-group diffusion.

      Core:      d^2 φ1 / dx^2 + B1^2 φ1 = 0
      Reflector: d^2 φ2 / dx^2 - φ2 / L_ref^2 = 0
    """
    # 1) Solve for B1 from matching condition
    B1 = find_B1(D_core, D_ref, L_ref, a, t)

    # 2) Ratio phi(a) / J(a) from reflector side
    z0 = 2.13 * D_ref
    A_dimless = (t + z0) / L_ref
    K = (L_ref / D_ref) * np.tanh(A_dimless)  # phi(a) / J(a)

    # 3) Core amplitude A from: phi_1(a) = A cos(B1 a) = J0 * K
    A_amp = J0 * K / np.cos(B1 * a)

    # 4) Core solution (0 <= x <= a)
    x_core = np.linspace(0.0, a, n_core)
    phi_core = A_amp * np.cos(B1 * x_core)
    J_core = A_amp * D_core * B1 * np.sin(B1 * x_core)

    # 5) Reflector solution (local coord y = x - a)
    x_refl_local = np.linspace(0.0, t, n_refl)
    x_refl = a + x_refl_local

    phi_a = J0 * K          # flux at interface
    C = phi_a               # φ2(0) = C
    D_const = -J0 * L_ref / D_ref  # from current continuity at interface

    phi_refl = C * np.cosh(x_refl_local / L_ref) + \
               D_const * np.sinh(x_refl_local / L_ref)

    J_refl = -D_ref * (C / L_ref * np.sinh(x_refl_local / L_ref) +
                       D_const / L_ref * np.cosh(x_refl_local / L_ref))

    # Concatenate core + reflector (right side)
    x_pos = np.concatenate((x_core, x_refl))
    phi_pos = np.concatenate((phi_core, phi_refl))
    J_pos = np.concatenate((J_core, J_refl))

    return x_pos, phi_pos, J_pos


# ------------------------------------------------------------
# Build full-domain profiles by symmetry for a given core/reflector pair
# ------------------------------------------------------------
def full_domain_profiles(core_name, refl_name):
    """
    Compute flux and current for x in [-(a+t), +(a+t)]
    for a given core material and reflector material.
    """
    D_core = materials[core_name]["D"]
    D_ref = materials[refl_name]["D"]
    L_ref = materials[refl_name]["L"]

    x_pos, phi_pos, J_pos = core_plus_reflector_profile(
        D_core, D_ref, L_ref,
        a_core, t_refl, J_interface
    )

    # Mirror: φ is even, J is odd
    x_full = np.concatenate((-x_pos[::-1], x_pos))
    phi_full = np.concatenate((phi_pos[::-1], phi_pos))
    J_full = np.concatenate((-J_pos[::-1], J_pos))

    return x_full, phi_full, J_full


# ------------------------------------------------------------
# Main: make requested plots
# ------------------------------------------------------------
def main():
    core_materials = ["H2O", "D2O", "Be", "Graphite"]

    # 1) All core materials with Zr as reflector
    refl_Zr = "Zr"

    plt.figure(figsize=(10, 6))
    for core in core_materials:
        x, phi, J = full_domain_profiles(core, refl_Zr)
        plt.plot(x, phi, label=f"Core: {core}")
    plt.xlabel("x (cm)")
    plt.ylabel("Flux ϕ(x) (a.u.)")
    plt.title("Core + Zr reflector: flux vs x (D from Σ_s, μ̄=2/3A, z₀=2.13D)")
    plt.legend()
    plt.grid(True)

    plt.figure(figsize=(10, 6))
    for core in core_materials:
        x, phi, J = full_domain_profiles(core, refl_Zr)
        plt.plot(x, J, label=f"Core: {core}")
    plt.xlabel("x (cm)")
    plt.ylabel("Current J(x) (n/cm$^2$-s, relative)")
    plt.title("Core + Zr reflector: current vs x (D from Σ_s, μ̄=2/3A, z₀=2.13D)")
    plt.legend()
    plt.grid(True)

    # 2) All core materials with Fe as reflector
    refl_Fe = "Fe"

    plt.figure(figsize=(10, 6))
    for core in core_materials:
        x, phi, J = full_domain_profiles(core, refl_Fe)
        plt.plot(x, phi, label=f"Core: {core}")
    plt.xlabel("x (cm)")
    plt.ylabel("Flux ϕ(x) (a.u.)")
    plt.title("Core + Fe reflector: flux vs x (D from Σ_s, μ̄=2/3A, z₀=2.13D)")
    plt.legend()
    plt.grid(True)

    plt.figure(figsize=(10, 6))
    for core in core_materials:
        x, phi, J = full_domain_profiles(core, refl_Fe)
        plt.plot(x, J, label=f"Core: {core}")
    plt.xlabel("x (cm)")
    plt.ylabel("Current J(x) (n/cm$^2$-s, relative)")
    plt.title("Core + Fe reflector: current vs x (D from Σ_s, μ̄=2/3A, z₀=2.13D)")
    plt.legend()
    plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
