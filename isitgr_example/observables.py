#!/usr/bin/env python3
"""
ISiTGR demo script (Python) – only:
  1) (mu, eta) parameterization  [= "mueta"]  (CMB TT and lensing phi-phi)
  2) (mu, Sigma) parameterization [= "muSigma"] (CMB TT and TE)
  3) binning methods example
  
Usage:
  python isitgr_demo_mueta_musigma.py
  python isitgr_demo_mueta_musigma.py --only mueta
  python isitgr_demo_mueta_musigma.py --only musigma
  python isitgr_demo_mueta_musigma.py --save-dir ./figs

Notes:
  - This script assumes you have built ISiTGR (Fortran shared library) and that
    `import isitgr` works *in your terminal environment* (modules, libgfortran, etc).
  - If you are running from inside the ISiTGR repo, it will also try to import
    from the local checkout by adding a sensible path to sys.path.
"""

from __future__ import annotations

import os
import sys
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# SECTION 0 — Import ISiTGR
# -----------------------------------------------------------------------------
def import_isitgr() -> "module":
    """
    Try to import ISiTGR in a way that works both for:
      - an installed package, and
      - a local checkout (e.g., running this script inside the repo).
    """
    # 1) If user provided an override, try that first
    env_path = os.environ.get("ISITGR_PATH", "").strip()
    if env_path:
        sys.path.insert(0, env_path)

    # 2) Otherwise, try local repo root (directory containing this script),
    #    and its parent (in case script lives in docs/ or examples/)
    here = Path(__file__).resolve().parent
    for p in [here, here.parent]:
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))

    try:
        import isitgr  # type: ignore
        print(f"Using CAMB-ISiTGR {isitgr.__version__} at {Path(isitgr.__file__).resolve().parent}")
        return isitgr
    except Exception as exc:
        raise ImportError(
            "Could not import `isitgr`.\n"
            "Make sure you built ISiTGR (e.g. `python setup.py make` in the repo)\n"
            "and that your environment has compatible GCC/libgfortran loaded.\n"
            "If you want to import from a local checkout, set ISITGR_PATH=/path/to/ISiTGR.\n"
        ) from exc


# -----------------------------------------------------------------------------
# SECTION 1 — (mu, eta) parameterization (mueta): reproduce Planck collaboration, Fig. 1 of 1502.01590v2
# -----------------------------------------------------------------------------
def run_mueta(isitgr, save_dir: Path | None = None) -> None:
    # Defining array for different MG parameter values.
    E11 = [1, -1, 0.5, 0]
    E22 = [1, -1, 0.5, 1]

    # Set up different set of parameters for CAMB (including MG)
    pars_GR = isitgr.CAMBparams()
    pars_MG1 = isitgr.CAMBparams()
    pars_MG2 = isitgr.CAMBparams()
    pars_MG3 = isitgr.CAMBparams()
    pars_MG4 = isitgr.CAMBparams()

    # For GR
    pars_GR.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09)
    pars_GR.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_GR.set_for_lmax(2500, lens_potential_accuracy=0)

    # For MG models (mueta)
    pars_MG1.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
                           MG_parameterization="mueta", E11=E11[0], E22=E22[0])
    pars_MG1.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG1.set_for_lmax(2500, lens_potential_accuracy=0)

    pars_MG2.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
                           MG_parameterization="mueta", E11=E11[1], E22=E22[1])
    pars_MG2.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG2.set_for_lmax(2500, lens_potential_accuracy=0)

    pars_MG3.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
                           MG_parameterization="mueta", E11=E11[2], E22=E22[2])
    pars_MG3.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG3.set_for_lmax(2500, lens_potential_accuracy=0)

    pars_MG4.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
                           MG_parameterization="mueta", E11=E11[3], E22=E22[3])
    pars_MG4.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG4.set_for_lmax(2500, lens_potential_accuracy=0)

    # Calculate results for different models.
    results_GR = isitgr.get_results(pars_GR)
    results_MG1 = isitgr.get_results(pars_MG1)
    results_MG2 = isitgr.get_results(pars_MG2)
    results_MG3 = isitgr.get_results(pars_MG3)
    results_MG4 = isitgr.get_results(pars_MG4)

    # Get dictionary of CAMB power spectra for different models, including GR and MG models.
    powers_GR = results_GR.get_cmb_power_spectra(pars_GR, CMB_unit="muK")
    powers_MG1 = results_MG1.get_cmb_power_spectra(pars_MG1, CMB_unit="muK")
    powers_MG2 = results_MG2.get_cmb_power_spectra(pars_MG2, CMB_unit="muK")
    powers_MG3 = results_MG3.get_cmb_power_spectra(pars_MG3, CMB_unit="muK")
    powers_MG4 = results_MG4.get_cmb_power_spectra(pars_MG4, CMB_unit="muK")

    # Total lensed CMB power spectra (TT is column 0)
    totCL_GR = powers_GR["total"]
    totCL_MG1 = powers_MG1["total"]
    totCL_MG2 = powers_MG2["total"]
    totCL_MG3 = powers_MG3["total"]
    totCL_MG4 = powers_MG4["total"]

    # Lensing potential cl's (phi-phi is column 0)
    pars_GR.set_dark_energy(w=-1, wa=0, dark_energy_model="fluid")
    cl_GR = results_GR.get_lens_potential_cls(lmax=2550)

    pars_MG1.set_dark_energy(w=-1, wa=0, dark_energy_model="fluid")
    cl_MG1 = results_MG1.get_lens_potential_cls(lmax=2550)

    pars_MG2.set_dark_energy(w=-1, wa=0, dark_energy_model="fluid")
    cl_MG2 = results_MG2.get_lens_potential_cls(lmax=2550)

    pars_MG3.set_dark_energy(w=-1, wa=0, dark_energy_model="fluid")
    cl_MG3 = results_MG3.get_lens_potential_cls(lmax=2550)

    pars_MG4.set_dark_energy(w=-1, wa=0, dark_energy_model="fluid")
    cl_MG4 = results_MG4.get_lens_potential_cls(lmax=2550)

    Tspectra = [totCL_GR[:, 0], totCL_MG1[:, 0], totCL_MG2[:, 0], totCL_MG3[:, 0], totCL_MG4[:, 0]]
    Lspectra = [cl_GR[:, 0], cl_MG1[:, 0], cl_MG2[:, 0], cl_MG3[:, 0], cl_MG4[:, 0]]
    ell = np.arange(2551)

    # Plot
    import matplotlib.ticker as mtick
    from matplotlib import gridspec

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)

    ax0.semilogx(ell, Tspectra[0], label=r"$\Lambda$CDM")
    ax0.semilogx(ell, Tspectra[1], label=r"$E_{11}=1,\ E_{22}=1$")
    ax0.semilogx(ell, Tspectra[2], linestyle="dashdot", label=r"$E_{11}=-1,\ E_{22}=-1$")
    ax0.semilogx(ell, Tspectra[3], linestyle="dotted", label=r"$E_{11}=0.5,\ E_{22}=0.5$")
    ax0.semilogx(ell, Tspectra[4], linestyle="dashed", label=r"$E_{11}=0,\ E_{22}=1$")

    ax1.semilogx(ell, Lspectra[0], label=r"$\Lambda$CDM")
    ax1.semilogx(ell, Lspectra[1], label=r"$E_{11}=1,\ E_{22}=1$")
    ax1.semilogx(ell, Lspectra[2], linestyle="dashdot", label=r"$E_{11}=-1,\ E_{22}=-1$")
    ax1.semilogx(ell, Lspectra[3], linestyle="dotted", label=r"$E_{11}=0.5,\ E_{22}=0.5$")
    ax1.semilogx(ell, Lspectra[4], linestyle="dashed", label=r"$E_{11}=0,\ E_{22}=1$")

    fig.subplots_adjust(hspace=0.03)

    ax0.set_title("ISiTGR demo: (mu, eta) parameterization")
    ax0.set_ylabel(r"$\ell(\ell+1)C_{\ell}^{TT}/2\pi$ [$\mu K^2$]", fontsize=12)
    ax1.set_ylabel(r"$[\ell(\ell+1)]^2C_{\ell}^{\phi\phi}/2\pi$", fontsize=12)
    ax1.set_xlabel(r"$\ell$", fontsize=12)

    ax0.set_xlim(2, 2100)
    ax1.set_xlim(2, 2100)
    ax0.set_ylim(0, 8000)
    ax1.set_ylim(0, 4e-7)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.e"))
    ax1.legend(loc="upper right", fontsize=10)

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / "isitgr_mueta_TT_phiphi.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[mueta] saved: {out}")

    plt.show()


# -----------------------------------------------------------------------------
# SECTION 2 — (mu, Sigma) parameterization (muSigma): reproduce CFHTLens collaboration Fig. 4 of 1212.3339v2
# -----------------------------------------------------------------------------
def run_musigma(isitgr, save_dir: Path | None = None) -> None:
    # Defining array for different MG parameter values.
    mu0 = [0.5, 0]
    Sigma0 = [0, -0.5]

    pars_GR = isitgr.CAMBparams()
    pars_MG1 = isitgr.CAMBparams()
    pars_MG2 = isitgr.CAMBparams()

    # GR
    pars_GR.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09)
    pars_GR.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_GR.set_for_lmax(500, lens_potential_accuracy=0)

    # MG1 (muSigma)
    pars_MG1.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
                           MG_parameterization="muSigma", mu0=mu0[0], Sigma0=Sigma0[0])
    pars_MG1.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG1.set_for_lmax(500, lens_potential_accuracy=0)

    # MG2 (muSigma)
    pars_MG2.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
                           MG_parameterization="muSigma", mu0=mu0[1], Sigma0=Sigma0[1])
    pars_MG2.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG2.set_for_lmax(500, lens_potential_accuracy=0)

    # Results
    results_GR = isitgr.get_results(pars_GR)
    results_MG1 = isitgr.get_results(pars_MG1)
    results_MG2 = isitgr.get_results(pars_MG2)

    powers_GR = results_GR.get_cmb_power_spectra(pars_GR, CMB_unit="muK")
    powers_MG1 = results_MG1.get_cmb_power_spectra(pars_MG1, CMB_unit="muK")
    powers_MG2 = results_MG2.get_cmb_power_spectra(pars_MG2, CMB_unit="muK")

    totCL_GR = powers_GR["total"]
    totCL_MG1 = powers_MG1["total"]
    totCL_MG2 = powers_MG2["total"]

    Tspectra = [totCL_GR[:, 0], totCL_MG1[:, 0], totCL_MG2[:, 0]]  # TT
    TEspectra = [totCL_GR[:, 3], totCL_MG1[:, 3], totCL_MG2[:, 3]]  # TE

    ell = np.arange(551)
    ell_safe = ell.copy()
    ell_safe[0] = 1
    TEspectra_divided_by_ell = np.array(TEspectra) / ell_safe

    from matplotlib import gridspec
    fig = plt.figure(figsize=(10, 15))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)

    ax0.loglog(ell, Tspectra[0], label=r"$\Lambda$CDM")
    ax0.loglog(ell, Tspectra[1], linestyle="dashed", label=r"$\mu_0=0.5,\ \Sigma_0=0$")
    ax0.loglog(ell, Tspectra[2], linestyle="dashdot", label=r"$\mu_0=0,\ \Sigma_0=-0.5$")

    ax1.semilogx(ell, TEspectra_divided_by_ell[0], label=r"$\Lambda$CDM")
    ax1.semilogx(ell, TEspectra_divided_by_ell[1], linestyle="dashed", label=r"$\mu_0=0.5,\ \Sigma_0=0$")
    ax1.semilogx(ell, TEspectra_divided_by_ell[2], linestyle="dashdot", label=r"$\mu_0=0,\ \Sigma_0=-0.5$")

    fig.subplots_adjust(hspace=0.03)

    ax0.set_title("ISiTGR demo: (mu, Sigma) parameterization")
    ax0.set_ylabel(r"$\ell(\ell+1)C_{\ell}^{TT}/2\pi$ [$\mu K^2$]", fontsize=12)
    ax1.set_ylabel(r"$(\ell+1)C_{\ell}^{TE}/2\pi$ [$\mu K^2$]", fontsize=12)
    ax1.set_xlabel(r"$\ell$", fontsize=12)

    ax0.set_xlim(2, 500)
    ax1.set_xlim(2, 500)
    ax0.set_ylim(100, 100000)
    ax1.set_ylim(-1.5, 3)
    ax0.legend(loc="upper left", fontsize=10)

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / "isitgr_musigma_TT_TE.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[musigma] saved: {out}")

    plt.show()


# -----------------------------------------------------------------------------
# SECTION 3 — mu(z) from redshift-binned muSigma (example "binning method")
# -----------------------------------------------------------------------------
def run_musigma_redshift_binning_mu(isitgr, save_dir: Path | None = None) -> None:
    # Defining z,k arrays
    z = np.linspace(0.0, 2.5, 200)
    k = np.array([0.01])  # MG is scale independent here so this doesn't matter

    # -------------------------------------------------------------------------
    # Build 3 explicit parameter sets
    # -------------------------------------------------------------------------
    pars_MG1 = isitgr.CAMBparams()
    pars_MG1.set_cosmology(
        H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
        MG_parameterization="muSigma",
        z_TGR=2.0, z_tw=0.01,
        redshift_bins=True,
        mu1=1.4, mu2=1.3, mu3=1.2, mu4=1.1,
        Sigma1=1.0, Sigma2=1.0, Sigma3=1.0, Sigma4=1.0,
    )
    pars_MG1.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG1.set_for_lmax(2500, lens_potential_accuracy=0)
    results_MG1 = isitgr.get_results(pars_MG1)
    mu_MG1 = results_MG1.mu_MG(pars_MG1, z, k)[:, 0]

    pars_MG2 = isitgr.CAMBparams()
    pars_MG2.set_cosmology(
        H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
        MG_parameterization="muSigma",
        z_TGR=2.0, z_tw=0.05,
        redshift_bins=True,
        mu1=1.4, mu2=1.3, mu3=1.2, mu4=1.1,
        Sigma1=1.0, Sigma2=1.0, Sigma3=1.0, Sigma4=1.0,
    )
    pars_MG2.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG2.set_for_lmax(2500, lens_potential_accuracy=0)
    results_MG2 = isitgr.get_results(pars_MG2)
    mu_MG2 = results_MG2.mu_MG(pars_MG2, z, k)[:, 0]

    pars_MG3 = isitgr.CAMBparams()
    pars_MG3.set_cosmology(
        H0=70, ombh2=0.0226, omch2=0.112, mnu=0, omk=0, tau=0.09,
        MG_parameterization="muSigma",
        z_TGR=2.0, z_tw=0.15,
        redshift_bins=True,
        mu1=1.4, mu2=1.3, mu3=1.2, mu4=1.1,
        Sigma1=1.0, Sigma2=1.0, Sigma3=1.0, Sigma4=1.0,
    )
    pars_MG3.InitPower.set_params(As=2.1e-9, ns=0.96, r=0)
    pars_MG3.set_for_lmax(2500, lens_potential_accuracy=0)
    results_MG3 = isitgr.get_results(pars_MG3)
    mu_MG3 = results_MG3.mu_MG(pars_MG3, z, k)[:, 0]

    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(8, 5))
    plt.xlabel(r"$z$", fontsize=16)
    plt.ylabel(r"$\mu(z)$", fontsize=16)
    plt.hlines(1.0, 0, 2.5, color="gray", linestyles="--")

    plt.plot(z, mu_MG1, label="z_tw = 0.01")
    plt.plot(z, mu_MG2, linestyle="--", label="z_tw = 0.05")
    plt.plot(z, mu_MG3, linestyle="-.", label="z_tw = 0.15")

    plt.legend(loc="upper right", fontsize=12)
    plt.title("obtained using ISiTGR-python wrapper")

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / "isitgr_musigma_redshift_binning_mu.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[musigma-binning] saved: {out}")

    plt.show()

# -----------------------------------------------------------------------------
# SECTION 4 — CLI
# -----------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--only",
        choices=["mueta", "musigma", "musigma_binning", "both"],
        default="both",
        help="Which demo(s) to run.",
    )
    p.add_argument(
        "--save-dir",
        type=Path,
        default="./",
        help="If set, save figures to this directory (PNG).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    isitgr = import_isitgr()

    if args.only in ("mueta", "both"):
        run_mueta(isitgr, save_dir=args.save_dir)
    if args.only in ("musigma", "both"):
        run_musigma(isitgr, save_dir=args.save_dir)
    if args.only in ("musigma_binning", "both"):
        run_musigma_redshift_binning_mu(isitgr, save_dir=args.save_dir)


if __name__ == "__main__":
    main()
