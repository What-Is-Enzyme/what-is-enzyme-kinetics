#!/usr/bin/env python3
"""
enzyme_kinetics_fitter.py

Fit Michaelis–Menten enzyme kinetics to experimental data.

Takes a CSV file with substrate concentration [S] and initial velocity v,
fits Vmax and Km using non-linear least squares, and generates plots.

Usage examples:
  python enzyme_kinetics_fitter.py --input data.csv

  python enzyme_kinetics_fitter.py \
      --input data.csv \
      --substrate-col Substrate_mM \
      --velocity-col Rate_uM_s \
      --output-plot kinetics_fit.png
"""

import argparse
import json
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def michaelis_menten(S, Vmax, Km):
    """Michaelis–Menten equation: v = (Vmax * [S]) / (Km + [S])"""
    return (Vmax * S) / (Km + S)


@dataclass
class FitResult:
    Vmax: float
    Km: float
    Vmax_std: float | None
    Km_std: float | None
    r_squared: float | None


def fit_michaelis_menten(S: np.ndarray, v: np.ndarray) -> FitResult:
    """Fit Michaelis–Menten parameters to data."""
    # Clean data: remove NaN and negative values
    mask = (~np.isnan(S)) & (~np.isnan(v)) & (S > 0) & (v > 0)
    S_clean = S[mask]
    v_clean = v[mask]

    if len(S_clean) < 3:
        raise ValueError("Not enough valid data points after cleaning (need ≥ 3).")

    # Initial guesses
    Vmax_guess = np.max(v_clean)
    Km_guess = np.median(S_clean)

    try:
        popt, pcov = curve_fit(
            michaelis_menten,
            S_clean,
            v_clean,
            p0=[Vmax_guess, Km_guess],
            maxfev=10000,
        )
        Vmax_fit, Km_fit = popt
        # Standard deviations from covariance matrix
        perr = np.sqrt(np.diag(pcov))
        Vmax_std, Km_std = perr

        # R^2 calculation
        v_pred = michaelis_menten(S_clean, Vmax_fit, Km_fit)
        ss_res = np.sum((v_clean - v_pred) ** 2)
        ss_tot = np.sum((v_clean - np.mean(v_clean)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else None

        return FitResult(
            Vmax=Vmax_fit,
            Km=Km_fit,
            Vmax_std=Vmax_std,
            Km_std=Km_std,
            r_squared=r_squared,
        )
    except Exception as e:
        raise RuntimeError(f"Curve fitting failed: {e}") from e


def plot_fit(S: np.ndarray,
             v: np.ndarray,
             fit: FitResult,
             output_path: str | None = None,
             title: str | None = None):
    """Generate plot with data + fitted curve + residuals."""
    # Clean data same way as fitting
    mask = (~np.isnan(S)) & (~np.isnan(v)) & (S > 0) & (v > 0)
    S_clean = S[mask]
    v_clean = v[mask]

    S_range = np.linspace(0, S_clean.max() * 1.1, 200)
    v_fit_curve = michaelis_menten(S_range, fit.Vmax, fit.Km)
    v_pred = michaelis_menten(S_clean, fit.Vmax, fit.Km)
    residuals = v_clean - v_pred

    fig = plt.figure(figsize=(8, 8))

    # Top: data + fit
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.scatter(S_clean, v_clean, label="Data")
    ax1.plot(S_range, v_fit_curve, label="Fit (Michaelis–Menten)")
    ax1.set_xlabel("[S]")
    ax1.set_ylabel("v")
    ax1.set_title(title or "Enzyme Kinetics Fit (Michaelis–Menten)")
    ax1.legend()
    ax1.grid(True, linestyle="--", alpha=0.4)

    # Bottom: residuals
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.axhline(0, linestyle="--")
    ax2.scatter(S_clean, residuals)
    ax2.set_xlabel("[S]")
    ax2.set_ylabel("Residuals (v_obs - v_fit)")
    ax2.set_title("Residuals")
    ax2.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"[INFO] Plot saved to: {output_path}")
    else:
        plt.show()

    plt.close(fig)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Fit Michaelis–Menten kinetics (Vmax & Km) from CSV data."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to CSV file containing substrate and velocity data.",
    )
    parser.add_argument(
        "--substrate-col",
        default=None,
        help="Name of the substrate concentration column (default: first numeric column).",
    )
    parser.add_argument(
        "--velocity-col",
        default=None,
        help="Name of the velocity column (default: second numeric column).",
    )
    parser.add_argument(
        "--output-plot",
        "-o",
        default="kinetics_fit.png",
        help="Output path for the plot PNG file (default: kinetics_fit.png).",
    )
    parser.add_argument(
        "--output-json",
        default=None,
        help="Optional: output fit parameters as JSON to this file.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="If set, do not generate a plot.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title for the plot.",
    )
    return parser.parse_args(argv)


def guess_numeric_columns(df: pd.DataFrame):
    """Return names of first two numeric columns, or raise error."""
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) < 2:
        raise ValueError("CSV must contain at least two numeric columns.")
    return numeric_cols[0], numeric_cols[1]


def main(argv=None):
    args = parse_args(argv)

    # Load CSV
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"[ERROR] Failed to read CSV: {e}")
        sys.exit(1)

    # Detect columns
    if args.substrate_col is None or args.velocity_col is None:
        s_col, v_col = guess_numeric_columns(df)
        print(f"[INFO] Using numeric columns: [S] = '{s_col}', v = '{v_col}'")
    else:
        s_col, v_col = args.substrate_col, args.velocity_col

    if s_col not in df.columns or v_col not in df.columns:
        print(f"[ERROR] Columns not found in CSV. Available columns: {list(df.columns)}")
        sys.exit(1)

    S = df[s_col].to_numpy(dtype=float)
    v = df[v_col].to_numpy(dtype=float)

    # Fit
    try:
        fit = fit_michaelis_menten(S, v)
    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    # Print results
    print("\n=== Michaelis–Menten Fit Results ===")
    print(f"Vmax      = {fit.Vmax:.4g}")
    if fit.Vmax_std is not None:
        print(f"Vmax_std  = {fit.Vmax_std:.4g}")
    print(f"Km        = {fit.Km:.4g}")
    if fit.Km_std is not None:
        print(f"Km_std    = {fit.Km_std:.4g}")
    if fit.r_squared is not None:
        print(f"R²        = {fit.r_squared:.4f}")
    print("====================================\n")

    # Save JSON if requested
    if args.output_json:
        payload = {
            "Vmax": fit.Vmax,
            "Km": fit.Km,
            "Vmax_std": fit.Vmax_std,
            "Km_std": fit.Km_std,
            "r_squared": fit.r_squared,
            "input_csv": args.input,
            "substrate_col": s_col,
            "velocity_col": v_col,
        }
        with open(args.output_json, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"[INFO] Fit parameters saved to JSON: {args.output_json}")

    # Plot
    if not args.no_plot:
        plot_fit(
            S,
            v,
            fit,
            output_path=args.output_plot,
            title=args.title,
        )


if __name__ == "__main__":
    main()
