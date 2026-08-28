#!/usr/bin/env python3
"""Run the PARMA flux executable over an altitude/threshold grid."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Iterable, Sequence

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_EXECUTABLE = PROJECT_ROOT / "electron_fluxes"
RESULT_PREFIX = "RESULT_CSV"


def resolve_executable(path: Path) -> Path:
    """Return the executable path, allowing the Windows '.exe' suffix."""
    if path.is_file():
        return path
    candidate = path.parent / (path.name + ".exe")
    if candidate.is_file():
        return candidate
    return path

PARTICLE_NAMES = {
    0: "neutron",
    29: "muon+",
    30: "muon-",
    31: "electron",
    32: "positron",
    33: "photon",
}


@dataclass(frozen=True)
class FluxResult:
    particle_id: int
    min_energy: float
    max_energy: float
    altitude_km: float
    latitude_deg: float
    longitude_deg: float
    year: int
    month: int
    day: int
    solar_w: float
    solar_status: int
    rigidity_gv: float
    atmospheric_depth_g_cm2: float
    geometry: float
    atmosphere_model: int
    continuum_flux_cm2_s: float
    line_511_flux_cm2_s: float
    total_flux_cm2_s: float


def particle_name(particle_id: int) -> str:
    if 1 <= particle_id <= 28:
        return "H-Ni ion"
    return PARTICLE_NAMES.get(particle_id, "unknown")


def parse_result(stdout: str) -> FluxResult:
    """Parse the stable RESULT_CSV record emitted by the Fortran program."""
    result_lines = [line for line in stdout.splitlines() if line.startswith(f"{RESULT_PREFIX},")]
    if len(result_lines) != 1:
        raise ValueError(
            f"Expected exactly one {RESULT_PREFIX} line, found {len(result_lines)}.\n"
            f"Program output:\n{stdout}"
        )

    row = [value.strip() for value in next(csv.reader([result_lines[0]]))]
    if len(row) != 19 or row[0] != RESULT_PREFIX:
        raise ValueError(f"Malformed {RESULT_PREFIX} record: {result_lines[0]}")

    try:
        result = FluxResult(
            particle_id=int(row[1]),
            min_energy=float(row[2]),
            max_energy=float(row[3]),
            altitude_km=float(row[4]),
            latitude_deg=float(row[5]),
            longitude_deg=float(row[6]),
            year=int(row[7]),
            month=int(row[8]),
            day=int(row[9]),
            solar_w=float(row[10]),
            solar_status=int(row[11]),
            rigidity_gv=float(row[12]),
            atmospheric_depth_g_cm2=float(row[13]),
            geometry=float(row[14]),
            atmosphere_model=int(row[15]),
            continuum_flux_cm2_s=float(row[16]),
            line_511_flux_cm2_s=float(row[17]),
            total_flux_cm2_s=float(row[18]),
        )
    except ValueError as exc:
        raise ValueError(f"Non-numeric value in {RESULT_PREFIX} record: {result_lines[0]}") from exc

    numeric_values = (
        result.min_energy,
        result.max_energy,
        result.altitude_km,
        result.latitude_deg,
        result.longitude_deg,
        result.solar_w,
        result.rigidity_gv,
        result.atmospheric_depth_g_cm2,
        result.geometry,
        result.continuum_flux_cm2_s,
        result.line_511_flux_cm2_s,
        result.total_flux_cm2_s,
    )
    if not all(math.isfinite(value) for value in numeric_values):
        raise ValueError(f"Non-finite value in {RESULT_PREFIX} record: {result_lines[0]}")
    if not 0 <= result.particle_id <= 33:
        raise ValueError(f"Invalid particle ID in {RESULT_PREFIX} record: {result_lines[0]}")
    if result.atmosphere_model not in (0, 1):
        raise ValueError(f"Invalid atmosphere model in {RESULT_PREFIX} record: {result_lines[0]}")
    if result.solar_status not in (0, 1, 2, 3):
        raise ValueError(f"Invalid solar status in {RESULT_PREFIX} record: {result_lines[0]}")
    if result.min_energy <= 0.0 or result.max_energy <= result.min_energy:
        raise ValueError(f"Invalid energy range in {RESULT_PREFIX} record: {result_lines[0]}")
    nonnegative_values = (
        result.solar_w,
        result.rigidity_gv,
        result.atmospheric_depth_g_cm2,
        result.continuum_flux_cm2_s,
        result.line_511_flux_cm2_s,
        result.total_flux_cm2_s,
    )
    if any(value < 0.0 for value in nonnegative_values):
        raise ValueError(f"Negative physical value in {RESULT_PREFIX} record: {result_lines[0]}")
    expected_total = result.continuum_flux_cm2_s + result.line_511_flux_cm2_s
    if not math.isclose(result.total_flux_cm2_s, expected_total, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError(f"Inconsistent flux total in {RESULT_PREFIX} record: {result_lines[0]}")
    return result


def run_executable_with_parameters(
    particle_id: int,
    energy_threshold: float,
    altitude: float,
    latitude: float,
    longitude: float,
    *,
    year: int = 2019,
    month: int = 5,
    day: int = 27,
    geometry: float = 0.15,
    atmosphere: str = "standard",
    solar_w: float | None = None,
    executable: Path = DEFAULT_EXECUTABLE,
    timeout: float = 60.0,
) -> FluxResult:
    """Execute one model calculation and return its parsed result."""
    executable = resolve_executable(executable.expanduser().resolve())
    if not executable.is_file():
        raise FileNotFoundError(f"Executable not found: {executable}. Run 'make' first.")
    if not math.isfinite(timeout) or timeout <= 0.0:
        raise ValueError("Timeout must be a finite number greater than zero.")

    command = [
        str(executable),
        str(particle_id),
        format(energy_threshold, ".17g"),
        format(altitude, ".17g"),
        format(latitude, ".17g"),
        format(longitude, ".17g"),
        str(year),
        str(month),
        str(day),
        format(geometry, ".17g"),
        atmosphere,
    ]
    if solar_w is not None:
        command.append(format(solar_w, ".17g"))

    try:
        process = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"Flux calculation exceeded the {timeout:g}-second timeout.") from exc
    except OSError as exc:
        raise RuntimeError(f"Could not execute {executable}: {exc}") from exc

    if process.returncode != 0:
        details = "\n".join(part for part in (process.stderr.strip(), process.stdout.strip()) if part)
        raise RuntimeError(
            f"Flux executable failed with exit status {process.returncode}."
            + (f"\n{details}" if details else "")
        )

    return parse_result(process.stdout)


def inclusive_decimal_range(start: float, stop: float, step: float) -> list[float]:
    """Return an inclusive decimal grid without binary floating-point drift."""
    start_d = Decimal(str(start))
    stop_d = Decimal(str(stop))
    step_d = Decimal(str(step))
    if step_d <= 0:
        raise ValueError("Altitude step must be greater than zero.")
    if stop_d < start_d:
        raise ValueError("Maximum altitude must be at least the minimum altitude.")

    values: list[float] = []
    current = start_d
    while current <= stop_d:
        values.append(float(current))
        current += step_d
    return values


def output_path(path: Path) -> Path:
    expanded = path.expanduser()
    return expanded if expanded.is_absolute() else PROJECT_ROOT / expanded


def write_results(path: Path, results: Sequence[FluxResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "Particle ID",
                "Particle",
                "Energy Threshold",
                "Energy Unit",
                "Altitude (km)",
                "Latitude (deg)",
                "Longitude (deg)",
                "Year",
                "Month",
                "Day",
                "Solar W Index",
                "Solar Status",
                "Cut-off Rigidity (GV)",
                "Atmospheric Depth (g/cm2)",
                "Geometry",
                "Atmosphere Model",
                "Continuum Flux (cm^-2 s^-1)",
                "511 keV Line Flux (cm^-2 s^-1)",
                "Total Flux (cm^-2 s^-1)",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.particle_id,
                    particle_name(result.particle_id),
                    format(result.min_energy, ".16g"),
                    "MeV/nucleon" if 1 <= result.particle_id <= 28 else "MeV",
                    format(result.altitude_km, ".16g"),
                    format(result.latitude_deg, ".16g"),
                    format(result.longitude_deg, ".16g"),
                    result.year,
                    result.month,
                    result.day,
                    format(result.solar_w, ".16g"),
                    result.solar_status,
                    format(result.rigidity_gv, ".16g"),
                    format(result.atmospheric_depth_g_cm2, ".16g"),
                    format(result.geometry, ".16g"),
                    "standard" if result.atmosphere_model == 0 else "msis",
                    format(result.continuum_flux_cm2_s, ".16g"),
                    format(result.line_511_flux_cm2_s, ".16g"),
                    format(result.total_flux_cm2_s, ".16g"),
                ]
            )


def plot_results(path: Path, results: Sequence[FluxResult], show: bool) -> None:
    try:
        import matplotlib

        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib; install requirements.txt.") from exc

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots()
    thresholds = sorted({result.min_energy for result in results})
    for threshold in thresholds:
        subset = sorted(
            (result for result in results if result.min_energy == threshold),
            key=lambda result: result.altitude_km,
        )
        axis.plot(
            [result.altitude_km for result in subset],
            [result.total_flux_cm2_s for result in subset],
            label=f"Minimum energy = {threshold:g} "
            + ("MeV/nucleon" if 1 <= subset[0].particle_id <= 28 else "MeV"),
        )

    first = results[0]
    axis.set_xlabel("Altitude (km)")
    axis.set_ylabel("Angular- and energy-integrated flux (cm$^{-2}$ s$^{-1}$)")
    axis.set_title(
        f"{particle_name(first.particle_id).capitalize()}, "
        f"latitude {first.latitude_deg:g}°, longitude {first.longitude_deg:g}°"
    )
    axis.grid(True)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path)
    if show:
        plt.show()
    plt.close(figure)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particle-id", type=int, default=31)
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.3, 1.0])
    parser.add_argument("--altitude-min", type=float, default=10.0)
    parser.add_argument("--altitude-max", type=float, default=20.0)
    parser.add_argument("--altitude-step", type=float, default=0.2)
    parser.add_argument("--latitude", type=float, default=28.7)
    parser.add_argument("--longitude", type=float, default=-80.8)
    parser.add_argument("--year", type=int, default=2019)
    parser.add_argument("--month", type=int, default=5)
    parser.add_argument("--day", type=int, default=27)
    parser.add_argument("--geometry", type=float, default=0.15)
    parser.add_argument("--atmosphere", choices=("standard", "msis"), default="standard")
    parser.add_argument("--solar-w", type=float, default=None, help="Override the date-based solar W index.")
    parser.add_argument("--output", type=Path, default=Path("results.csv"))
    parser.add_argument("--plot", type=Path, default=Path("flux_vs_altitude_electron.png"))
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--show", action="store_true", help="Display the plot interactively after saving it.")
    parser.add_argument("--executable", type=Path, default=DEFAULT_EXECUTABLE)
    parser.add_argument("--timeout", type=float, default=60.0)
    return parser


def calculate_grid(args: argparse.Namespace) -> list[FluxResult]:
    altitudes = inclusive_decimal_range(args.altitude_min, args.altitude_max, args.altitude_step)
    results: list[FluxResult] = []
    for threshold in args.thresholds:
        for altitude in altitudes:
            result = run_executable_with_parameters(
                args.particle_id,
                threshold,
                altitude,
                args.latitude,
                args.longitude,
                year=args.year,
                month=args.month,
                day=args.day,
                geometry=args.geometry,
                atmosphere=args.atmosphere,
                solar_w=args.solar_w,
                executable=args.executable,
                timeout=args.timeout,
            )
            results.append(result)
            print(
                f"threshold={threshold:g}, altitude={altitude:g} km, "
                f"flux={result.total_flux_cm2_s:.10e} cm^-2 s^-1"
            )
    return results


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        results = calculate_grid(args)
        csv_path = output_path(args.output)
        write_results(csv_path, results)
        print(f"Wrote {len(results)} rows to {csv_path}")
        if not args.no_plot:
            plot_path = output_path(args.plot)
            plot_results(plot_path, results, args.show)
            print(f"Wrote plot to {plot_path}")
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
