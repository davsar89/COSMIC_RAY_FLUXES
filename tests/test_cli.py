from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from run_on_grid import main, parse_result, resolve_executable, run_executable_with_parameters  # noqa: E402

EXECUTABLE = resolve_executable(ROOT / "electron_fluxes")


class FluxCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not EXECUTABLE.is_file():
            raise RuntimeError("electron_fluxes is missing; run make before the tests")

    def run_model(self, *arguments: object) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(EXECUTABLE), *(str(argument) for argument in arguments)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )

    def test_default_calculation_is_finite_and_uses_valid_daily_data(self) -> None:
        process = self.run_model(31, 0.3, 15, 20, -80)
        self.assertEqual(process.returncode, 0, process.stderr)
        result = parse_result(process.stdout)
        self.assertEqual((result.year, result.month, result.day), (2019, 5, 27))
        self.assertEqual(result.solar_status, 1)
        self.assertTrue(math.isfinite(result.total_flux_cm2_s))
        self.assertGreater(result.total_flux_cm2_s, 0.0)

    def test_all_documented_particle_ids_run(self) -> None:
        for particle_id in range(34):
            with self.subTest(particle_id=particle_id):
                process = self.run_model(particle_id, 0.3, 15, 20, -80)
                self.assertEqual(process.returncode, 0, process.stderr)
                result = parse_result(process.stdout)
                self.assertEqual(result.particle_id, particle_id)
                self.assertGreaterEqual(result.total_flux_cm2_s, 0.0)

    def test_invalid_inputs_fail_with_nonzero_status(self) -> None:
        invalid_argument_sets = (
            (-1, 0.3, 15, 20, -80),
            (34, 0.3, 15, 20, -80),
            (31, 0, 15, 20, -80),
            (31, -1, 15, 20, -80),
            (31, 1e15, 15, 20, -80),
            (31, 0.3, 15, -91, -80),
            (31, 0.3, 15, 20, -181),
            (31, 0.3, 100, 20, -80),
            (31, 0.3, 15, 20, -80, 2019, 2, 30),
            (31, 0.3, 15, 20, -80, 2019, 5, 27, -10),
            (31, 0.3, 15, 20, -80, 2019, 5, 27, 5),
            (31, 0.3, 15, 20, -80, 2019, 5, 27, 0.15, "bogus"),
            (31, 0.3, 15, 20, -80, 2019, 5, 27, 0.15, "standard", -1),
        )
        for arguments in invalid_argument_sets:
            with self.subTest(arguments=arguments):
                process = self.run_model(*arguments)
                self.assertNotEqual(process.returncode, 0)
                self.assertNotIn("RESULT_CSV,", process.stdout)

    def test_unavailable_date_fails_but_solar_override_succeeds(self) -> None:
        unavailable = self.run_model(31, 0.3, 15, 20, -80, 2024, 7, 29)
        self.assertNotEqual(unavailable.returncode, 0)

        overridden = self.run_model(
            31, 0.3, 15, 20, -80, 2024, 7, 29, 0.15, "standard", 50
        )
        self.assertEqual(overridden.returncode, 0, overridden.stderr)
        result = parse_result(overridden.stdout)
        self.assertEqual(result.solar_status, 0)
        self.assertAlmostEqual(result.solar_w, 50.0)

    def test_historical_annual_fallback_is_used(self) -> None:
        process = self.run_model(31, 0.3, 15, 20, -80, 1900, 7, 1)
        self.assertEqual(process.returncode, 0, process.stderr)
        result = parse_result(process.stdout)
        self.assertEqual(result.solar_status, 2)
        self.assertAlmostEqual(result.solar_w, 19.2)

    def test_photon_511_line_is_included_only_when_threshold_covers_it(self) -> None:
        included = parse_result(self.run_model(33, 0.3, 15, 20, -80).stdout)
        self.assertGreater(included.line_511_flux_cm2_s, 0.0)
        self.assertAlmostEqual(
            included.total_flux_cm2_s,
            included.continuum_flux_cm2_s + included.line_511_flux_cm2_s,
            places=12,
        )

        excluded_process = self.run_model(33, 1.0, 15, 20, -80)
        self.assertEqual(excluded_process.returncode, 0, excluded_process.stderr)
        excluded = parse_result(excluded_process.stdout)
        self.assertEqual(excluded.line_511_flux_cm2_s, 0.0)
        self.assertAlmostEqual(excluded.total_flux_cm2_s, excluded.continuum_flux_cm2_s, places=12)

    def test_atmosphere_and_geometry_controls_change_the_model(self) -> None:
        standard = parse_result(self.run_model(31, 0.3, 15, 20, -80).stdout)
        msis_process = self.run_model(31, 0.3, 15, 20, -80, 2019, 5, 27, 0.15, "msis")
        self.assertEqual(msis_process.returncode, 0, msis_process.stderr)
        msis = parse_result(msis_process.stdout)
        self.assertEqual(msis.atmosphere_model, 1)
        self.assertNotEqual(msis.atmospheric_depth_g_cm2, standard.atmospheric_depth_g_cm2)

        ground = parse_result(self.run_model(0, 0.3, 15, 20, -80).stdout)
        ideal_process = self.run_model(0, 0.3, 15, 20, -80, 2019, 5, 27, 10, "standard")
        self.assertEqual(ideal_process.returncode, 0, ideal_process.stderr)
        ideal = parse_result(ideal_process.stdout)
        self.assertNotEqual(ground.total_flux_cm2_s, ideal.total_flux_cm2_s)

    def test_python_wrapper_rejects_invalid_timeout(self) -> None:
        with self.assertRaises(ValueError):
            run_executable_with_parameters(
                31,
                0.3,
                15,
                20,
                -80,
                executable=EXECUTABLE,
                timeout=0,
            )

    def test_python_wrapper_is_cwd_independent(self) -> None:
        original_cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as workdir:
            os.chdir(workdir)
            try:
                # A relative executable path must anchor to the project root,
                # and the run must succeed away from the repository.
                result = run_executable_with_parameters(
                    31,
                    0.3,
                    15,
                    20,
                    -80,
                    executable=Path("electron_fluxes"),
                    timeout=60,
                )
            finally:
                os.chdir(original_cwd)
        self.assertGreater(result.total_flux_cm2_s, 0.0)

    def test_grid_pipeline_writes_csv(self) -> None:
        with tempfile.TemporaryDirectory() as workdir:
            output = Path(workdir) / "grid.csv"
            exit_code = main(
                [
                    "--thresholds", "0.3",
                    "--altitude-min", "15",
                    "--altitude-max", "15.2",
                    "--altitude-step", "0.2",
                    "--no-plot",
                    "--jobs", "2",
                    "--output", str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
        self.assertEqual(len(rows), 3)
        self.assertEqual(len(rows[0]), 19)
        self.assertEqual([float(row[4]) for row in rows[1:]], [15.0, 15.2])


if __name__ == "__main__":
    unittest.main()
