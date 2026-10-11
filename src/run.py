"""Command-line entry point: python -m src.run --help."""

import argparse
from pathlib import Path
import sys

from .core import ROOT
from .report import RunConfig, run_pipeline


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description="Generate and analyze synthetic data; entirely offline.")
    result.add_argument("--n", type=int, default=3_000, help="Number of generated people")
    result.add_argument("--seed", type=int, default=20260907)
    result.add_argument("--bootstrap", type=int, default=1_000, help="Resamples (minimum 100)")
    result.add_argument("--bootstrap-seed", type=int, default=8675309)
    result.add_argument("--simulations", type=int, default=100, help="Consecutive simulation seeds; 0 to skip")
    result.add_argument("--simulation-seed", type=int, default=20260908)
    result.add_argument("--output-dir", type=Path, default=ROOT / "results" / "local")
    result.add_argument("--overwrite", action="store_true", help="Explicitly replace generated files in the output folder")
    return result


def main() -> None:
    # Windows redirected output may default to a code page that cannot encode
    # non-ASCII workspace paths. The CLI emits a consistent UTF-8 text stream.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")
    arguments = parser().parse_args()
    try:
        config = RunConfig(arguments.n, arguments.seed, arguments.bootstrap,
                           arguments.bootstrap_seed, arguments.simulations, arguments.simulation_seed)
        metadata = run_pipeline(arguments.output_dir, config, arguments.overwrite)
    except (ValueError, OSError) as exc:
        raise SystemExit(f"Analysis stopped: {exc}") from exc
    print(f"Synthetic analysis complete: {arguments.output_dir.resolve() / 'REPORT.md'}")
    print(f"IPTW risk difference: {metadata['summary']['weighted_risk_difference']:.6f}")
    print("No clinical claim; inspect diagnostics, limitations, and the manifest.")


if __name__ == "__main__":
    main()
