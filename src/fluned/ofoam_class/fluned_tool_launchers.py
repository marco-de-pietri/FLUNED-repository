"""foam_utils.py

Thin wrappers around common FLUNED/OpenFOAM execution commands.

"""

from __future__ import annotations

import logging
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Sequence, Union

# -----------------------------------------------------------------------------
# Logging - WARNING by default (silent). Elevate to INFO when verbose=True.
# -----------------------------------------------------------------------------
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
if not logger.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(_handler)

# -----------------------------------------------------------------------------
# Internal helper
# -----------------------------------------------------------------------------


def _run_case_command(
    case_dir: Union[str, Path],
    cmd: Sequence[str],
    *,
    log_file: bool = False,
    log_path: Union[str, Path, None] = None,
    verbose: bool = False,
) -> None:
    """Run *cmd* inside *case_dir*.

    Parameters
    ----------
    case_dir
        Path to the OpenFOAM case directory.
    cmd
        Command and arguments as a sequence of strings.
    log_file
        Write combined stdout/stderr to ``<command>_<timestamp>.log`` if *True*.
    log_path
        Append combined stdout/stderr to the provided log file path if given.
    verbose
        Stream child output to parent stdout/stderr and emit INFO-level log
        messages if *True*.
    """
    case_dir = Path(case_dir).expanduser().resolve()
    if not case_dir.is_dir():
        raise FileNotFoundError(f"{case_dir} is not a directory")

    if log_file and log_path is not None:
        raise ValueError("log_file and log_path cannot be used together")

    # Promote logger level when verbose output is requested.
    if verbose and logger.level > logging.INFO:
        logger.setLevel(logging.INFO)

    logger.info("Running `%s` in %s", " ".join(cmd), case_dir)

    if log_file:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = case_dir / f"{cmd[0]}_{ts}.log"

    if log_path is not None:
        resolved_log_path = Path(log_path)
        if not resolved_log_path.is_absolute():
            resolved_log_path = case_dir / resolved_log_path
        with resolved_log_path.open("a", encoding="utf-8") as out:
            subprocess.run(
                cmd, cwd=case_dir, stdout=out, stderr=subprocess.STDOUT, check=True
            )
        return

    stdout_target = None if verbose else subprocess.DEVNULL
    stderr_target = None if verbose else subprocess.DEVNULL
    subprocess.run(
        cmd, cwd=case_dir, stdout=stdout_target, stderr=stderr_target, check=True
    )


# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def launch_volume_func_object(
    path: Union[str, Path], *, log: bool = False, verbose: bool = False
) -> None:
    """Compute cell volumes via ``postProcess -func writeCellVolumes``."""
    _run_case_command(
        Path(path),
        ["postProcess", "-func", "writeCellVolumes"],
        log_file=log,
        verbose=verbose,
    )


def launch_centroid_func_object(
    path: Union[str, Path], *, log: bool = False, verbose: bool = False
) -> None:
    """Compute cell centres via ``postProcess -func writeCellCentres``."""
    _run_case_command(
        Path(path),
        ["postProcess", "-func", "writeCellCentres"],
        log_file=log,
        verbose=verbose,
    )


def launch_grad_func_object(
    path: Union[str, Path], *, log: bool = False, verbose: bool = False
) -> None:
    """Compute ?T via ``postProcess -func grad(T)``."""
    _run_case_command(
        Path(path), ["postProcess", "-func", "grad(T)"], log_file=log, verbose=verbose
    )


def launch_fluned_solver(
    path: Union[str, Path], *, log: bool = True, verbose: bool = False
) -> None:
    """Run ``FLUNED-solver`` inside the case directory."""
    _run_case_command(
        Path(path),
        ["FLUNED-solver"],
        log_path="simulation_log" if log else None,
        verbose=verbose,
    )


def launch_foam_to_vtk(
    path: Union[str, Path], *, log: bool = False, verbose: bool = False
) -> None:
    """Export latest timestep fields to VTK using ``foamToVTK``.

    OpenFOAM expects the regex after ``-excludePatches`` to be wrapped in
    parentheses. Passing ``".*"`` without the outer pair results in exit status
    1, hence the ``(".*")`` token.
    """
    _run_case_command(
        Path(path),
        [
            "foamToVTK",
            "-latestTime",
            "-noFaceZones",
            "-noFunctionObjects",
            "-fields",
            "(T Ta Td)",
            "-excludePatches",
            '(".*")',  # keep parentheses to satisfy foamToVTK
        ],
        log_file=log,
        verbose=verbose,
    )


# -----------------------------------------------------------------------------
# Re-export list
# -----------------------------------------------------------------------------
__all__: list[str] = [
    "launch_fluned_solver",
    "launch_volume_func_object",
    "launch_centroid_func_object",
    "launch_grad_func_object",
    "launch_foam_to_vtk",
]
