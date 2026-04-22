"""foam_utils.py

Thin wrappers around common FLUNED/OpenFOAM execution commands.

"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from contextlib import contextmanager
from datetime import datetime
from hashlib import sha1
from pathlib import Path
from typing import Iterator, Sequence, Union

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


def _path_contains_openfoam_invalid_chars(path: Path) -> bool:
    """Return True when OpenFOAM is likely to reject *path* as a fileName."""
    return any(char.isspace() for char in str(path))


@contextmanager
def _case_command_context(
    case_dir: Path, cmd: Sequence[str], *, openfoam_safe_case: bool = False
) -> Iterator[tuple[Path, list[str]]]:
    """Yield the working directory and command to use for an OpenFOAM case.

    OpenFOAM 12 aborts when the process working directory contains whitespace.
    For tools that support ``-case``, run from a temporary no-whitespace alias
    instead and pass the alias as the case path. Writes still land in the real
    case directory because the alias is a symlink to it.
    """
    if not openfoam_safe_case or not _path_contains_openfoam_invalid_chars(case_dir):
        yield case_dir, list(cmd)
        return

    digest = sha1(str(case_dir).encode("utf-8")).hexdigest()[:12]
    with tempfile.TemporaryDirectory(prefix="fluned_openfoam_case_") as tmp_dir:
        alias = Path(tmp_dir) / f"case_{digest}"
        alias.symlink_to(case_dir, target_is_directory=True)
        yield Path(tmp_dir), [cmd[0], "-case", str(alias), *cmd[1:]]


def _run_case_command(
    case_dir: Union[str, Path],
    cmd: Sequence[str],
    *,
    log_file: bool = False,
    log_path: Union[str, Path, None] = None,
    verbose: bool = False,
    openfoam_safe_case: bool = False,
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
    openfoam_safe_case
        For OpenFOAM tools that support ``-case``, avoid running from a case
        path with whitespace because OpenFOAM rejects those paths.
    """
    case_dir = Path(case_dir).expanduser().resolve()
    if not case_dir.is_dir():
        raise FileNotFoundError(f"{case_dir} is not a directory")

    if log_file and log_path is not None:
        raise ValueError("log_file and log_path cannot be used together")

    # Promote logger level when verbose output is requested.
    if verbose and logger.level > logging.INFO:
        logger.setLevel(logging.INFO)

    if log_file:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = case_dir / f"{cmd[0]}_{ts}.log"

    with _case_command_context(
        case_dir, cmd, openfoam_safe_case=openfoam_safe_case
    ) as (run_dir, run_cmd):
        logger.info("Running `%s` in %s", " ".join(run_cmd), run_dir)

        if log_path is not None:
            resolved_log_path = Path(log_path)
            if not resolved_log_path.is_absolute():
                resolved_log_path = case_dir / resolved_log_path
            with resolved_log_path.open("a", encoding="utf-8") as out:
                subprocess.run(
                    run_cmd,
                    cwd=run_dir,
                    stdout=out,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            return

        stdout_target = None if verbose else subprocess.DEVNULL
        stderr_target = None if verbose else subprocess.DEVNULL
        subprocess.run(
            run_cmd,
            cwd=run_dir,
            stdout=stdout_target,
            stderr=stderr_target,
            check=True,
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
        openfoam_safe_case=True,
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
        openfoam_safe_case=True,
    )


def launch_grad_func_object(
    path: Union[str, Path], *, log: bool = False, verbose: bool = False
) -> None:
    """Compute ?T via ``postProcess -func grad(T)``."""
    _run_case_command(
        Path(path),
        ["postProcess", "-func", "grad(T)"],
        log_file=log,
        verbose=verbose,
        openfoam_safe_case=True,
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
        openfoam_safe_case=True,
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
        openfoam_safe_case=True,
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
