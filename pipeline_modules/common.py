from __future__ import annotations

import shlex
import subprocess
import sys


def run_cmd(cmd: list[str]) -> None:
    print(f"\n>>> Running: {' '.join(shlex.quote(x) for x in cmd)}")
    subprocess.run(cmd, check=True)


def run_python_script(script: str, args: list[str] | None = None) -> None:
    cmd = [sys.executable, script]
    if args:
        cmd.extend(args)
    run_cmd(cmd)


def run_python_module(module_name: str, args: list[str] | None = None) -> None:
    cmd = [sys.executable, "-m", module_name]
    if args:
        cmd.extend(args)
    run_cmd(cmd)

