"""Run a long task in a child process, streaming its output.

Both slow operations in qaview — recruitment and figure export — are minutes-scale,
must not block the Panel server, and are far easier to reason about out of process:
SQLite connections stay on the thread that opened them, matplotlib never renders on a
server thread, progress is just the child's stdout, and Cancel is a real terminate().
"""

import os
import subprocess
import sys
import threading
from collections import deque


class SubprocessJob:
    """Poll `lines`, `done` and `error`; override `_on_success` to collect a result."""

    def __init__(self, argv: list, python: str | None = None):
        self.argv = argv
        self.lines: deque = deque(maxlen=400)
        self.done = False
        self.error: str | None = None
        self._proc: subprocess.Popen | None = None
        self._python = python or sys.executable

    def start(self) -> "SubprocessJob":
        env = dict(os.environ)
        env["PYTHONUNBUFFERED"] = "1"        # so progress arrives while it runs
        self._proc = subprocess.Popen(
            [self._python, *self.argv],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env,
        )
        threading.Thread(target=self._pump, daemon=True).start()
        return self

    def _pump(self):
        assert self._proc is not None and self._proc.stdout is not None
        for line in self._proc.stdout:
            self.lines.append(line.rstrip())
        code = self._proc.wait()
        if code != 0:
            self.error = f"{self.label} exited with code {code}"
        else:
            try:
                self._on_success()
            except Exception as exc:
                self.error = f"{self.label} produced no usable result: {exc}"
        self.done = True

    def _on_success(self):
        """Collect whatever the child wrote.  Raise to report failure."""

    def cancel(self):
        if self._proc is not None and self._proc.poll() is None:
            self._proc.terminate()
            self.lines.append("[cancelled]")
        self.done = True

    @property
    def label(self) -> str:
        return type(self).__name__

    @property
    def log(self) -> str:
        return "\n".join(self.lines)
