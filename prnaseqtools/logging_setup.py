"""
Logging setup for pRNASeqTools.
Mirrors the Perl IO::Tee behavior: log to file + stderr simultaneously.
"""

import sys
import logging
import os


class TeeHandler:
    """Custom handler that writes to both a file and stderr."""

    def __init__(self, log_path):
        self.file = open(log_path, 'w')
        self.stderr = sys.stderr

    def write(self, message):
        self.file.write(message)
        self.file.flush()
        self.stderr.write(message)
        self.stderr.flush()

    def flush(self):
        self.file.flush()
        self.stderr.flush()

    def close(self):
        self.file.close()


def setup_logging(prefix, start_time):
    """Initialize logging: redirect stdout to log file + stderr tee."""
    log_path = os.path.join(os.getcwd(), f"log_{start_time}.txt")
    tee = TeeHandler(log_path)

    import prnaseqtools.cli as cli_mod
    cli_mod.LOG_FILE = log_path
    cli_mod.TEE = tee

    # Write header (mirrors Perl version)
    import time as _time
    localtime = _time.asctime(_time.localtime(start_time))
    tee.write(f"\npRNASeqTools\nVersion 1.0\tPython3 rewrite\n")
    tee.write(f"Command: {' '.join(sys.argv)}\n")
    tee.write(f"Start: {localtime}\n")
    tee.write("##########\n")

    return tee


def finalize_logging(tee, start_time, cleanup=False):
    """Write completion message and close log.  If *cleanup* is True the log
    file is deleted afterwards — use this for help-only or empty invocations
    where no meaningful analysis output was produced."""
    import time as _time
    log_path = tee.file.name if hasattr(tee, 'file') else None

    if not cleanup:
        localtime = _time.asctime(_time.localtime(_time.time()))
        tee.write(f"\nMission Completed!\nEnd: {localtime}\n")

    tee.close()

    if cleanup and log_path and os.path.exists(log_path):
        os.unlink(log_path)
