"""Repository locations used by scripts (independent of the current working directory)."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECS_DIR = os.path.join(ROOT, "specs")              # sweep specifications (tracked)
SWEEPS_DIR = os.path.join(ROOT, "sweeps")            # sweep outputs (ignored)
REPORTS_DIR = os.path.join(ROOT, "reports")          # figures / tables / write-ups (ignored)
EXTERNAL_DIR = os.path.join(ROOT, "external")        # downloaded third-party data (ignored)
SUMMARY_DIR = os.path.join(ROOT, "summary_plots")    # single-run outputs of run_sim (ignored)
