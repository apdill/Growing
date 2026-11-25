"""Run one simulation with the default configuration (growing.sim.run_sim.main): python scripts/run_single.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from growing.sim import run_sim  # noqa: E402

if __name__ == "__main__":
    run_sim.main()
