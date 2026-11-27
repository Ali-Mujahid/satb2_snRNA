"""Compatibility entrypoint for the current CA1 and inhibitory marker readers.

The earlier top-two-only builder is preserved in analysis/ca1_top_markers_20261127/source_v2/.
"""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).with_name('build_marker_explorers.py')), run_name='__main__')
