#!/usr/bin/env python3
"""
Streamlit GUI for plotlyMol3D - Visual Testing & Demo App.

Run with:
    streamlit run examples/gui_app.py

Requirements:
    pip install streamlit
"""
from pathlib import Path
import sys

try:
    from plotlymol3d.app import main
except ModuleNotFoundError:
    repo_root = Path(__file__).resolve().parents[1]
    src_path = repo_root / "src"
    sys.path.insert(0, str(src_path))
    from plotlymol3d.app import main

main()
