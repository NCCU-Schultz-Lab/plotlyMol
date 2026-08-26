"""
Render harness for viewing plotlyMol output without a browser.

plotlyMol figures are normally viewed interactively (`fig.show()`), which
doesn't work in a headless/cloud session. This script rasterizes a figure
to a PNG via Plotly's kaleido image export, so it can be inspected with an
image viewer (or Claude Code's Read tool) instead.

Usage:
    # Single molecule
    python examples/render_preview.py --smiles "c1ccccc1" --mode ball+stick

    # Side-by-side comparison, e.g. checking a resolution or radius change
    python examples/render_preview.py --compare \\
        --smiles "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C" --resolution 16 --title "res=16" \\
        --smiles "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C" --resolution 32 --title "res=32"

    # From a file instead of SMILES
    python examples/render_preview.py --xyz path/to/molecule.xyz

Requires kaleido (already a project dependency). In sandboxed/cloud
environments without system Chrome, point BROWSER_PATH at a local Chromium
binary before running, or let this script auto-detect one under
/opt/pw-browsers (the layout used by Claude Code's cloud sessions).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from plotly.subplots import make_subplots  # noqa: E402

from plotlymol3d import draw_3D_rep  # noqa: E402

DEFAULT_CAMERA_EYE = (1.6, 1.6, 1.0)


def _ensure_browser_path() -> None:
    """Point kaleido at a pre-installed Chromium if available.

    Kaleido v1+ requires a real Chrome/Chromium binary and otherwise tries
    to download one, which fails without network access. If BROWSER_PATH
    isn't already set, look for the Chromium bundled in Claude Code's cloud
    sessions (PLAYWRIGHT_BROWSERS_PATH, conventionally /opt/pw-browsers).
    """
    if os.environ.get("BROWSER_PATH"):
        return

    search_roots = []
    env_root = os.environ.get("PLAYWRIGHT_BROWSERS_PATH")
    if env_root:
        search_roots.append(Path(env_root))
    search_roots.append(Path("/opt/pw-browsers"))

    for root in search_roots:
        if not root.is_dir():
            continue
        candidates = sorted(root.glob("chromium-*/chrome-linux/chrome"))
        if candidates:
            os.environ["BROWSER_PATH"] = str(candidates[-1])
            return


def _panel_layout(bgcolor: str, camera_eye) -> dict[str, Any]:
    return {
        "bgcolor": bgcolor,
        "camera": {"eye": {"x": camera_eye[0], "y": camera_eye[1], "z": camera_eye[2]}},
        # Without aspectmode="data", each scene stretches x/y/z independently
        # to fill its (non-square) subplot domain, badly distorting geometry.
        "aspectmode": "data",
        "xaxis": {"visible": False},
        "yaxis": {"visible": False},
        "zaxis": {"visible": False},
    }


def render_molecule(
    out: str,
    width: int = 900,
    height: int = 700,
    bgcolor: str = "white",
    camera_eye=DEFAULT_CAMERA_EYE,
    **draw_kwargs: Any,
) -> Path:
    """Render a single molecule (via draw_3D_rep) to a PNG file."""
    _ensure_browser_path()
    fig = draw_3D_rep(**draw_kwargs)
    fig.update_layout(
        width=width,
        height=height,
        paper_bgcolor=bgcolor,
        scene=_panel_layout(bgcolor, camera_eye),
        margin={"l": 0, "r": 0, "t": 10, "b": 0},
    )
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(str(out_path))
    return out_path


def render_comparison(
    panels: list[dict[str, Any]],
    out: str,
    panel_width: int = 480,
    height: int = 520,
    bgcolor: str = "white",
    camera_eye=DEFAULT_CAMERA_EYE,
) -> Path:
    """Render several draw_3D_rep configurations side by side in one PNG.

    Args:
        panels: one dict per panel, e.g. {"title": "res=16", "smiles": "CCO",
            "resolution": 16, "mode": "ball+stick"}. All keys except
            "title" are forwarded to draw_3D_rep.
        out: output PNG path.
    """
    _ensure_browser_path()
    n = len(panels)
    fig = make_subplots(
        rows=1,
        cols=n,
        specs=[[{"type": "scene"}] * n],
        subplot_titles=[p.get("title", f"panel {i + 1}") for i, p in enumerate(panels)],
        horizontal_spacing=0.02,
    )

    for i, panel in enumerate(panels):
        kwargs = {k: v for k, v in panel.items() if k != "title"}
        sub_fig = draw_3D_rep(**kwargs)
        for trace in sub_fig.data:
            fig.add_trace(trace, row=1, col=i + 1)

        scene_key = "scene" if i == 0 else f"scene{i + 1}"
        fig.layout[scene_key].update(_panel_layout(bgcolor, camera_eye))

    fig.update_layout(
        width=panel_width * n,
        height=height,
        paper_bgcolor=bgcolor,
        margin={"l": 0, "r": 0, "t": 40, "b": 0},
    )
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(str(out_path))
    return out_path


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--smiles",
        action="append",
        default=[],
        help="SMILES string. Repeat with --compare for multiple panels.",
    )
    parser.add_argument(
        "--xyz", action="append", default=[], help="Path to an XYZ file."
    )
    parser.add_argument(
        "--mol", action="append", default=[], help="Path to a MOL file."
    )
    parser.add_argument(
        "--mode", action="append", default=[], help="ball+stick, ball, stick, or vdw."
    )
    parser.add_argument(
        "--resolution", action="append", type=int, default=[], help="Mesh resolution."
    )
    parser.add_argument(
        "--title", action="append", default=[], help="Panel title (compare mode only)."
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Render one panel per --smiles/--xyz/--mol entry, side by side.",
    )
    parser.add_argument("--out", default=None, help="Output PNG path.")
    parser.add_argument(
        "--width", type=int, default=900, help="Image width (single-render mode)."
    )
    parser.add_argument("--height", type=int, default=700, help="Image height.")
    parser.add_argument("--bgcolor", default="white", help="Background color.")
    return parser


def _pad(values: list[Any], n: int, default: Any) -> list[Any]:
    return [values[i] if i < len(values) else default for i in range(n)]


def main(argv: list[str] | None = None) -> None:
    args = _build_arg_parser().parse_args(argv)

    inputs: list[dict[str, Any]] = []
    for s in args.smiles:
        inputs.append({"smiles": s})
    for x in args.xyz:
        inputs.append({"xyzfile": x})
    for m in args.mol:
        inputs.append({"molfile": m})

    if not inputs:
        raise SystemExit("Provide at least one --smiles, --xyz, or --mol.")

    n = len(inputs)
    modes = _pad(args.mode, n, "ball+stick")
    resolutions = _pad(args.resolution, n, 32)
    titles = _pad(args.title, n, None)

    for i, inp in enumerate(inputs):
        inp["mode"] = modes[i]
        inp["resolution"] = resolutions[i]
        if titles[i] is not None:
            inp["title"] = titles[i]

    if args.compare or n > 1:
        out = args.out or "compare.png"
        for i, inp in enumerate(inputs):
            inp.setdefault(
                "title",
                inp.get("smiles")
                or inp.get("xyzfile")
                or inp.get("molfile")
                or f"panel {i + 1}",
            )
        out_path = render_comparison(inputs, out=out, bgcolor=args.bgcolor)
    else:
        out = args.out or "preview.png"
        inp = inputs[0]
        inp.pop("title", None)
        out_path = render_molecule(
            out=out, width=args.width, height=args.height, bgcolor=args.bgcolor, **inp
        )

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
