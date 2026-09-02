"""Command line interface for Phylo3D-Trait.

Commands:
  plot: Render interactive 3D phylogenetic tree with continuous trait evolution to HTML.
  template-values: Generate a CSV template of all tips and ancestral clades with stable IDs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo3d_trait.io import load_trait_values
from phylo3d_trait.renderer import build_figure, build_plot_data
from phylo3d_trait.template import generate_template_csv
from phylo3d_trait.tree import parse_tree


def build_parser() -> argparse.ArgumentParser:
    """Construct command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="phylo3d-trait",
        description="Universal interactive 3D phylogenetic tree tool with continuous trait evolution.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommand to run")

    # 1. 'plot' command
    plot_parser = subparsers.add_parser("plot", help="Generate interactive 3D HTML plot")
    plot_parser.add_argument(
        "--tree", "-t", required=True, type=str, help="Path to Newick or Nexus tree file"
    )
    plot_parser.add_argument(
        "--values", "-v", required=True, type=str, help="Path to CSV/TSV containing trait values"
    )
    plot_parser.add_argument(
        "--output", "-o", required=True, type=str, help="Path to save output standalone HTML file"
    )
    plot_parser.add_argument(
        "--title", type=str, default=None, help="Plot title"
    )
    plot_parser.add_argument(
        "--colorscale", type=str, default="Turbo", help="Plotly colorscale name (default: Turbo)"
    )
    plot_parser.add_argument(
        "--segments", "-s", type=int, default=10, help="Number of interpolation segments per branch (default: 10)"
    )
    plot_parser.add_argument(
        "--branch-width", type=float, default=1.0, help="Line width for 3D branch top outlines (default: 1.0)"
    )
    plot_parser.add_argument(
        "--no-labels", action="store_true", help="Do not render text labels for tips"
    )
    plot_parser.add_argument(
        "--baseline-y", type=float, default=None, help="Custom baseline Y trait plane height"
    )
    plot_parser.add_argument(
        "--baseline-raw-value", type=float, default=None,
        help="Custom numeric scientific trait label to display at baseline Y (default: raw_trait_max + 2 for reverse transform)"
    )
    plot_parser.add_argument(
        "--trait-display-range", nargs=2, type=float, default=None, metavar=("START", "END"),
        help="Linearly remap raw trait values [min, max] to custom display range [START, END] (e.g. 13 5)"
    )
    plot_parser.add_argument(
        "--opacity", type=float, default=1.0, help="Opacity for curtain meshes (0.0 - 1.0, default: 1.0)"
    )
    plot_parser.add_argument(
        "--no-mesh", action="store_true", help="Disable continuous curtain mesh surfaces"
    )
    plot_parser.add_argument(
        "--no-centerline", action="store_true", help="Disable branch top centerline outlines"
    )
    plot_parser.add_argument(
        "--centerline-color", type=str, default="dark", help="Centerline color ('dark', 'trait', or CSS color)"
    )
    plot_parser.add_argument(
        "--background", type=str, default="white", choices=["white", "transparent"],
        help="Background style ('white' or 'transparent', default: 'white')"
    )
    plot_parser.add_argument(
        "--camera-preset", type=str, default="elife", choices=["elife", "root_front", "tips_front"],
        help="Initial camera viewing preset (default: 'elife')"
    )
    plot_parser.add_argument(
        "--show-node-markers", action="store_true", help="Render diamond markers at ancestral nodes (default: False)"
    )

    # 2. 'template-values' command
    tpl_parser = subparsers.add_parser(
        "template-values", help="Generate CSV template listing all tip and ancestral node IDs"
    )
    tpl_parser.add_argument(
        "--tree", "-t", required=True, type=str, help="Path to Newick or Nexus tree file"
    )
    tpl_parser.add_argument(
        "--output", "-o", required=True, type=str, help="Path to save template CSV file"
    )
    tpl_parser.add_argument(
        "--default-val", type=str, default="", help="Default placeholder for trait column (default: '')"
    )

    return parser


def run_plot(args: argparse.Namespace) -> int:
    """Execute plot subcommand."""
    tree_path = Path(args.tree)
    if not tree_path.exists():
        print(f"Error: Tree file not found at '{tree_path}'", file=sys.stderr)
        return 1

    values_path = Path(args.values)
    if not values_path.exists():
        print(f"Error: Trait values file not found at '{values_path}'", file=sys.stderr)
        return 1

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        tree = parse_tree(tree_path)
        traits = load_trait_values(values_path)
        plot_data = build_plot_data(
            tree_input=tree,
            trait_values=traits,
            num_segments=args.segments,
            colorscale=args.colorscale,
            title=args.title if args.title else "3D Phylogenetic Tree with Continuous Trait Evolution",
            baseline_y=args.baseline_y,
            trait_display_range=args.trait_display_range,
            baseline_raw_value=args.baseline_raw_value,
        )
        fig = build_figure(
            plot_data=plot_data,
            title=args.title,
            branch_width=args.branch_width,
            show_tip_labels=not args.no_labels,
            show_mesh=not args.no_mesh,
            mesh_opacity=args.opacity,
            show_centerline=not args.no_centerline,
            centerline_color=args.centerline_color,
            baseline_y=args.baseline_y,
            baseline_raw_value=args.baseline_raw_value,
            show_node_markers=args.show_node_markers,
            background=args.background,
            camera_preset=args.camera_preset,
        )
        fig.write_html(str(out_path), include_plotlyjs="cdn", full_html=True)
        print(f"Successfully generated 3D phylogenetic visualization: {out_path}")
        return 0
    except Exception as e:
        print(f"Error during plotting: {e}", file=sys.stderr)
        return 1


def run_template(args: argparse.Namespace) -> int:
    """Execute template-values subcommand."""
    tree_path = Path(args.tree)
    if not tree_path.exists():
        print(f"Error: Tree file not found at '{tree_path}'", file=sys.stderr)
        return 1

    out_path = Path(args.output)
    try:
        generate_template_csv(
            tree_input=tree_path,
            output_path=out_path,
            default_value=args.default_val,
        )
        print(f"Successfully generated template CSV: {out_path}")
        return 0
    except Exception as e:
        print(f"Error generating template: {e}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    """Main CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "plot":
        return run_plot(args)
    elif args.command == "template-values":
        return run_template(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
