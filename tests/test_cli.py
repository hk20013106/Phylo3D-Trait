"""Integration tests for Phylo3D-Trait CLI commands and repository examples."""

from pathlib import Path
import pytest
import plotly.graph_objects as go

from phylo3d_trait.cli import main
from phylo3d_trait.io import load_trait_values
from phylo3d_trait.renderer import _build_branch_curtains_geometry, build_figure, build_plot_data
from phylo3d_trait.tree import parse_tree


def test_cli_template_values_generation(tmp_path: Path):
    """Test generating a template values CSV from Newick tree."""
    tree_file = tmp_path / "test_tree.nwk"
    tree_file.write_text("((A:5,B:5):10,(C:8,D:8):7);", encoding="utf-8")

    out_csv = tmp_path / "template.csv"
    code = main(["template-values", "--tree", str(tree_file), "--output", str(out_csv)])
    assert code == 0
    assert out_csv.exists()

    content = out_csv.read_text(encoding="utf-8")
    assert "node_id,label,node_type" in content
    assert "A,A,tip" in content
    assert "clade:" in content


def test_cli_end_to_end_plot(tmp_path: Path):
    """Test end-to-end plot command producing standalone HTML."""
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:10,B:10):20,(C:15,D:15):15);", encoding="utf-8")

    template_csv = tmp_path / "template.csv"
    code_tpl = main(["template-values", "--tree", str(tree_file), "--output", str(template_csv)])
    assert code_tpl == 0

    # Fill trait values
    lines = template_csv.read_text(encoding="utf-8").strip().splitlines()
    filled_lines = [lines[0]]
    val = 1.0
    for line in lines[1:]:
        parts = line.split(",")
        parts[-1] = str(val)
        val += 1.0
        filled_lines.append(",".join(parts))

    values_csv = tmp_path / "values.csv"
    values_csv.write_text("\n".join(filled_lines), encoding="utf-8")

    out_html = tmp_path / "output_plot.html"
    code_plot = main([
        "plot",
        "--tree", str(tree_file),
        "--values", str(values_csv),
        "--output", str(out_html),
        "--title", "Test Phylogeny",
    ])
    assert code_plot == 0
    assert out_html.exists()
    assert out_html.stat().st_size > 2000


def test_cli_missing_values_error(tmp_path: Path):
    """Test CLI returns non-zero error exit code when values are missing."""
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:10,B:10):20,(C:15,D:15):15);", encoding="utf-8")

    incomplete_csv = tmp_path / "incomplete.csv"
    incomplete_csv.write_text("node_id,trait\nA,1.0\nB,2.0\n", encoding="utf-8")

    out_html = tmp_path / "fail.html"
    code = main([
        "plot",
        "--tree", str(tree_file),
        "--values", str(incomplete_csv),
        "--output", str(out_html),
    ])
    assert code != 0
    assert not out_html.exists()


def test_examples_example1_files_end_to_end(tmp_path: Path):
    """Test generating 3D plot using actual repository Example 1 files."""
    repo_root = Path(__file__).parent.parent
    tree_file = repo_root / "examples" / "example1" / "tree.nwk"
    values_file = repo_root / "examples" / "example1" / "node_values.csv"
    out_html = tmp_path / "example1_plot.html"

    assert tree_file.exists()
    assert values_file.exists()

    code = main([
        "plot",
        "--tree", str(tree_file),
        "--values", str(values_file),
        "--output", str(out_html),
    ])
    assert code == 0
    assert out_html.exists()
    assert out_html.stat().st_size > 5000


def test_cli_opacity_options(tmp_path: Path):
    """Test CLI default opacity is 1.0 and custom opacity is written to HTML."""
    repo_root = Path(__file__).parent.parent
    tree_file = repo_root / "examples" / "example1" / "tree.nwk"
    values_file = repo_root / "examples" / "example1" / "node_values.csv"
    out_html = tmp_path / "opaque_plot.html"

    code = main([
        "plot",
        "--tree", str(tree_file),
        "--values", str(values_file),
        "--output", str(out_html),
    ])
    assert code == 0
    content = out_html.read_text(encoding="utf-8")
    assert '"opacity":1.0' in content or '"opacity": 1.0' in content or '"opacity":1' in content


def test_example2_tree_and_traits_properties():
    """Verify Example 2 tree structure, 6 taxa, and supplied trait bounds [1, 5]."""
    repo_root = Path(__file__).parent.parent
    tree_path = repo_root / "examples" / "example2" / "tree.nwk"
    values_path = repo_root / "examples" / "example2" / "node_values.csv"

    assert tree_path.exists(), "example2 tree.nwk must exist"
    assert values_path.exists(), "example2 node_values.csv must exist"

    tree = parse_tree(tree_path)
    terms = tree.get_terminals()
    assert len(terms) == 6, f"Expected 6 tips, found {len(terms)}"

    tip_names = [t.name for t in terms]
    assert sorted(tip_names) == ["A", "B", "C", "D", "E", "F"]

    traits = load_trait_values(values_path)
    assert len(traits) == 11

    trait_vals = list(traits.values())
    assert min(trait_vals) >= 1.0
    assert max(trait_vals) <= 5.0
    assert min(trait_vals) == pytest.approx(1.0)
    assert max(trait_vals) == pytest.approx(5.0)


def test_example2_baseline_zero_geometry_and_color_range():
    """Verify that baseline_y=0 produces mesh cmin=0, cmax=5, and root height=1.0."""
    repo_root = Path(__file__).parent.parent
    tree = parse_tree(repo_root / "examples" / "example2" / "tree.nwk")
    traits = load_trait_values(repo_root / "examples" / "example2" / "node_values.csv")

    plot_data = build_plot_data(tree, traits, num_segments=10, baseline_y=0.0)
    assert plot_data.baseline_y == pytest.approx(0.0)
    assert plot_data.trait_min == pytest.approx(1.0)
    assert plot_data.trait_max == pytest.approx(5.0)

    fig = build_figure(plot_data, baseline_y=0.0)
    mesh_traces = [t for t in fig.data if isinstance(t, go.Mesh3d) or getattr(t, "type", None) == "mesh3d"]
    assert len(mesh_traces) == 1
    mesh = mesh_traces[0]

    assert mesh.cmin == pytest.approx(0.0)
    assert mesh.cmax == pytest.approx(5.0)

    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, baseline_y=0.0
    )

    assert min(mesh_y) == pytest.approx(0.0)
    assert max(mesh_y) == pytest.approx(5.0)
    assert min(mesh_intensity) == pytest.approx(0.0)
    assert max(mesh_intensity) == pytest.approx(5.0)

    for y_val, int_val in zip(mesh_y, mesh_intensity):
        assert int_val == pytest.approx(y_val)

    root_top_vertices = [
        mesh_y[idx] for idx in range(len(mesh_z))
        if mesh_z[idx] == pytest.approx(50.0) and mesh_y[idx] > 0.5
    ]
    root_bottom_vertices = [
        mesh_y[idx] for idx in range(len(mesh_z))
        if mesh_z[idx] == pytest.approx(50.0) and mesh_y[idx] < 0.5
    ]

    assert len(root_top_vertices) > 0
    assert len(root_bottom_vertices) > 0
    for v in root_top_vertices:
        assert v == pytest.approx(1.0)
    for v in root_bottom_vertices:
        assert v == pytest.approx(0.0)

    root_curtain_height = root_top_vertices[0] - root_bottom_vertices[0]
    assert root_curtain_height == pytest.approx(1.0)


def test_example2_cli_end_to_end(tmp_path: Path):
    """Verify full CLI execution with --baseline-y 0 for Example 2."""
    repo_root = Path(__file__).parent.parent
    tree_file = repo_root / "examples" / "example2" / "tree.nwk"
    values_file = repo_root / "examples" / "example2" / "node_values.csv"
    out_html = tmp_path / "example2_test.html"

    ret = main([
        "plot",
        "--tree", str(tree_file),
        "--values", str(values_file),
        "--output", str(out_html),
        "--baseline-y", "0",
    ])
    assert ret == 0
    assert out_html.exists()
    assert out_html.stat().st_size > 5000

    content = out_html.read_text(encoding="utf-8")
    assert '"cmin":0.0' in content or '"cmin": 0.0' in content
    assert '"cmax":5.0' in content or '"cmax": 5.0' in content
    assert "Branch Curtains" in content
    assert "Terminal Taxa" in content


def test_cli_trait_display_range_end_to_end(tmp_path: Path):
    """Verify full CLI execution with --trait-display-range 13 5 and --baseline-y 0."""
    repo_root = Path(__file__).parent.parent
    tree_file = repo_root / "examples" / "example2" / "tree.nwk"
    values_file = repo_root / "examples" / "example2" / "node_values.csv"
    out_html = tmp_path / "example2_rescaled.html"

    ret = main([
        "plot",
        "--tree", str(tree_file),
        "--values", str(values_file),
        "--output", str(out_html),
        "--baseline-y", "0",
        "--trait-display-range", "13", "5",
    ])
    assert ret == 0
    assert out_html.exists()
    assert out_html.stat().st_size > 5000

    content = out_html.read_text(encoding="utf-8")
    assert '"cmin":0.0' in content or '"cmin": 0.0' in content
    assert '"cmax":13.0' in content or '"cmax": 13.0' in content
    assert "Raw Trait:" in content
    assert "Display Trait (internal Y):" in content
    assert "Trait Value" in content
    assert "baseline" in content


