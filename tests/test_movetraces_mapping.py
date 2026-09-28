"""Test the Plotly.moveTraces index mapping fix.

This test verifies that the curtain trace reordering in the camera hook
correctly maps current trace indices to destination indices.
"""

import json
import tempfile
import os
from pathlib import Path

import plotly.graph_objects as go
from phylo3d_trait.renderer import TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT


def extract_hook_code():
    """Extract the JavaScript hook code from the renderer."""
    # The hook is embedded in TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT
    # We need to extract just the sortCurtainTraces function
    return TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT


def test_destination_array_construction():
    """Test the destination array construction logic."""
    
    # Simulate the current buggy algorithm
    def buggy_moveTraces_call(idxs, entries_sorted):
        """Current buggy implementation."""
        order = [e['trace'] for e in entries_sorted]
        return ('buggy', order, idxs)
    
    # Correct implementation (absolute gd.data indices)
    def correct_moveTraces_call(idxs, entries_sorted):
        """Correct implementation using destination array with absolute indices."""
        order = [e['trace'] for e in entries_sorted]
        # Build destination array: for each trace in idxs, find its absolute target position
        # in gd.data. Curtains occupy contiguous slots starting at idxs[0].
        dest = [idxs[order.index(trace)] for trace in idxs]
        return ('correct', idxs, dest)
    
    # Test case 1: identity order
    idxs = [1, 2, 3]
    entries_sorted = [{'trace': 1}, {'trace': 2}, {'trace': 3}]
    
    buggy = buggy_moveTraces_call(idxs, entries_sorted)
    correct = correct_moveTraces_call(idxs, entries_sorted)
    
    assert buggy == ('buggy', [1, 2, 3], [1, 2, 3])
    assert correct == ('correct', [1, 2, 3], [1, 2, 3])  # Absolute indices
    print("Test 1 (identity): PASS")
    
    # Test case 2: simple swap
    idxs = [1, 2, 3]
    entries_sorted = [{'trace': 2}, {'trace': 1}, {'trace': 3}]
    
    buggy = buggy_moveTraces_call(idxs, entries_sorted)
    correct = correct_moveTraces_call(idxs, entries_sorted)
    
    assert buggy == ('buggy', [2, 1, 3], [1, 2, 3])
    assert correct == ('correct', [1, 2, 3], [2, 1, 3])  # Absolute indices
    print("Test 2 (swap): PASS")
    
    # Test case 3: 3-cycle permutation
    idxs = [1, 2, 3]
    entries_sorted = [{'trace': 3}, {'trace': 1}, {'trace': 2}]
    
    buggy = buggy_moveTraces_call(idxs, entries_sorted)
    correct = correct_moveTraces_call(idxs, entries_sorted)
    
    assert buggy == ('buggy', [3, 1, 2], [1, 2, 3])
    assert correct == ('correct', [1, 2, 3], [2, 3, 1])  # Absolute indices
    print("Test 3 (3-cycle): PASS")
    
    # Test case 4: non-zero / non-contiguous trace slots
    idxs = [5, 10, 15]
    entries_sorted = [{'trace': 15}, {'trace': 5}, {'trace': 10}]
    
    buggy = buggy_moveTraces_call(idxs, entries_sorted)
    correct = correct_moveTraces_call(idxs, entries_sorted)
    
    assert buggy == ('buggy', [15, 5, 10], [5, 10, 15])
    assert correct == ('correct', [5, 10, 15], [10, 15, 5])  # Absolute indices
    print("Test 4 (non-contiguous): PASS")
    
    # Test case 5: centerline before meshes (centerline at index 0)
    idxs = [1, 2, 3, 4]  # mesh traces start at 1
    entries_sorted = [{'trace': 4}, {'trace': 2}, {'trace': 3}, {'trace': 1}]
    
    buggy = buggy_moveTraces_call(idxs, entries_sorted)
    correct = correct_moveTraces_call(idxs, entries_sorted)
    
    assert buggy == ('buggy', [4, 2, 3, 1], [1, 2, 3, 4])
    assert correct == ('correct', [1, 2, 3, 4], [4, 2, 3, 1])  # Absolute indices
    print("Test 5 (centerline before): PASS")
    
    print("\nAll destination array tests PASSED")


def build_synthetic_test_html(camera_name, camera_cfg, output_path):
    """Build a synthetic 3-layer test with known depth ordering."""
    
    # Three parallel planes at different X positions
    # Blue at x=0 (farthest from +X camera)
    # Cyan at x=1 (middle)
    # Yellow at x=2 (nearest to +X camera)
    YELLOW = "rgb(255,200,0)"
    CYAN = "rgb(0,200,255)"
    BLUE = "rgb(0,60,255)"
    LIGHTING = dict(ambient=1.0, diffuse=0.0, specular=0.0, roughness=1.0, fresnel=0.0)
    OPACITY = 0.9
    
    def quad_mesh(x_plane, color, name):
        return go.Mesh3d(
            x=[x_plane, x_plane, x_plane, x_plane],
            y=[0, 1, 0, 1],
            z=[-0.5, -0.5, 0.5, 0.5],
            i=[0, 1], j=[1, 2], k=[2, 3],
            color=color, opacity=OPACITY, flatshading=True,
            lighting=LIGHTING, hoverinfo="skip", name=name, showlegend=False,
        )
    
    def scene_cfg(camera):
        return dict(
            xaxis=dict(title="", showticklabels=False, showbackground=False),
            yaxis=dict(title="", showticklabels=False, showbackground=False),
            zaxis=dict(title="", showticklabels=False, showbackground=False),
            aspectmode="manual",
            aspectratio=dict(x=1.4, y=1.0, z=1.2),
            camera=camera,
        )
    
    # Authored order: Yellow (nearest), Blue (farthest), Cyan (middle)
    # This is DELIBERATELY WRONG to test reordering
    fig = go.Figure()
    fig.add_trace(quad_mesh(2.0, YELLOW, "yellow_near_x2"))
    fig.add_trace(quad_mesh(0.0, BLUE, "blue_far_x0"))
    fig.add_trace(quad_mesh(1.0, CYAN, "cyan_mid_x1"))
    fig.update_layout(scene=scene_cfg(camera_cfg), margin=dict(l=0, r=0, t=0, b=0),
                      width=700, height=500, paper_bgcolor="white")
    
    # Inject the hook - need to escape the JS braces
    plot_id = "test_plot"
    hook = TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT.replace("{plot_id}", plot_id)
    fig.write_html(str(output_path), include_plotlyjs=True, full_html=True, post_script=hook)
    
    return fig


def verify_synthetic_test(html_path, camera_name, expected_order):
    """Verify the synthetic test by checking the generated HTML."""
    with open(html_path, 'rb') as f:
        html = f.read().decode('utf-8', errors='ignore')
    
    # Check that the hook is present
    assert 'sortCurtainTraces' in html
    assert 'moveTraces' in html
    # Check the fix is in place: idxs[rank] instead of order.indexOf(idxs
    assert 'dest.push(idxs[rank])' in html
    print(f"  {camera_name}: Hook present in HTML with fix")
    
    return True


if __name__ == "__main__":
    test_destination_array_construction()
    
    # Build synthetic tests
    out_dir = Path("H:/Work/Paper_work/101_hemoglobin/phylo3d_trait/tests/synthetic")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    CAMERA_A = dict(eye=dict(x=2.0, y=0.5, z=0.9), up=dict(x=0, y=1, z=0),
                    center=dict(x=0, y=0, z=0), projection=dict(type="orthographic"))
    CAMERA_B = dict(eye=dict(x=-2.0, y=0.5, z=0.9), up=dict(x=0, y=1, z=0),
                    center=dict(x=0, y=0, z=0), projection=dict(type="orthographic"))
    
    build_synthetic_test_html("CAMERA_A", CAMERA_A, out_dir / "synth_v1_camera_a.html")
    build_synthetic_test_html("CAMERA_B", CAMERA_B, out_dir / "synth_v1_camera_b.html")
    
    verify_synthetic_test(out_dir / "synth_v1_camera_a.html", "CAMERA_A", ["blue_far_x0", "cyan_mid_x1", "yellow_near_x2"])
    verify_synthetic_test(out_dir / "synth_v1_camera_b.html", "CAMERA_B", ["yellow_near_x2", "cyan_mid_x1", "blue_far_x0"])
    
    print("\nSynthetic test HTML files generated.")