"""Fixed four-layer per-fragment transparency renderer.

Peels the nearest four fragments at each pixel. Layers 1-3 use the requested
standard opacity; layer 4 is the terminal opaque layer. Deeper fragments are
not evaluated. Supported opacity is 0.5..1.0 (transparency 0.0..0.5).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from plotly.colors import get_colorscale, sample_colorscale, unlabel_rgb

from phylo3d_trait.models import PlotData
from phylo3d_trait.renderer import (
    CAMERA_PRESETS,
    DEFAULT_TIP_LABEL_OFFSET_FRACTION,
    _build_branch_curtains_geometry,
    _generate_rescaled_ticks,
)

FOUR_LAYER_PEELS = 4
MIN_SUPPORTED_OPACITY = 0.5
MAX_SUPPORTED_TRANSPARENCY = 0.5


def _validate_opacity(opacity: float) -> float:
    try:
        value = float(opacity)
    except (TypeError, ValueError) as err:
        raise ValueError(f"opacity must be numeric, got {opacity!r}") from err
    if not math.isfinite(value) or not 0.5 <= value <= 1.0:
        raise ValueError(
            "four-layer renderer supports opacity in [0.5, 1.0] "
            "(transparency in [0.0, 0.5])"
        )
    return value


def _rgb01(rgb_text: str) -> List[float]:
    r, g, b = unlabel_rgb(rgb_text)
    return [float(r) / 255.0, float(g) / 255.0, float(b) / 255.0]


def _sample_colors(
    values: List[float],
    colorscale: str,
    cmin: float,
    cmax: float,
    reverse: bool,
) -> List[List[float]]:
    scale = get_colorscale(colorscale)
    span = cmax - cmin or 1.0
    positions = []
    for value in values:
        t = max(0.0, min(1.0, (float(value) - cmin) / span))
        positions.append(1.0 - t if reverse else t)
    return [_rgb01(x) for x in sample_colorscale(scale, positions, colortype="rgb")]


def _css_color(color: str) -> List[float]:
    color = color.strip()
    if color.startswith("#") and len(color) == 7:
        return [int(color[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]
    if color.lower().startswith("rgb"):
        return _rgb01(color)
    raise ValueError(
        "custom centerline color must be '#RRGGBB' or 'rgb(r,g,b)' "
        "for the four-layer renderer"
    )


def _scene_transform(
    plot_data: PlotData,
    baseline_y: float,
    trait_axis_scale: float,
):
    mins = [
        float(plot_data.x_min),
        min(float(plot_data.trait_min), baseline_y),
        float(plot_data.time_min),
    ]
    maxs = [
        float(plot_data.x_max),
        max(float(plot_data.trait_max), baseline_y),
        float(plot_data.time_max),
    ]
    centers = [(a + b) / 2.0 for a, b in zip(mins, maxs)]
    spans = [max(b - a, 1e-9) for a, b in zip(mins, maxs)]
    scales = [2.8 / spans[0], 2.0 * trait_axis_scale / spans[1], 2.4 / spans[2]]

    def convert(x: float, y: float, z: float) -> List[float]:
        return [
            (float(x) - centers[0]) * scales[0],
            (float(y) - centers[1]) * scales[1],
            (float(z) - centers[2]) * scales[2],
        ]

    return convert


def _build_payload(
    plot_data: PlotData,
    *,
    opacity: float,
    baseline_y: Optional[float],
    reverse_colorscale: bool,
    curtain_color_mode: str,
    trait_axis_scale: float,
    tip_label_offset: Optional[float],
    show_tip_labels: bool,
    show_centerline: bool,
    centerline_color: str,
    background: str,
    camera_preset: str,
    show_x_axis: bool = True,
    show_y_axis: bool = True,
    show_z_axis: bool = True,
    show_tip_hover: bool = True,
    show_internal_hover: bool = True,
) -> Dict[str, Any]:
    opacity = _validate_opacity(opacity)
    trait_axis_scale = float(trait_axis_scale)
    if not math.isfinite(trait_axis_scale) or trait_axis_scale <= 0:
        raise ValueError("trait_axis_scale must be finite and > 0")
    if curtain_color_mode not in {"height", "branch"}:
        raise ValueError("curtain_color_mode must be 'height' or 'branch'")
    if background not in {"white", "transparent"}:
        raise ValueError("background must be 'white' or 'transparent'")

    baseline = float(
        baseline_y
        if baseline_y is not None
        else (
            plot_data.baseline_y
            if plot_data.baseline_y is not None
            else plot_data.trait_min
        )
    )
    trait_cmin = float(plot_data.trait_min)
    trait_cmax = float(plot_data.trait_max)
    if trait_cmin == trait_cmax:
        trait_cmin -= 0.5
        trait_cmax += 0.5

    if curtain_color_mode == "height":
        mesh_cmin = min(float(plot_data.trait_min), baseline)
        mesh_cmax = max(float(plot_data.trait_max), baseline)
        if mesh_cmin == mesh_cmax:
            mesh_cmin -= 0.5
            mesh_cmax += 0.5
    else:
        mesh_cmin, mesh_cmax = trait_cmin, trait_cmax

    mx, my, mz, mi, mj, mk, intensity = _build_branch_curtains_geometry(
        plot_data, baseline, curtain_color_mode
    )
    convert = _scene_transform(plot_data, baseline, trait_axis_scale)

    positions: List[float] = []
    for x, y, z in zip(mx, my, mz):
        positions.extend(convert(x, y, z))
    colors = _sample_colors(
        [float(v) for v in intensity],
        plot_data.colorscale,
        mesh_cmin,
        mesh_cmax,
        reverse_colorscale,
    )
    indices = [
        q for a, b, c in zip(mi, mj, mk) for q in (int(a), int(b), int(c))
    ]

    line_positions: List[float] = []
    line_traits: List[float] = []
    if show_centerline:
        for seg in plot_data.segments:
            for x, y, z, trait in (
                (seg.x0, seg.y0, seg.z0, seg.trait0),
                (seg.x1, seg.y1, seg.z1, seg.trait1),
            ):
                line_positions.extend(convert(x, y, z))
                line_traits.append(float(trait))

    if centerline_color == "trait":
        line_colors = _sample_colors(
            line_traits,
            plot_data.colorscale,
            trait_cmin,
            trait_cmax,
            reverse_colorscale,
        )
    else:
        rgb = [43 / 255.0] * 3 if centerline_color == "dark" else _css_color(centerline_color)
        line_colors = [rgb for _ in line_traits]

    labels = []
    if show_tip_labels:
        offset = (
            DEFAULT_TIP_LABEL_OFFSET_FRACTION
            if tip_label_offset is None
            else float(tip_label_offset)
        )
        if not math.isfinite(offset) or offset < 0:
            raise ValueError("tip_label_offset must be finite and >= 0")
        span = abs(float(plot_data.time_max) - float(plot_data.time_min))
        sign = -1.0 if plot_data.time_max >= plot_data.time_min else 1.0
        label_z = float(plot_data.time_min) + sign * offset * span
        for node in sorted(
            (n for n in plot_data.nodes.values() if n.is_tip),
            key=lambda n: n.x,
        ):
            labels.append(
                {"text": node.label, "position": convert(node.x, baseline, label_z)}
            )

    eye = CAMERA_PRESETS.get(
        camera_preset.lower(), CAMERA_PRESETS["elife"]
    )["eye"]

    gradient_values = [
        mesh_cmin + (mesh_cmax - mesh_cmin) * i / 10.0 for i in range(11)
    ]
    gradient = _sample_colors(
        gradient_values,
        plot_data.colorscale,
        mesh_cmin,
        mesh_cmax,
        reverse_colorscale,
    )

    nodes = []
    for node in plot_data.nodes.values():
        descendants = ""
        if not node.is_tip and node.descendant_tips:
            desc_sample = ", ".join(node.descendant_tips[:3])
            if len(node.descendant_tips) > 3:
                desc_sample += "..."
            descendants = f"Descendants ({len(node.descendant_tips)} tips): {desc_sample}"

        raw_val = float(node.raw_trait if node.raw_trait is not None else node.trait)
        disp_val = float(node.display_trait if node.display_trait is not None else node.y)

        nodes.append({
            "node_id": str(node.node_id),
            "label": str(node.label),
            "is_tip": bool(node.is_tip),
            "position": convert(node.x, node.y, node.z),
            "baseline_pos": convert(node.x, baseline, node.z),
            "x": float(node.x),
            "time": float(node.time),
            "raw_trait": raw_val,
            "display_trait": disp_val,
            "descendants": descendants,
        })

    y_tickvals, y_ticktext = _generate_rescaled_ticks(
        plot_data=plot_data,
        baseline_y=baseline,
        num_ticks=5,
        include_baseline=True,
    )
    if not y_tickvals:
        d_min = min(float(plot_data.trait_min), baseline)
        d_max = float(plot_data.trait_max)
        if d_max == d_min:
            y_tickvals = [d_min]
            y_ticktext = [f"{d_min:.4f}".rstrip("0").rstrip(".")]
        else:
            step = (d_max - d_min) / 4.0
            y_tickvals = [d_min + i * step for i in range(5)]
            y_ticktext = [
                str(int(round(v))) if abs(v - round(v)) < 1e-6 else f"{v:.4f}".rstrip("0").rstrip(".")
                for v in y_tickvals
            ]

    z_min = float(plot_data.time_min)
    z_max = float(plot_data.time_max)
    if z_max == z_min:
        z_tickvals = [z_min]
        z_ticktext = [f"{z_min:.4f}".rstrip("0").rstrip(".")]
    else:
        step_z = (z_max - z_min) / 4.0
        z_tickvals = [z_min + i * step_z for i in range(5)]
        z_ticktext = [
            str(int(round(v))) if abs(v - round(v)) < 1e-4 else f"{v:.1f}"
            for v in z_tickvals
        ]

    bx_min = float(plot_data.x_min)
    bx_max = float(plot_data.x_max)
    by_min = min(float(plot_data.trait_min), baseline)
    by_max = max([float(v) for v in y_tickvals] + [float(plot_data.trait_max)])
    bz_min = float(plot_data.time_min)
    bz_max = float(plot_data.time_max)

    # The four-layer renderer intentionally draws only the two scientific axes:
    # Y = Trait value and Z = Time before present. Their screen-facing corner is
    # selected dynamically in JavaScript from the current camera orientation.
    # X remains part of the scientific coordinate system and hover metadata but
    # has no visible axis/ticks in this backend.
    n_xmin = convert(bx_min, by_min, bz_min)[0]
    n_xmax = convert(bx_max, by_min, bz_min)[0]
    n_ymin = convert(bx_min, by_min, bz_min)[1]
    n_ymax = convert(bx_min, by_max, bz_min)[1]
    n_zmin = convert(bx_min, by_min, bz_min)[2]
    n_zmax = convert(bx_min, by_min, bz_max)[2]

    axes = {
        "show_x": bool(show_x_axis),
        "show_y": bool(show_y_axis),
        "show_z": bool(show_z_axis),
        "bounds": {
            "x_min": n_xmin,
            "x_max": n_xmax,
            "y_min": n_ymin,
            "y_max": n_ymax,
            "z_min": n_zmin,
            "z_max": n_zmax,
        },
        "y_axis": {
            "visible": bool(show_y_axis),
            "title": "Trait value",
            "ticks": [
                {
                    "val": float(v),
                    "text": str(t),
                    "coord": convert(bx_min, v, bz_min)[1],
                }
                for v, t in zip(y_tickvals, y_ticktext)
            ],
        },
        "z_axis": {
            "visible": bool(show_z_axis),
            "title": "Time before present",
            "ticks": [
                {
                    "val": float(v),
                    "text": str(t),
                    "coord": convert(bx_min, by_min, v)[2],
                }
                for v, t in zip(z_tickvals, z_ticktext)
            ],
        },
        "x_axis": {
            "visible": False,
            "show_numeric_labels": False,
        },
    }

    return {
        "schema": "phylo3d-four-layer-v1",
        "peel_layers": 4,
        "terminal_layer": 4,
        "opacity": opacity,
        "transparency": 1.0 - opacity,
        "max_omitted_transmittance": (1.0 - opacity) ** 4,
        "mesh": {
            "positions": positions,
            "colors": [v for rgb in colors for v in rgb],
            "indices": indices,
        },
        "centerline": {
            "positions": line_positions,
            "colors": [v for rgb in line_colors for v in rgb],
        },
        "labels": labels,
        "nodes": nodes,
        "axes": axes,
        "hover": {
            "show_tip": bool(show_tip_hover),
            "show_internal": bool(show_internal_hover),
        },
        "background": background,
        "camera_eye": [float(eye[k]) for k in ("x", "y", "z")],
        "title": plot_data.title,
        "colorbar": {
            "colors": gradient,
            "raw_min": float(plot_data.raw_trait_min),
            "raw_max": float(plot_data.raw_trait_max),
        },
        "stats": {
            "vertices": len(positions) // 3,
            "triangles": len(indices) // 3,
            "line_vertices": len(line_positions) // 3,
            "nodes": len(nodes),
        },
    }


_HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Phylo3D four-layer renderer</title>
<style>
html,body,#wrap{margin:0;width:100%;height:100%;overflow:hidden;font-family:Arial,sans-serif}
#wrap{position:relative;background:white}canvas{position:absolute;inset:0;width:100%;height:100%;display:block}
#axes{position:absolute;inset:0;width:100%;height:100%;pointer-events:none;overflow:visible}
.axis-line{stroke:#888;stroke-width:1.2;stroke-linecap:round}
.axis-tick{stroke:#888;stroke-width:1.2}
.axis-label{font-family:Arial,sans-serif;font-size:11px;fill:#444}
.axis-title{font-family:Arial,sans-serif;font-size:13px;font-weight:bold;fill:#222}
#labels{position:absolute;inset:0;pointer-events:none}.tip{position:absolute;font-size:12px;color:#222;white-space:nowrap}
#tooltip{display:none;position:absolute;pointer-events:none;background:rgba(255,255,255,.95);border:1px solid #ccc;border-radius:4px;padding:8px 10px;font:12px/1.4 Arial,sans-serif;color:#222;box-shadow:0 2px 8px rgba(0,0,0,.15);white-space:nowrap;z-index:100}
#title{position:absolute;left:50%;top:10px;transform:translateX(-50%);font-size:18px;color:#222}
#badge{position:absolute;left:10px;bottom:10px;background:rgba(255,255,255,.86);border:1px solid #ddd;padding:6px 8px;font:12px monospace}
#cb{position:absolute;right:16px;top:18%;height:64%;width:62px}#grad{position:absolute;right:0;top:18px;width:18px;height:calc(100% - 36px);border:1px solid #aaa}
#mx,#mn{position:absolute;right:24px;font-size:11px;color:#333}#mx{top:8px}#mn{bottom:8px}
#toolbar{position:absolute;right:14px;top:12px;display:flex;gap:4px;z-index:120}
#toolbar button{font:12px Arial,sans-serif;padding:5px 8px;border:1px solid #c9c9c9;border-radius:4px;background:rgba(255,255,255,.94);color:#333;cursor:pointer}
#toolbar button:hover{background:#f2f2f2}
#err{display:none;position:absolute;inset:20px;background:#fff4f4;border:1px solid #b00;padding:16px;color:#900;white-space:pre-wrap}
</style></head>
<body><div id="wrap"><canvas id="gl"></canvas><svg id="axes"></svg><div id="labels"></div><div id="tooltip"></div><div id="title"></div><div id="toolbar"><button id="reset-view" type="button">Reset</button><button id="download-png" type="button">PNG</button><button id="download-svg" type="button">SVG</button></div>
<div id="cb"><div id="grad"></div><div id="mx"></div><div id="mn"></div></div>
<div id="badge"></div><div id="err"></div></div>
<script>
"use strict";
const DATA=__PHYLO3D_PAYLOAD__;
const MAX_PEEL_LAYERS=4;
function fail(m){const e=document.getElementById("err");e.style.display="block";e.textContent=m;throw new Error(m);}
const canvas=document.getElementById("gl"),gl=canvas.getContext("webgl2",{alpha:true,antialias:true,premultipliedAlpha:false});
if(!gl)fail("WebGL2 is required.");
document.getElementById("title").textContent=DATA.title||"";
document.getElementById("mx").textContent=String(DATA.colorbar.raw_max);
document.getElementById("mn").textContent=String(DATA.colorbar.raw_min);
const stops=DATA.colorbar.colors.map((c,i)=>"rgb("+c.map(x=>Math.round(x*255)).join(",")+") "+(i*10)+"%");
document.getElementById("grad").style.background="linear-gradient(to top,"+stops.join(",")+")";
document.getElementById("badge").textContent="4-layer peel | opacity="+DATA.opacity.toFixed(2)+" | omitted <= "+(100*DATA.max_omitted_transmittance).toFixed(2)+"%";

function S(type,src){const s=gl.createShader(type);gl.shaderSource(s,src);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))fail(gl.getShaderInfoLog(s));return s}
function P(vs,fs){const p=gl.createProgram();gl.attachShader(p,S(gl.VERTEX_SHADER,vs));gl.attachShader(p,S(gl.FRAGMENT_SHADER,fs));gl.linkProgram(p);if(!gl.getProgramParameter(p,gl.LINK_STATUS))fail(gl.getProgramInfoLog(p));return p}

const VS=["#version 300 es","precision highp float;","layout(location=0) in vec3 aPosition;","layout(location=1) in vec3 aColor;","uniform mat4 uMVP;","out vec3 vColor;","void main(){vColor=aColor;gl_Position=uMVP*vec4(aPosition,1.0);}"].join("\n");
const PEEL=["#version 300 es","precision highp float;","in vec3 vColor;","uniform sampler2D uPrevDepth;","uniform vec2 uResolution;","uniform bool uHasPrev;","uniform float uOpacity;","out vec4 outColor;","void main(){if(uHasPrev){float p=texture(uPrevDepth,gl_FragCoord.xy/uResolution).r;if(gl_FragCoord.z<=p+0.0000005)discard;}outColor=vec4(vColor,uOpacity);}"].join("\n");
const QVS=["#version 300 es","precision highp float;","out vec2 vUV;","void main(){vec2 p=vec2((gl_VertexID==1||gl_VertexID==2||gl_VertexID==4)?1.0:-1.0,(gl_VertexID==2||gl_VertexID==4||gl_VertexID==5)?1.0:-1.0);vUV=p*.5+.5;gl_Position=vec4(p,0,1);}"].join("\n");
const COMP=["#version 300 es","precision highp float;","in vec2 vUV;","uniform sampler2D uLayer0,uLayer1,uLayer2,uLayer3;","uniform vec4 uBackground;","out vec4 outColor;","vec4 over(vec4 s,vec4 d){float a=s.a;return vec4(s.rgb*a+d.rgb*(1.0-a),a+d.a*(1.0-a));}","void main(){vec4 a=texture(uLayer0,vUV),b=texture(uLayer1,vUV),c=texture(uLayer2,vUV),d=texture(uLayer3,vUV),o=uBackground;if(d.a>0.0)o=vec4(d.rgb,1.0);if(c.a>0.0)o=over(c,o);if(b.a>0.0)o=over(b,o);if(a.a>0.0)o=over(a,o);outColor=o;}"].join("\n");
const LFS=["#version 300 es","precision highp float;","in vec3 vColor;","uniform sampler2D uFrontDepth;","uniform vec2 uResolution;","uniform bool uHasCurtain;","out vec4 outColor;","void main(){if(uHasCurtain){float d=texture(uFrontDepth,gl_FragCoord.xy/uResolution).r;if(d<0.999999&&gl_FragCoord.z>d+0.00002)discard;}outColor=vec4(vColor,1.0);}"].join("\n");
const peel=P(VS,PEEL),comp=P(QVS,COMP),line=P(VS,LFS);

function vao(pos,col,idx){const v=gl.createVertexArray();gl.bindVertexArray(v);let b=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(pos),gl.STATIC_DRAW);gl.enableVertexAttribArray(0);gl.vertexAttribPointer(0,3,gl.FLOAT,false,0,0);b=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(col),gl.STATIC_DRAW);gl.enableVertexAttribArray(1);gl.vertexAttribPointer(1,3,gl.FLOAT,false,0,0);if(idx){b=gl.createBuffer();gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,b);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,new Uint32Array(idx),gl.STATIC_DRAW)}gl.bindVertexArray(null);return v}
const meshV=vao(DATA.mesh.positions,DATA.mesh.colors,DATA.mesh.indices),lineV=vao(DATA.centerline.positions,DATA.centerline.colors,null),nidx=DATA.mesh.indices.length,nline=DATA.centerline.positions.length/3;

let targets=[];
function target(w,h){const fb=gl.createFramebuffer();gl.bindFramebuffer(gl.FRAMEBUFFER,fb);const color=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,color);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA8,w,h,0,gl.RGBA,gl.UNSIGNED_BYTE,null);gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.COLOR_ATTACHMENT0,gl.TEXTURE_2D,color,0);const depth=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,depth);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.texImage2D(gl.TEXTURE_2D,0,gl.DEPTH_COMPONENT24,w,h,0,gl.DEPTH_COMPONENT,gl.UNSIGNED_INT,null);gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.DEPTH_ATTACHMENT,gl.TEXTURE_2D,depth,0);if(gl.checkFramebufferStatus(gl.FRAMEBUFFER)!==gl.FRAMEBUFFER_COMPLETE)fail("Depth-peeling framebuffer incomplete");return{fb,color,depth}}
function resize(){const d=Math.min(window.devicePixelRatio||1,2),w=Math.max(1,Math.floor(canvas.clientWidth*d)),h=Math.max(1,Math.floor(canvas.clientHeight*d));if(canvas.width!==w||canvas.height!==h){for(const t of targets){gl.deleteFramebuffer(t.fb);gl.deleteTexture(t.color);gl.deleteTexture(t.depth)}targets=[];canvas.width=w;canvas.height=h;for(let i=0;i<4;i++)targets.push(target(w,h))}gl.viewport(0,0,w,h)}

function norm(v){const n=Math.hypot(v[0],v[1],v[2])||1;return[v[0]/n,v[1]/n,v[2]/n]}function cross(a,b){return[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]}function dot(a,b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]}
function look(eye){const z=norm(eye),x=norm(cross([0,1,0],z)),y=cross(z,x);return new Float32Array([x[0],y[0],z[0],0,x[1],y[1],z[1],0,x[2],y[2],z[2],0,-dot(x,eye),-dot(y,eye),-dot(z,eye),1])}
function ortho(l,r,b,t,n,f){const o=new Float32Array(16);o[0]=2/(r-l);o[5]=2/(t-b);o[10]=-2/(f-n);o[12]=-(r+l)/(r-l);o[13]=-(t+b)/(t-b);o[14]=-(f+n)/(f-n);o[15]=1;return o}
function mul(a,b){const o=new Float32Array(16);for(let c=0;c<4;c++)for(let r=0;r<4;r++){let s=0;for(let k=0;k<4;k++)s+=a[k*4+r]*b[c*4+k];o[c*4+r]=s}return o}
function project(m,p){const x=p[0],y=p[1],z=p[2],w=m[3]*x+m[7]*y+m[11]*z+m[15];return[(m[0]*x+m[4]*y+m[8]*z+m[12])/w,(m[1]*x+m[5]*y+m[9]*z+m[13])/w,(m[2]*x+m[6]*y+m[10]*z+m[14])/w]}

const iv=norm(DATA.camera_eye),initialYaw=Math.atan2(iv[0],iv[2]),initialPitch=Math.asin(Math.max(-.98,Math.min(.98,iv[1]))),initialZoom=2.2;let yaw=initialYaw,pitch=initialPitch,zoom=initialZoom,drag=false,lx=0,ly=0;
function mvp(){const d=4,cp=Math.cos(pitch),eye=[d*cp*Math.sin(yaw),d*Math.sin(pitch),d*cp*Math.cos(yaw)],a=canvas.width/canvas.height;return mul(ortho(-zoom*a,zoom*a,-zoom,zoom,-20,20),look(eye))}

const axesSvg=document.getElementById("axes");
const svgYAxis=document.createElementNS("http://www.w3.org/2000/svg","line");svgYAxis.setAttribute("class","axis-line");axesSvg.appendChild(svgYAxis);
const svgZAxis=document.createElementNS("http://www.w3.org/2000/svg","line");svgZAxis.setAttribute("class","axis-line");axesSvg.appendChild(svgZAxis);
const svgYTicks=(DATA.axes&&DATA.axes.y_axis?DATA.axes.y_axis.ticks:[]).map(t=>{const line=document.createElementNS("http://www.w3.org/2000/svg","line");line.setAttribute("class","axis-tick");axesSvg.appendChild(line);const text=document.createElementNS("http://www.w3.org/2000/svg","text");text.setAttribute("class","axis-label");text.textContent=t.text;axesSvg.appendChild(text);return{line,text,data:t}});
const svgYTitle=document.createElementNS("http://www.w3.org/2000/svg","text");svgYTitle.setAttribute("class","axis-title");if(DATA.axes&&DATA.axes.y_axis)svgYTitle.textContent=DATA.axes.y_axis.title;axesSvg.appendChild(svgYTitle);
const svgZTicks=(DATA.axes&&DATA.axes.z_axis?DATA.axes.z_axis.ticks:[]).map(t=>{const line=document.createElementNS("http://www.w3.org/2000/svg","line");line.setAttribute("class","axis-tick");axesSvg.appendChild(line);const text=document.createElementNS("http://www.w3.org/2000/svg","text");text.setAttribute("class","axis-label");text.textContent=t.text;axesSvg.appendChild(text);return{line,text,data:t}});
const svgZTitle=document.createElementNS("http://www.w3.org/2000/svg","text");svgZTitle.setAttribute("class","axis-title");if(DATA.axes&&DATA.axes.z_axis)svgZTitle.textContent=DATA.axes.z_axis.title;axesSvg.appendChild(svgZTitle);
const hoverGuide=document.createElementNS("http://www.w3.org/2000/svg","line");hoverGuide.setAttribute("id","hover-guide");hoverGuide.setAttribute("stroke","#e65100");hoverGuide.setAttribute("stroke-width","1.5");hoverGuide.setAttribute("stroke-dasharray","3,3");hoverGuide.style.display="none";axesSvg.appendChild(hoverGuide);
const hoverMarker=document.createElementNS("http://www.w3.org/2000/svg","circle");hoverMarker.setAttribute("id","hover-marker");hoverMarker.setAttribute("r","5");hoverMarker.setAttribute("fill","#ff7f0e");hoverMarker.setAttribute("stroke","#ffffff");hoverMarker.setAttribute("stroke-width","1.5");hoverMarker.style.display="none";axesSvg.appendChild(hoverMarker);

function toScreen(p){return[(p[0]*.5+.5)*canvas.clientWidth,(-p[1]*.5+.5)*canvas.clientHeight,p[2]]}
function setSvgLine(el,a,b){el.setAttribute("x1",a[0]);el.setAttribute("y1",a[1]);el.setAttribute("x2",b[0]);el.setAttribute("y2",b[1])}
function frontAxisCorner(){
const b=DATA.axes.bounds,cp=Math.cos(pitch),eyeX=cp*Math.sin(yaw),eyeZ=cp*Math.cos(yaw);
return{x:eyeX>=0?b.x_max:b.x_min,z:eyeZ>=0?b.z_max:b.z_min}
}
function outward2D(M,corner){
const sc=toScreen(project(M,[0,0,0])),sa=toScreen(project(M,[corner.x,DATA.axes.bounds.y_min,corner.z]));
let dx=sa[0]-sc[0],dy=sa[1]-sc[1],n=Math.hypot(dx,dy);
if(n<1e-6){dx=1;dy=0;n=1}
return[dx/n,dy/n]
}
function placeAxisText(el,base,out,offset){
el.setAttribute("x",base[0]+out[0]*offset);el.setAttribute("y",base[1]+out[1]*offset);
el.setAttribute("text-anchor",out[0]<-.2?"end":out[0]>.2?"start":"middle");
el.setAttribute("dominant-baseline",out[1]<-.2?"auto":out[1]>.2?"hanging":"central")
}
function updateAxes(M){
if(!DATA.axes)return;
const b=DATA.axes.bounds,corner=frontAxisCorner(),out=outward2D(M,corner);
const showY=DATA.axes.show_y!==false,showZ=DATA.axes.show_z!==false;
const y0=toScreen(project(M,[corner.x,b.y_min,corner.z])),y1=toScreen(project(M,[corner.x,b.y_max,corner.z]));
if(showY){svgYAxis.style.display="block";setSvgLine(svgYAxis,y0,y1)}else svgYAxis.style.display="none";
const z0=toScreen(project(M,[corner.x,b.y_min,b.z_min])),z1=toScreen(project(M,[corner.x,b.y_min,b.z_max]));
if(showZ){svgZAxis.style.display="block";setSvgLine(svgZAxis,z0,z1)}else svgZAxis.style.display="none";

svgYTicks.forEach(item=>{
if(!showY){item.line.style.display="none";item.text.style.display="none";return}
const s=toScreen(project(M,[corner.x,item.data.coord,corner.z])),e=[s[0]+out[0]*7,s[1]+out[1]*7];
item.line.style.display="block";item.text.style.display="block";setSvgLine(item.line,s,e);placeAxisText(item.text,e,out,5)
});
if(showY){svgYTitle.style.display="block";const mid=toScreen(project(M,[corner.x,(b.y_min+b.y_max)/2,corner.z]));placeAxisText(svgYTitle,mid,out,48)}else svgYTitle.style.display="none";

svgZTicks.forEach(item=>{
if(!showZ){item.line.style.display="none";item.text.style.display="none";return}
const s=toScreen(project(M,[corner.x,b.y_min,item.data.coord])),e=[s[0]+out[0]*7,s[1]+out[1]*7];
item.line.style.display="block";item.text.style.display="block";setSvgLine(item.line,s,e);placeAxisText(item.text,e,out,5)
});
if(showZ){svgZTitle.style.display="block";const mid=toScreen(project(M,[corner.x,b.y_min,(b.z_min+b.z_max)/2]));placeAxisText(svgZTitle,mid,out,34)}else svgZTitle.style.display="none";
}

const lr=document.getElementById("labels"),le=DATA.labels.map(x=>{const d=document.createElement("div");d.className="tip";d.textContent=x.text;lr.appendChild(d);return d});
function labels(M){const eyeX=Math.cos(pitch)*Math.sin(yaw);const tx=eyeX<0?"translate(calc(-100% - 3px),-50%)":"translate(3px,-50%)";DATA.labels.forEach((x,i)=>{const p=project(M,x.position),e=le[i];e.style.display=(p[2]<-1||p[2]>1)?"none":"block";e.style.left=((p[0]*.5+.5)*canvas.clientWidth)+"px";e.style.top=((-p[1]*.5+.5)*canvas.clientHeight)+"px";e.style.transform=tx})}

let currentM=null,projectedNodes=[],lastMouse=null;
const tooltip=document.getElementById("tooltip");
function updateProjectedNodes(M){
currentM=M;
if(!DATA.nodes)return;
projectedNodes=DATA.nodes.map(n=>{const p=project(M,n.position);return{node:n,visible:(p[2]>=-1&&p[2]<=1),sx:(p[0]*.5+.5)*canvas.clientWidth,sy:(-p[1]*.5+.5)*canvas.clientHeight,depth:p[2]}});
}
function updateHover(mx,my){
lastMouse={x:mx,y:my};
if(!projectedNodes.length)return;
const showTip = !DATA.hover || DATA.hover.show_tip !== false;
const showInternal = !DATA.hover || DATA.hover.show_internal !== false;
if(!showTip && !showInternal){
tooltip.style.display="none";
hoverMarker.style.display="none";
hoverGuide.style.display="none";
return;
}
const radius=18;
let closest=null,minDist=radius;
for(const p of projectedNodes){
if(!p.visible)continue;
if(p.node.is_tip && !showTip)continue;
if(!p.node.is_tip && !showInternal)continue;
const d=Math.hypot(p.sx-mx,p.sy-my);
if(d<minDist){minDist=d;closest=p}
}
if(closest){
const n=closest.node;
if(n.is_tip){
tooltip.innerHTML="<b>Taxon: "+n.label+"</b><br>Node ID: "+n.node_id+"<br>Trait value (Y / Height): "+n.raw_trait.toFixed(4)+"<br>Time before present (Z): "+n.time.toFixed(4)+"<br>Tree Layout (X): "+n.x.toFixed(2);
}else{
tooltip.innerHTML="<b>Node: "+n.node_id+"</b><br>Type: Ancestral Node<br>Trait value (Y / Height): "+n.raw_trait.toFixed(4)+"<br>Time before present (Z): "+n.time.toFixed(4)+"<br>Tree Layout (X): "+n.x.toFixed(2)+(n.descendants?("<br>"+n.descendants):"");
}
tooltip.style.display="block";
const tw=tooltip.offsetWidth||180,th=tooltip.offsetHeight||100;
const left=Math.max(10,Math.min(window.innerWidth-tw-10,mx+12));
const top=Math.max(10,Math.min(window.innerHeight-th-10,my+12));
tooltip.style.left=left+"px";
tooltip.style.top=top+"px";
hoverMarker.setAttribute("cx",closest.sx);
hoverMarker.setAttribute("cy",closest.sy);
hoverMarker.style.display="block";
if(n.baseline_pos&&currentM){
const bp=project(currentM,n.baseline_pos);
const bx=(bp[0]*.5+.5)*canvas.clientWidth,by=(-bp[1]*.5+.5)*canvas.clientHeight;
hoverGuide.setAttribute("x1",closest.sx);
hoverGuide.setAttribute("y1",closest.sy);
hoverGuide.setAttribute("x2",bx);
hoverGuide.setAttribute("y2",by);
hoverGuide.style.display="block";
}else{hoverGuide.style.display="none"}
}else{
tooltip.style.display="none";
hoverMarker.style.display="none";
hoverGuide.style.display="none";
}
}

function escXml(v){return String(v).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;")}
function currentExportSvg(){
render();
const w=canvas.clientWidth,h=canvas.clientHeight,raster=canvas.toDataURL("image/png");
const axisMarkup=axesSvg.innerHTML;
let labelMarkup="";
for(const el of le){if(el.style.display==="none")continue;const r=el.getBoundingClientRect(),root=document.getElementById("wrap").getBoundingClientRect();labelMarkup+='<text x="'+(r.left-root.left).toFixed(2)+'" y="'+(r.top-root.top+12).toFixed(2)+'" font-family="Arial,sans-serif" font-size="12" fill="#222">'+escXml(el.textContent)+'</text>'}
const colors=DATA.colorbar.colors,stopsSvg=colors.map((c,i)=>'<stop offset="'+(i/(colors.length-1)*100).toFixed(2)+'%" stop-color="rgb('+c.map(x=>Math.round(x*255)).join(",")+')" />').join("");
const cbX=w-32,cbY=h*.18,cbH=h*.64;
const title=DATA.title?'<text x="'+(w/2)+'" y="28" text-anchor="middle" font-family="Arial,sans-serif" font-size="18" fill="#222">'+escXml(DATA.title)+'</text>':"";
const colorbar='<defs><linearGradient id="cbgrad" x1="0" y1="1" x2="0" y2="0">'+stopsSvg+'</linearGradient></defs><rect x="'+cbX+'" y="'+cbY+'" width="18" height="'+cbH+'" fill="url(#cbgrad)" stroke="#aaa"/><text x="'+(cbX-6)+'" y="'+(cbY+10)+'" text-anchor="end" font-family="Arial,sans-serif" font-size="11" fill="#333">'+escXml(DATA.colorbar.raw_max)+'</text><text x="'+(cbX-6)+'" y="'+(cbY+cbH)+'" text-anchor="end" font-family="Arial,sans-serif" font-size="11" fill="#333">'+escXml(DATA.colorbar.raw_min)+'</text>';
return '<svg xmlns="http://www.w3.org/2000/svg" width="'+w+'" height="'+h+'" viewBox="0 0 '+w+' '+h+'"><rect width="100%" height="100%" fill="'+(DATA.background==="transparent"?"none":"white")+'"/><image href="'+raster+'" x="0" y="0" width="'+w+'" height="'+h+'"/><g>'+axisMarkup+'</g>'+labelMarkup+title+colorbar+'</svg>'
}
function downloadBlob(blob,name){const a=document.createElement("a");a.href=URL.createObjectURL(blob);a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(a.href),1000)}
function exportSvg(){downloadBlob(new Blob([currentExportSvg()],{type:"image/svg+xml;charset=utf-8"}),"phylo3d_trait.svg")}
function exportPng(){
const svg=currentExportSvg(),blob=new Blob([svg],{type:"image/svg+xml;charset=utf-8"}),url=URL.createObjectURL(blob),img=new Image();
img.onload=()=>{const scale=2,c=document.createElement("canvas");c.width=Math.round(canvas.clientWidth*scale);c.height=Math.round(canvas.clientHeight*scale);const ctx=c.getContext("2d");ctx.scale(scale,scale);ctx.drawImage(img,0,0,canvas.clientWidth,canvas.clientHeight);URL.revokeObjectURL(url);c.toBlob(b=>downloadBlob(b,"phylo3d_trait.png"),"image/png")};
img.onerror=()=>{URL.revokeObjectURL(url);fail("PNG export failed")};img.src=url
}
document.getElementById("reset-view").addEventListener("click",()=>{yaw=initialYaw;pitch=initialPitch;zoom=initialZoom;render()});
document.getElementById("download-svg").addEventListener("click",exportSvg);
document.getElementById("download-png").addEventListener("click",exportPng);

canvas.addEventListener("pointerdown",e=>{drag=true;lx=e.clientX;ly=e.clientY;canvas.setPointerCapture(e.pointerId);tooltip.style.display="none";hoverMarker.style.display="none";hoverGuide.style.display="none"});
canvas.addEventListener("pointerup",()=>drag=false);
canvas.addEventListener("pointermove",e=>{if(drag){tooltip.style.display="none";hoverMarker.style.display="none";hoverGuide.style.display="none";yaw-=(e.clientX-lx)*.008;pitch+=(e.clientY-ly)*.008;pitch=Math.max(-1.45,Math.min(1.45,pitch));lx=e.clientX;ly=e.clientY;render();return}updateHover(e.clientX,e.clientY)});
canvas.addEventListener("wheel",e=>{e.preventDefault();zoom*=Math.exp(e.deltaY*.001);zoom=Math.max(.5,Math.min(8,zoom));render()},{passive:false});
canvas.addEventListener("pointerleave",()=>{tooltip.style.display="none";hoverMarker.style.display="none";hoverGuide.style.display="none"});

function render(){resize();const M=mvp();gl.disable(gl.BLEND);gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LESS);gl.depthMask(true);gl.disable(gl.CULL_FACE);gl.useProgram(peel);gl.bindVertexArray(meshV);gl.uniformMatrix4fv(gl.getUniformLocation(peel,"uMVP"),false,M);gl.uniform2f(gl.getUniformLocation(peel,"uResolution"),canvas.width,canvas.height);gl.uniform1f(gl.getUniformLocation(peel,"uOpacity"),DATA.opacity);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);
for(let pass=0;pass<4;pass++){gl.bindFramebuffer(gl.FRAMEBUFFER,targets[pass].fb);gl.clearColor(0,0,0,0);gl.clearDepth(1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.uniform1i(gl.getUniformLocation(peel,"uHasPrev"),pass>0?1:0);if(pass>0){gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,targets[pass-1].depth);gl.uniform1i(gl.getUniformLocation(peel,"uPrevDepth"),0)}if(nidx)gl.drawElements(gl.TRIANGLES,nidx,gl.UNSIGNED_INT,0)}
gl.bindFramebuffer(gl.FRAMEBUFFER,null);gl.disable(gl.DEPTH_TEST);gl.useProgram(comp);for(let i=0;i<4;i++){gl.activeTexture(gl.TEXTURE0+i);gl.bindTexture(gl.TEXTURE_2D,targets[i].color);gl.uniform1i(gl.getUniformLocation(comp,"uLayer"+i),i)}const bg=gl.getUniformLocation(comp,"uBackground");if(DATA.background==="transparent")gl.uniform4f(bg,0,0,0,0);else gl.uniform4f(bg,1,1,1,1);gl.drawArrays(gl.TRIANGLES,0,6);for(let i=0;i<4;i++){gl.activeTexture(gl.TEXTURE0+i);gl.bindTexture(gl.TEXTURE_2D,null)}gl.activeTexture(gl.TEXTURE0);
if(nline){gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);gl.useProgram(line);gl.bindVertexArray(lineV);gl.uniformMatrix4fv(gl.getUniformLocation(line,"uMVP"),false,M);gl.uniform2f(gl.getUniformLocation(line,"uResolution"),canvas.width,canvas.height);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,targets[0].depth);gl.uniform1i(gl.getUniformLocation(line,"uFrontDepth"),0);gl.uniform1i(gl.getUniformLocation(line,"uHasCurtain"),nidx?1:0);gl.drawArrays(gl.LINES,0,nline);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);gl.disable(gl.BLEND)}gl.bindVertexArray(null);gl.bindFramebuffer(gl.FRAMEBUFFER,null);
labels(M);
updateAxes(M);
updateProjectedNodes(M);
if(!drag&&lastMouse)updateHover(lastMouse.x,lastMouse.y);
}
window.addEventListener("resize",render);render();
</script></body></html>"""


def write_four_layer_html(
    plot_data: PlotData,
    output_path: Union[Path, str],
    *,
    opacity: float = 0.9,
    baseline_y: Optional[float] = None,
    reverse_colorscale: bool = False,
    curtain_color_mode: str = "height",
    trait_axis_scale: float = 1.0,
    tip_label_offset: Optional[float] = None,
    show_tip_labels: bool = True,
    show_centerline: bool = True,
    centerline_color: str = "dark",
    background: str = "white",
    camera_preset: str = "elife",
    show_x_axis: bool = True,
    show_y_axis: bool = True,
    show_z_axis: bool = True,
    show_tip_hover: bool = True,
    show_internal_hover: bool = True,
) -> Dict[str, Any]:
    """Write an interactive four-layer WebGL2 HTML visualization."""
    payload = _build_payload(
        plot_data,
        opacity=opacity,
        baseline_y=baseline_y,
        reverse_colorscale=reverse_colorscale,
        curtain_color_mode=curtain_color_mode,
        trait_axis_scale=trait_axis_scale,
        tip_label_offset=tip_label_offset,
        show_tip_labels=show_tip_labels,
        show_centerline=show_centerline,
        centerline_color=centerline_color,
        background=background,
        camera_preset=camera_preset,
        show_x_axis=show_x_axis,
        show_y_axis=show_y_axis,
        show_z_axis=show_z_axis,
        show_tip_hover=show_tip_hover,
        show_internal_hover=show_internal_hover,
    )
    html = _HTML.replace(
        "__PHYLO3D_PAYLOAD__",
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
    )
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return payload
