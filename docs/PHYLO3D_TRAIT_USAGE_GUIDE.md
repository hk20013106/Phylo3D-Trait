# Phylo3D-Trait 软件使用说明书
## 面向 AI Agent 与科研用户的标准操作手册

> **适用对象**：Antigravity、Codex、Claude Code、其他 AI Agent，以及人工使用者  
> **软件目标**：将一棵带 branch length 的系统发育树，以及所有 tip / internal node / root 的连续 Trait 数值，转换为一个可交互旋转的三维 rectangular phylogram，并将 Trait 同时映射为三维高度与连续颜色。  
> **重要原则**：本软件负责**可视化**，不负责 ancestral-state reconstruction（ASR）或 Trait 推断。

---

# 1. 软件做什么

`Phylo3D-Trait` 是一个 Python package，用于生成类似 eLife Figure 5 风格的三维系统发育可视化。

输入：

1. 一棵系统发育树：
   - Newick
   - Nexus
2. 所有节点的 Trait 值：
   - 所有 terminal tips
   - 所有 internal nodes
   - root / MRCA

输出：

- 一个 Plotly/WebGL 交互式 HTML 文件
- 支持鼠标自由旋转、缩放和平移
- 默认使用正交投影，避免透视导致的“近大远小”

核心科学坐标定义：

```text
X = Tree layout
Y = Trait value
Z = Time before present
```

因此核心不变量为：

```text
Y = Trait = Surface height = Surface color
```

Trait 越高，曲面在 Y 轴越高；相同 Trait 在整棵树上必须使用相同颜色。

---

# 2. AI Agent 使用本软件前的强制顺序

任何 AI Agent 处理真实数据前，必须：

1. **完整阅读本说明书。**
2. 检查当前仓库结构、`README.md`、`pyproject.toml`、`phylo3d_trait/` 和 CLI help。
3. 检查输入 tree 与 Trait table。
4. 确认 branch length 的含义。
5. 先生成 internal-node 模板，不人工猜 internal-node ID。
6. 验证所有 tip、internal node 和 root 均有 Trait。
7. 再运行 `plot`。
8. 验证 HTML 实际生成。
9. 仅当现有输入契约无法处理数据或确认存在 bug 时，才修改 Python 源码。

**默认原则：换数据 ≠ 改代码。**

---

# 3. 项目结构与模块职责

推荐结构：

```text
Phylo3D-Trait/
│
├── README.md
├── AGENTS.md
├── pyproject.toml
├── .gitignore
│
├── phylo3d_trait/
│   ├── __init__.py
│   ├── cli.py
│   ├── io.py
│   ├── models.py
│   ├── renderer.py
│   ├── template.py
│   └── tree.py
│
├── docs/
│   └── PHYLO3D_TRAIT_USAGE_GUIDE.md
│
├── examples/
│   ├── example1/
│   │   ├── tree.nwk
│   │   ├── node_values.csv
│   │   └── README.md
│   │
│   └── example2/
│       ├── tree.nwk
│       ├── node_values.csv
│       └── README.md
│
├── data/
│   └── .gitkeep
│
├── results/
│   ├── example1/
│   │   └── tree3d.html
│   └── example2/
│       └── tree3d.html
│
└── tests/
    ├── test_cli.py
    ├── test_tree.py
    ├── test_renderer.py
    ├── test_curtain_mesh.py
    ├── test_rectangular_phylogram.py
    └── test_stable_ids.py
```

模块职责：

- `cli.py`：正式命令行入口，提供 `template-values` 与 `plot`
- `tree.py`：树解析、时间映射、stable clade ID、rectangular geometry
- `models.py`：内部数据模型，如 `AnnotatedNode`、`EdgeSegment`、`PlotData`
- `template.py`：生成所有 tip/internal node 的 Trait 模板
- `io.py`：读取 CSV/TSV Trait 表
- `renderer.py`：Plotly Mesh3d、颜色、背景、相机、遮挡、标签

---

# 4. Rectangular 3D phylogram 几何

本软件不是 slanted cladogram。

对 biological edge：

```text
parent -> child
```

设：

```text
parent = (Xp, Yp, Zp)
child  = (Xc, Yc, Zc)
```

拆成两段：

## 4.1 Connector segment

```text
(Xp, Yp, Zp)
    ->
(Xc, Yp, Zp)
```

含义：

- `X` 变化
- `Y` 不变
- `Z` 不变
- 仅表示 rectangular tree 的布局连接
- 不代表时间演化

所以：

```text
Trait = Yp
Color = Yp
```

沿 connector 保持不变。

## 4.2 Lineage segment

```text
(Xc, Yp, Zp)
    ->
(Xc, Yc, Zc)
```

含义：

- `X` 不变
- `Y` 从 parent Trait 变化到 child Trait
- `Z` 从 parent time 变化到 child time

这才是 lineage-through-time 段。

---

# 5. Trait 与颜色

整个系统必须满足：

```text
Y == Trait == ColorValue
```

对于 Mesh3d 任意 vertex：

```text
intensity == vertex.y
```

例如：

```text
Y=1 -> Trait 1 color
Y=2 -> Trait 2 color
Y=3 -> Trait 3 color
Y=4 -> Trait 4 color
Y=5 -> Trait 5 color
```

颜色范围必须按整个 tree 的全局范围归一化，不能每条 branch 单独归一化。

---

# 6. Baseline

默认：

```text
baseline_y = trait_min
```

因此：

```text
minimum height = minimum Trait = minimum color value
```

如果用户显式指定更低 baseline，则 color range 也必须覆盖 baseline，避免颜色被错误裁剪。

---

# 7. 遮挡

默认 publication 模式：

```text
opacity = 1.0
```

所有 surface 完全不透明，由 WebGL depth buffer 根据当前 camera 自动决定前后遮挡。

不要：

- 根据 X 顺序手动隐藏后面 branch
- 根据 tip order 决定前后
- 用透明 surface 模拟 publication figure

---

# 8. Camera

坐标：

```text
X = Tree layout
Y = Trait
Z = Time before present
```

为了让 Trait 轴在屏幕上接近竖直，推荐：

```python
up = dict(x=0, y=1, z=0)
```

root/MRCA：

```text
Z = maximum time before present
```

tips：

```text
Z = 0
```

eLife-style 初始视角建议从 `+Z` 一侧观察，使 MRCA 更接近前景、tips 向后展开。

初始 camera 只是默认视角，HTML 打开后仍可自由旋转。

---

# 9. 背景与节点

推荐默认：

```text
background = white
showbackground = False
internal-node markers = off
tip markers = off
tip labels = on
centerline = thin
projection = orthographic
```

背景也可设为透明，但不建议浅灰色大面积 wall。

---

# 10. 安装

在仓库根目录：

```bash
pip install -e .
```

开发/测试环境：

```bash
pip install -e .[dev]
```

检查 CLI：

```bash
python -m phylo3d_trait.cli --help
```

如果已安装 entry point，也可：

```bash
phylo3d-trait --help
```

---

# 11. 真实数据推荐目录

```text
data/
  my_project/
    tree.nwk
    node_values_template.csv
    node_values.csv
    README.md

results/
  my_project/
    tree3d.html
```

这些路径只是推荐，不是硬编码要求。

---

# 12. 输入 1：系统发育树

支持：

```text
Newick
Nexus
```

常见扩展名：

```text
.nwk
.newick
.tre
.tree
.nex
.nexus
```

示例：

```newick
((A:10,B:10):20,(C:15,D:15):15);
```

---

# 13. Branch length 的含义必须先确认

这是最重要的科学检查之一。

如果 branch length 是时间，例如 Ma，那么：

```text
Z = Time before present
```

可以直接做时间解释。

如果 branch length 是：

```text
substitutions/site
```

则不能把 Z 轴标为：

```text
Time before present (Ma)
```

AI 如果不能确认 branch length 的单位，应停止并向用户确认，不能自行猜测。

---

# 14. 输入 2：Trait 表

标准格式至少包含：

```csv
node_id,trait
Species_A,1.25
Species_B,2.10
clade:xxxxxxxxxxxx,1.64
```

Tip：

```text
node_id = tree 中 exact tip name
```

不要自行修改大小写、空格、下划线或物种名格式。

Internal node：

```text
node_id = 程序生成的 stable clade ID
```

不要使用不稳定的 `node 51`、`node 52` 等临时编号。

---

# 15. Stable internal-node ID

Internal node identity 基于 descendant-tip set 的确定性哈希。

因此：

```text
(A,B)
```

和：

```text
(B,A)
```

必须识别为同一 clade。

该设计避免：

- reroot
- ladderize
- child-order rotation
- 不同软件 read/write

导致内部节点临时编号改变。

**不要人工猜 internal-node ID。**

---

# 16. 标准工作流：Step 1 生成节点模板

运行：

```bash
python -m phylo3d_trait.cli template-values \
  --tree data/my_project/tree.nwk \
  --output data/my_project/node_values_template.csv
```

生成表结构如下（由 `template-values` 自动生成）：

```csv
node_id,label,node_type,descendant_count,descendant_tips,trait
Species_A,Species_A,tip,1,Species_A,
Species_B,Species_B,tip,1,Species_B,
clade:b17c8419f544,clade:b17c8419f544,internal,2,Species_A;Species_B,
clade:17f5f129f4c7,clade:17f5f129f4c7,root,4,Species_A;Species_B;Species_C;Species_D,
```

其中 `trait` 列暂为空，供用户填入测量值与 ASR 重建值。

> [!NOTE]
> `Phylo3D-Trait` 读取表时仅需 `node_id` 与 `trait`（或 `value`）两列，模板中附带的辅助列（`label`、`node_type`、`descendant_count`、`descendant_tips`）仅用于方便人工核对后代构成，保留在表中不会影响程序正常运行。

---

# 17. Step 2 填入 Trait

用户自行提供：

- 所有 tips
- 所有 internal nodes
- root/MRCA

的 Trait。

示例：

```csv
node_id,label,node_type,descendant_count,descendant_tips,trait
Species_A,Species_A,tip,1,Species_A,1.5
Species_B,Species_B,tip,1,Species_B,3.0
Species_C,Species_C,tip,1,Species_C,4.5
Species_D,Species_D,tip,1,Species_D,5.0
clade:b17c8419f544,clade:b17c8419f544,internal,2,Species_A;Species_B,2.2
clade:6a5756530335,clade:6a5756530335,internal,2,Species_C;Species_D,4.0
clade:17f5f129f4c7,clade:17f5f129f4c7,root,4,Species_A;Species_B;Species_C;Species_D,1.8
```

保存为：

```text
data/my_project/node_values.csv
```

---

# 18. 软件不做 ASR

`Phylo3D-Trait` **不执行 ancestral-state reconstruction**。

它不会：

- 根据 tip 自动推断 ancestor
- 自动填补 missing internal nodes
- 默认 Brownian motion
- 自动做 ML/Bayesian ASR
- 自动从 sequence 算 Trait

缺少任何必要节点 Trait 时应报错停止，而不是静默填 0、复制 parent 或自行插值。

---

# 19. Step 3 生成 3D HTML

运行：

```bash
python -m phylo3d_trait.cli plot \
  --tree data/my_project/tree.nwk \
  --values data/my_project/node_values.csv \
  --output results/my_project/tree3d.html
```

然后打开：

```text
results/my_project/tree3d.html
```

---

# 20. AI 的标准真实数据流程

```text
READ README
↓
INSPECT TREE
↓
VERIFY BRANCH-LENGTH MEANING
↓
GENERATE TEMPLATE
↓
MAP ALL TIP + INTERNAL + ROOT TRAITS
↓
VALIDATE
↓
RUN PLOT
↓
VERIFY HTML
↓
REPORT
```

---

# 21. AI 运行前检查清单

## Tree

- [ ] 文件存在
- [ ] 可解析
- [ ] tip names 唯一
- [ ] branch length 存在
- [ ] branch length 单位已知
- [ ] rooting 已确认
- [ ] 若 Z 被解释为 Ma，确认 tree 确实 dated

## Trait

- [ ] Trait 文件存在
- [ ] 包含 `node_id`
- [ ] 包含 `trait`
- [ ] trait 为 numeric
- [ ] 每个 tip 有值
- [ ] 每个 internal node 有值
- [ ] root 有值
- [ ] 无重复 node_id
- [ ] internal-node ID 与该 tree 匹配

## Output

- [ ] output 路径有效
- [ ] HTML 成功生成
- [ ] 文件大小 > 0
- [ ] 浏览器可打开
- [ ] topology 合理
- [ ] Trait 高低与 Y 高度一致
- [ ] Trait 高低与颜色一致

---

# 22. AI 收到简单任务时应该如何操作

例如用户说：

> “用 `data/project1/` 里的数据画三维树。”

AI 应：

1. 阅读本说明书。
2. 检查 `data/project1/`。
3. 找到 tree 和 Trait 文件。
4. 检查 Trait 表是否已经使用当前 tree 的 stable clade IDs。
5. 如果没有，先运行 `template-values`。
6. 将用户提供的 ancestor values 映射到模板。
7. 验证完整性。
8. 运行 `plot`。
9. 检查 HTML。
10. 报告输出路径、实际命令和警告。

默认不要修改 Python package。

---

# 23. Python API

CLI 是推荐入口。

若需在 Python workflow 内调用，可按当前 package API 使用，例如：

```python
from phylo3d_trait import (
    parse_tree,
    load_trait_values,
    build_plot_data,
    build_figure,
)

tree = parse_tree("data/my_project/tree.nwk")
traits = load_trait_values("data/my_project/node_values.csv")

plot_data = build_plot_data(tree, traits)

fig = build_figure(
    plot_data,
    camera_preset="elife",
    background="white",
)

fig.write_html("results/my_project/tree3d.html")
```

**注意**：函数名必须以当前 `__init__.py` 和实际代码为准。AI 不应凭记忆假设 API。

---

# 24. CLI Help 是当前接口的权威来源

如果文档与当前版本可能不一致，运行：

```bash
python -m phylo3d_trait.cli --help
python -m phylo3d_trait.cli template-values --help
python -m phylo3d_trait.cli plot --help
```

优先服从当前 CLI 实际支持的参数。

---

# 25. 常用显示选项

正式 `plot` 命令支持的核心参数：

| 参数 | 类型 / 可选值 | 默认值 | 作用说明 |
|---|---|---|---|
| `--colorscale` | string | `Turbo` | 连续色标（如 `Turbo`、`Viridis`、`Plasma`、`Spectral`） |
| `--camera-preset` | `elife`, `root_front`, `tips_front` | `elife` | 初始视角（`elife`: MRCA 在前景、Trait 轴竖直、正交投影） |
| `--background` | `white`, `transparent` | `white` | 背景风格（默认纯白，关闭 3D 墙壁；或全透明） |
| `--opacity` | float (0.0 - 1.0) | `1.0` | 幕帘曲面不透明度（默认 1.0 启用原生 WebGL 深度遮挡） |
| `--branch-width` | float | `1.0` | 分支顶缘轮廓线线宽 |
| `--baseline-y` | float | `trait_min` | 自定义基准平面 $Y$ 高度 |
| `--trait-display-range` | float float (例如 `13 5`) | `None` | 线性显示重标定（将原始 Trait `[min, max]` 线性映射至 `[START, END]` 显示空间；不修改科研原始数据，hover 保留 raw 与 display 两套数值） |
| `--segments, -s` | int | `10` | 每一个分支段内的线性插值步数 |
| `--show-node-markers` | flag | `False` | 是否在内部祖先节点处渲染菱形 marker |
| `--no-mesh` | flag | `False` | 是否关闭连续垂直幕帘曲面 |
| `--no-centerline` | flag | `False` | 是否关闭分支顶缘轮廓线 |
| `--centerline-color` | `dark`, `trait`, 或 CSS color | `dark` | 顶缘轮廓线颜色模式 |
| `--no-labels` | flag | `False` | 是否隐藏末端物种名称文字标签 |
| `--title` | string | 默认标题 | 图像顶部标题文本 |

可以通过 `python -m phylo3d_trait.cli plot --help` 查看所有参数的实时官方说明。

---

# 26. 推荐 publication 默认值

```text
tree geometry = rectangular
mesh = on
opacity = 1.0
background = white
camera preset = elife
internal-node markers = off
tip markers = off
tip labels = on
centerline = thin
projection = orthographic
baseline = trait_min
```

---

# 27. 大树注意事项

当 tips 很多时，可能出现：

- Mesh 数量增加
- WebGL 压力上升
- tip labels 重叠
- branch overlap 增多
- HTML 文件变大

优先处理顺序：

1. 不改科学几何。
2. 关闭 tip labels。
3. 降低 edge sampling density。
4. 关闭非必要 centerline。
5. 保持 opacity=1.0 的真实遮挡。
6. 先 benchmark，再修改 renderer。

不要通过删除 branch 来解决拥挤。

---

# 28. 非 ultrametric tree

如果 tree 非 ultrametric，tips 可能不会自然位于：

```text
Z = 0
```

此时不能无条件把 Z 解释为 time before present。

AI 应先检查：

- tree 是否 dated
- tips 是否 contemporaneous
- branch length 是否时间

不要静默改时间轴。

---

# 29. Polytomy

如果 tree 有 polytomy，不要未经允许自动 resolve。

自动 resolve 会改变 topology。

---

# 30. Rooting

如果 tree 未定根，而图需要解释 MRCA / time before present：

```text
STOP
```

不要未经用户同意自行 midpoint-root。

---

# 31. Missing Trait

缺失任何：

- tip Trait
- internal-node Trait
- root Trait

默认应：

```text
FAIL LOUDLY
```

不要自动填补。

---

# 32. Trait 元数据

每个真实项目最好记录：

```text
Trait name:
Biological meaning:
Unit:
Tip calculation:
Ancestral reconstruction method:
Reference:
```

这对于科研可重复性很重要。

---

# 33. 每个真实数据项目建议自带 README

例如：

```text
data/my_project/README.md
```

推荐内容：

```markdown
# Dataset description

## Tree
- source:
- rooting:
- number of tips:
- branch-length unit:
- dated/undated:

## Trait
- name:
- unit:
- tip value source:
- ancestral reconstruction method:

## Files
- tree.nwk
- node_values.csv

## Notes
- ...
```

这个 project-specific README 是数据说明，不是软件配置文件。

---

# 34. README 不是机器输入

程序不应：

```text
parse README
infer files from README
read Markdown as config
```

README 的角色：

```text
human/AI documentation
```

真正 machine-readable 输入是：

```text
tree file
trait table
CLI parameters
```

---

# 35. 修改代码前的强制原则

如果 AI 判断必须修改代码：

```text
Understand first
→ Search before coding
→ Reuse
→ Extend
→ Refactor
→ Create
```

禁止：

- 第二套 renderer
- 第二套 tree parser
- 第二套 node-ID system
- 第二套 coordinate engine
- 为单一数据集 hard-code species
- 为单一图临时改 topology

---

# 36. 禁止 hard-code 新数据

错误示例：

```python
if species == "Mogera_wogura":
    ...
```

```python
ROOT_ID = 73
```

```python
if project == "myoglobin":
    ...
```

软件必须 Trait-agnostic。

Trait 可以是：

- protein charge
- body mass
- gene expression
- GC content
- PC score
- physiological index
- morphological trait

只要是 numeric node value 即可。

---

# 37. 测试

修改 package 后：

```bash
pytest tests/ -v
```

仅换数据时通常无需修改代码，但至少应运行：

```bash
python -m phylo3d_trait.cli --help
```

并完成真实 CLI smoke test。

---

# 38. AI 最终汇报模板

运行完成后应报告：

```text
Input tree:
<path>

Trait table:
<path>

Tips:
<N>

Internal nodes:
<N>

Branch length interpretation:
<time / substitutions / unknown>

Trait:
<name if known>

Trait range:
<min> – <max>

Output:
<path/to/output.html>

Command used:
<exact command>

Validation:
- all node values present: yes/no
- HTML generated: yes/no
- output size: ...
- warnings: ...
```

不要只说“完成”。

---

# 39. 最小示例

```bash
# 步骤 1: 提取模板
python -m phylo3d_trait.cli template-values \
  --tree data/demo/tree.nwk \
  --output data/demo/template.csv

# 步骤 2: 填入 Trait 并保存为 data/demo/node_values.csv

# 步骤 3: 绘制 3D HTML
python -m phylo3d_trait.cli plot \
  --tree data/demo/tree.nwk \
  --values data/demo/node_values.csv \
  --output results/demo/tree3d.html
```

---

# 40. 核心概念总结

```text
Tree topology
+
branch length / time
+
user-provided node Trait
↓
rectangular 3D phylogeny
```

坐标：

```text
X = Tree layout
Y = Trait
Z = Time before present
```

颜色：

```text
Color = Trait = Y
```

几何：

```text
biological edge
↓
connector + lineage
↓
rectangular phylogram
```

surface：

```text
baseline -> Trait height
```

遮挡：

```text
opaque Mesh3d
+
WebGL depth buffer
```

camera：

```text
orthographic
+
Y-up
+
eLife-style initial view
```

---

# 41. AI 与科研人员极简速查卡片 (Quick Reference & Cheat Sheet)

### 核心坐标与不变量
- **$X$**：水平树布局（Lineage 横向排列）
- **$Y$**：Trait 连续性状高度（**$Y = \text{Trait} = \text{Color}$**）
- **$Z$**：距今演化时间（Tips 为 0，Root 为最大演化年龄）

### 极简执行命令
```bash
# 1. 提取带稳定 Clade ID 的节点模板
python -m phylo3d_trait.cli template-values -t tree.nwk -o values_template.csv

# 2. 用户在 CSV 的 trait 列填入测定值与 ASR 估值，保存为 values.csv

# 3. 渲染 3D 交互 HTML
python -m phylo3d_trait.cli plot -t tree.nwk -v values.csv -o tree3d.html
```

### 关键检查红线
- [x] **不做 ASR**：所有 Tips、Internal Nodes、Root 必须在 CSV 中有数值，不猜测不插值。
- [x] **不猜 Clade ID**：内部节点 ID 一律使用 `template-values` 生成的确定性 `clade:<hash>`。
- [x] **换数据不改代码**：通用 3D 引擎适配任意 Newick/Nexus 树与连续性状。

---

# 42. 软件定位

`Phylo3D-Trait` 是：

> 一个将已知系统发育树与已知节点连续性状映射为交互式三维 rectangular phylogram 的科研可视化工具。

它不是：

- phylogenetic inference software
- ASR software
- sequence analysis software
- statistical model-fitting software
- tree dating software

它负责：

```text
正确读取
→ 正确映射
→ 正确几何
→ 正确颜色
→ 正确交互式三维显示
```

---

# 43. 最终原则

对于 AI：

> **先读说明书，再检查输入，再运行 CLI。新数据默认不改代码。**

对于用户：

> **换 tree + 换 node Trait table，即可复用同一套绘图引擎。**

对于科研解释：

> **只有在 branch length 的意义、节点 Trait 的来源和 ancestral-state reconstruction 方法明确时，图中的时间轴和 Trait 演化轨迹才具有对应的生物学解释。**
