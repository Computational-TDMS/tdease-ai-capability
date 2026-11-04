# TopPIC workflow execution —— 论文标准执行链（Wilson & Zhou 2022）

> 来源：Wilson JW, Zhou M. **Discovery of Unknown Posttranslational Modifications by Top-Down
> Mass Spectrometry.** In: Sun L, Liu X (eds), *Proteoform Identification: Methods and
> Protocols*, Methods in Molecular Biology vol. 2500, Springer 2022, pp 181–199.
> DOI: 10.1007/978-1-0716-2325-1_13. 细节以原论文为准。
>
> 定位：`parameter-decisions.md` 要求 TopFD/TopPIC 的参数 "set these from the acquisition
> method or cited workflow" —— 本文件就是那个 **cited workflow** 的全参数忠实记录（论文
> Materials + Methods 3.1/3.2），含全部执行侧 Notes（1–15、18）。配置/运行/诊断 TopPIC
> 工具链时用它对照参数与产物；参数语义（单位、节点契约）仍以 `parameter-decisions.md`
> 与实际 Node 为准。

## 1. 链路总览

```text
样品制备(注2) → LC-MS 采集(3.1) → MSConvert(raw→mzML) → TopFD 去卷积 → TopPIC open search
  → 输出表(*_prsm.csv / *_proteoform.csv / Topview *_html) → 目标选择(步骤7)
  → [论文可选半程] Informed-Proteomics: PBF → ProMex(ms1ft) → focused FASTA
     + targeted PTM 列表 → MSPathFinder → LcMsSpectator 人工确认
```

TDEase 分工：MSConvert/TopFD/TopPIC 段 = 本 skill 执行；TopPIC 表解读 = 本 skill 的解读段；
目标选择与家族分诊 = `tdp-opensearch`（总览/attention）；人工确认循环 = `tdp-copilot`。
Informed-Proteomics 的可执行半程（PbfGen → ProMex → MSPathFinder）由
`templates/informed-proteomics-pipeline` 执行，ms1ft feature 列表接总览；LcMsSpectator
的目检 UI 由总览 + tdp-copilot 取代（对照见 §4）。

## 2. 样品制备与采集参数（论文 Materials + Methods 3.1 + Notes 1–9）

### 2.1 仪器与耗材

| 项 | 规格 |
|---|---|
| LC | Waters NanoAquity（或同类）双二元纳流泵：捕获/脱盐 3 μL/min，梯度 0.3 μL/min |
| 捕获柱 | 自装短 C18（Phenomenex Aeris widepore 3.6 μm，ID 150 μm / OD 360 μm × 5 cm）或 C2（3 μm 300 Å，ID 150 μm / OD 360 μm × 5 cm） |
| 分析柱 | 自装短 C18（Phenomenex 3 μm 300 Å，ID 75 μm / OD 360 μm × 70 cm） |
| 流动相 | MPA = 0.1% FA/水；MPB = 0.1% FA/乙腈（全部 MS 级） |
| MS | Thermo Orbitrap Fusion Lumos 或 Eclipse（或同类），配 HCD 与 ETD，高分辨 MS1/MS2 |

公开数据：小鼠 `MZ20160222DS_histone48.raw`（Mendeley Data doi:10.17632/k72vc6vxyv.1）；
高粱 `Sorghum-Histone0810162L01.raw`（PRIDE **PXD014660**）。
软件：MSConvert（ProteoWizard）、Informed-Proteomics（PNNL GitHub）、TopPIC Suite、
（可选）Molecular Weight Calculator（omics.pnl.gov，理论同位素分布）。

### 2.2 采集方法（3.1 全步骤）

1. 进样 1–2 μg 纯化组蛋白。
2. 1% MPB 上样捕获柱 5–10 min @ 3 μL/min，冲盐。
3. 切换到分析柱，启动梯度 + MS 采集。组蛋白典型在 **30–50% MPB** 洗脱；≥100 min 浅梯度跨 30–50% MPB 优化核心组蛋白分离。
4. DDA：**同一前体交替 HCD 与 ETD**。
5. MS1：分辨率 **120k**、**4 microscans**、AGC **1E6**、max injection time **50 ms**。
6. MS2：分辨率 **120k**、1 microscan、AGC **1E6**、隔离窗 **0.6 Da**；ETD 反应时间 **~20 ms**、max IT **500 ms**；HCD **28% NCE、±5% stepped energy**、max IT **100 ms**。动态排除 **120 s**、质量容差 **±0.7 Da**。排除**电荷 <5** 与**未定电荷**态；每个前体**只对一个电荷态**做 MS2。

### 2.3 采集层物理（Notes 1–9）

- **注 1**：组蛋白 11–16 kDa，C18 RPLC 分离良好；方法适用于其它小蛋白。**>20–30 kDa 或更疏水**改用短链 C2（同 MPA/MPB）；先例：蓝藻 phycobilisome 蛋白 17–20 kDa 发色团 PTM [ref 23]。
- **注 2（样品制备层）**：组蛋白纯化 = 核分离 + 组蛋白提取 [ref 24]（酸提取、离子交换——依赖组蛋白强正电荷；核分离的去垢剂浓度按样本优化）。**部分 PTM 酸不稳定：以酸不稳定 PTM 为目标时避免强酸方案** [ref 25]。商品化试剂盒可用（小鼠脑用 ActiveMotif kit；HeLa 组蛋白标准品可做方法开发）。高粱叶：蔗糖密度梯度富集细胞核 [refs 19,26] → 通用**弱阳离子交换（WCX）**纯化；WCX 还去除组蛋白样品常见的高聚物污染、改善 top-down 分析。
- **注 3**：上样量取决于柱尺寸与仪器灵敏度；防 carryover：不过载 + 样品间空白进样；单蛋白/高纯样品减载 5–10 倍。
- **注 4**：该系统死体积导致蛋白出峰在 MS 端**延迟 ~20 min**——梯度开始后延迟启动采集。
- **注 5**：长梯度（如 400 min）增加 MS2 谱数 → 更高 proteoform 鉴定量，尤其低速/低灵敏仪器。
- **注 7**：多 microscans 提升低丰度 S/N，太多拖慢采集、可能降低覆盖。
- **注 8**：**窄隔离窗（0.6 Da）的关键用途**：分离只差 **14 Da** 的甲基化组蛋白形式（多数检出峰上折合 m/z 差 <1）。→ 这就是 TopFD precursor window 降到一个 m/z 的原因。
- **注 9**：多数组蛋白可同位素分辨、机上实时定电荷；更大/更低强度蛋白同位素分布差 → 电荷未定。**靶向高质量/极低丰度蛋白时关掉"排除未定电荷"** 提高选中率。

## 3. 数据处理链（论文 Methods 3.2 步骤 1–13 + Notes 10–15、18）

1. 取样库 FASTA（示例：UniProt 的 *M. musculus* 与 *S. bicolor*）。
2. MSConvert：raw → mzML。
3. TopPIC Suite（GUI 或命令行均可）。
4. **TopFD** 去卷积：默认参数，唯独把 **"precursor window"（−w）降到 1 m/z**——匹配 0.6 Da 窄隔离窗。
5. **TopPIC** 鉴定：默认参数；spectrum 与 proteoform 的 cutoff 类型设 **FDR**、cutoff **0.01**（或按需）；**"proteoform error tolerance" = 5 Da**。输入 = FASTA + TopFD 输出的 `_ms2.msalign`。
   【注 11】5 Da 合并的目的 = 合并相似质量 proteoform、消去单纯去同位素误差造成的冗余；**但不同位点的位置异构体（isomers）也会被合并——慎用**。
6. 看结果：`*_prsm.csv`（全部 PrSM）或 `*_proteoform.csv`（相似质量合并后）；序列覆盖度看输出目录 Topview 模块 `*_html`。
7. **目标选择**：proteoform 列表把 PTM 注释为**方括号内的质量偏移**，定位残基（或残基范围）写在圆括号里。选人工确认目标的判据 = **低 FDR（或低 E-value/P-value）且质量偏移无法用常见 PTM 解释**。【注 12】查身份：UniMod（unimod.org）guest 登录，按质量排序的 PTM 列表里找偏移的候选指认。
8–13. **Informed-Proteomics 预处理（论文可选但推荐；TDEase 中由总览接管，见 §4）**：
   - raw/mzML → **PBF**（LcMsSpectator "Open Raw/MzML" → "XIC View" → "Create PBF now"；或命令行 `PbfGen.exe`）。【注 14】PBF 只为 XIC/feature map 所需，加载 raw 本身不需要。
   - **ProMex** 命令行去卷积 MS1 → **ms1ft** feature 列表（每个 feature = 质量 × RT 的唯一组合）。
   - 用 TopPIC 已鉴定蛋白建 **focused FASTA**，缩小 MSPathFinder 搜索空间。
   - 建 **targeted PTM 列表**。组蛋白常见集：**K 乙酰化、K 单/二/三甲基化、S/T/Y 磷酸化、N 端乙酰化、Met/Cys 氧化**（原文 cystine）；**植物/藻类组蛋白加 N 端单/二/三甲基化**；再加人工分析（步骤 7）发现的 putative PTM。【注 15】**必须精简 FASTA 与 PTM 列表**，否则算不动；算不动就删罕见 PTM——罕见 PTM 靠人工看目标谱确认。**MSPathFinder 的角色是展示主要/常见 proteoform 供人工分析参考，不负责罕见 PTM 确认**。
   - MSPathFinder 默认参数跑，输入 = PBF + ms1ft + focused FASTA + PTM 列表。
   - 【注 10】替代 GUI：MashExplorer [ref 27] 可做去卷积 + TopPIC proteoform 鉴定。【注 13】单谱人工验证工具：ProSightLite、TDValidator；LcMsSpectator 可直接加载 raw 手输序列确认指定谱；但至少用 ProMex 建 feature 文件启用 feature map——**没有 MS2 也能靠 feature map 指认 proteoform**。LcMsSpectator 仅 Windows；非 Windows 用 **PSpecteR**（EMSL-Computing/PSpecteR，R 实现）。

**注 18（两条搜索器的系统偏差）**：
TopPIC open search 把 PTM 注释为质量偏移——发现未知 PTM 的关键；**默认只允许 1 个 unexpected 偏移（上限 2）**，多个 PTM 的**和**会被注释成一个数。**Informed-Proteomics 倾向把未知 PTM 强行匹配成已知 PTM + 端截断的组合**——两者交叉比对提高置信度。TopPIC 的 targeted PTM 列表文件同样可能强行匹配。

## 4. 与 TDEase 的对应（论文环节 → TDEase 组件）

| 论文环节 | 论文工具/参数 | TDEase 对应 |
|---|---|---|
| raw 转换 / 去卷积 / 搜索 | MSConvert / TopFD（−w 1 m/z）/ TopPIC（FDR 0.01, 5 Da） | 本 skill 的 Node 执行；语义与单位见 `parameter-decisions.md`，产物清单见 `suite-and-files.md` |
| PrSM/proteoform 表解读 | `*_prsm.csv`、`*_proteoform.csv`、Topview | 本 skill（表语法、评分语义） |
| feature map（RT×质量）、家族/梯子分诊 | ProMex ms1ft + LcMsSpectator | `templates/informed-proteomics-pipeline` 产出 ms1ft → `tdp-opensearch` 总览（inspect/ingest/`families`、`attention_score`、familyRefs） |
| 目标选择（步骤 7 判据） | 低 FDR + 常见 PTM 闭不开 | tdp-opensearch attention 面 + E-value zone；交接给 tdp-copilot |
| MSPathFinder targeted 搜索 | focused FASTA + 精简 PTM 列表 | `templates/informed-proteomics-pipeline` 执行；`tdp-copilot` 的 targeted 边界——其已知 PTM 拟合只作交叉核对，不替换 TopPIC open-search shift |
| 单谱人工验证 | ProSightLite / TDValidator / LcMsSpectator 手改拟合 | `tdp-copilot` 确认循环 + `spectral-evidence`；对称 manifest 比较取代自由拟合 |
| 理论同位素分布 | Molecular Weight Calculator | Rust 化学权威 `tdeasenexus-toolbox theory` / `fit-envelope` |

## 5. 本文件承载的 Notes 清单（1–15、18）

| # | 一句话 | 本文位置 |
|---|---|---|
| 1 | C18 适合 11–16 kDa；>20–30 kDa/疏水用 C2；phycobilisome 先例 | §2.3 |
| 2 | 组蛋白纯化：核分离+提取；酸不稳定 PTM 避强酸；WCX 去高聚物 | §2.3 |
| 3 | 上样量优化；防 carryover；纯样减载 5–10× | §2.3 |
| 4 | 死体积 ~20 min 延迟；采集延后启动 | §2.3 |
| 5 | 长梯度（400 min）增 MS2 数与鉴定量 | §2.3 |
| 6 | ETD 端覆盖+不稳定 PTM 完整；HCD 序列特异+诊断丢失 | §2.2 |
| 7 | microscans：S/N vs 循环时间权衡 | §2.3 |
| 8 | 0.6 Da 窄窗隔离 14 Da 甲基化系列 | §2.3 |
| 9 | 电荷<5/未定排除；高质量/低丰度靶向时关闭 | §2.3 |
| 10 | MashExplorer 替代 GUI | §3（步骤8–13） |
| 11 | TopPIC 5 Da 合并容差：去同位素冗余**但连位置异构体也合并**——慎用 | §3（步骤5） |
| 12 | UniMod guest 按质量查 | §3（步骤7） |
| 13 | ProSightLite/TDValidator 单谱验证；PSpecteR 非 Windows；ProMex feature map 无 MS2 也可指认 | §3（步骤8–13） |
| 14 | PBF 只为 XIC/feature map；PbfGen.exe 可替代 | §3（步骤8–13） |
| 15 | MSPathFinder 必须精简 FASTA/PTM 列表；其角色是展示常见 proteoform | §3（步骤8–13） |
| 18 | open search 1 个 unexpected（上限2）→和值注释；Informed-Proteomics 强行匹配已知组合；交叉比对；targeted 列表同样会强行匹配 | §3（注18） |

## 6. OCR 修复说明（对全文提取文本适用）

PDF 提取文本的连字/符号乱码规律：`/g415` = "ti"（Sec/g415on=Section）；`/C6` = "±"（±0.7 Da、±5% stepped）；`/C0` = "−"（−18.01、−w）。
据此修复的三处参数：HCD = **28% NCE ±5% stepped energy**；动态排除容差 = **±0.7 Da**；TopFD **precursor window（−w）= 1 m/z**。
