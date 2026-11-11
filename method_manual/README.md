# 复合湿热浪 · 方法论学习手册

一套**已经被完整执行过**的流水线：从逐日气象数据到可投稿的方法章，
用 14 个 Jupyter notebook 覆盖全部方法、原理与常见陷阱。

> ⚠️ 所有结果基于**合成数据**，仅用于验证流程。任何数字都不得作为科学结论引用。

---

## 从哪里开始

| 我想…… | 打开 |
|---|---|
| 先了解路线与环境 | **`notebooks/00_index.ipynb`** |
| 看选题是否值得做 | **`00_topic_evaluation.md`** |
| 直接看结果长什么样 | **`results/RESULTS_TABLE.md`** |
| 开始学习 | `notebooks/S01_data_audit.ipynb` → 按顺序 |

---

## 目录

```
method_manual/
├── 00_topic_evaluation.md    选题的 SCI 标准评估（含风险登记表）
├── src/                      讲义原文（Markdown，改这里）
│   ├── 00_index.md
│   └── S01 … S14 .md
├── notebooks/                生成的 notebook（已全部执行）
│   ├── 00_setup.R            共享模块：所有函数、配置、路径
│   └── S01 … S14 .ipynb
├── tools/md2nb.py            Markdown → ipynb 转换器
├── figs/                     31 张论文级图
└── results/                  结果表（CSV）、RDS、RESULTS_TABLE.md
```

---

## 运行方式

```bash
# 1. 一次性：安装 IRkernel
R -e 'install.packages(c("IRkernel")); IRkernel::installspec()'

# 2. 打开
jupyter lab notebooks/00_index.ipynb

# 3. 从头跑某一讲（Kernel → Restart & Run All）
```

### 修改讲义

讲义以 Markdown 编写，由 `tools/md2nb.py` 转成 notebook。
**永远改 `src/*.md`，不要直接改 `.ipynb`。**

```bash
python tools/md2nb.py S07                       # 只转 S07
cd notebooks && jupyter nbconvert --to notebook --execute --inplace S07_ro2_tests.ipynb
```

---

## 14 天路线

| Day | Notebook | 主题 | 论文章节 |
|---|---|---|---|
| 1 | S01 | 数据审计、质控、湿球温度 | 3.1 |
| 2 | S02 | 日历日滑动窗百分位阈值 | 3.2 |
| 3 | S03 | 游程理论、去丛聚、事件属性 | 3.3 |
| 4 | S04 | 三套定义（Tmax / Tmin / Twmax） | 3.4 |
| 5 | S05 | 事件目录与样本量现实检查 | 4.1 |
| 6 | S06 | 前期降雨状态 ARP / CDD | 3.5 |
| 7 | S07 | RO2：相对风险、置换检验、FDR | 3.6 + 4.2 |
| 8 | S08 | Copula 框架、Sklar、τ 的不变性 | 3.7.1 |
| 9 | S09 | 边缘分布、截断、三种 PIT | 3.7.2 |
| 10 | S10 | copula 拟合、AIC、尾部依赖 | 3.7.3 + 4.3 |
| 11 | S11 | RO3：条件依赖与零结果 | 3.7.4 |
| 12 | S12 | 联合重现期（AND / OR / Kendall） | 3.8 + 4.4 |
| 13 | S13 | 定义敏感性、RDI、排序分歧 | 3.9 + 4.5 |
| 14 | S14 | 非平稳性、敏感性、写作模板 | 3.10 + 4.6 + 5 |

约 70 小时。时间紧时的核心链：**S02 → S03 → S07 → S08 → S09 → S10 → S12**。

---

## 换成真实数据

只需改 `notebooks/00_setup.R` 中的 `read_raw_station()` —— 把原始列名映射到：

```
station, date, tmax, tmin, twmax, rain, rh
```

湿球温度若缺失，用 `stull_tw(T, RH)` 计算（已在 S01 实现并讨论其高原偏差）。
换数据后**先跑 S01 → S02 → S05**，确认质控与样本量，再往后推进。

---

## 状态

- 14/14 notebook 已执行、无误
- 19 个结果文件已生成（13 个 CSV 表 + 6 个 RDS）
- 核心结果汇总见 `results/RESULTS_TABLE.md`
