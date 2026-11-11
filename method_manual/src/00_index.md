# 复合湿热浪研究 · 方法论学习手册（14 天）

> **从原始日值到可投稿的方法章** —— 用一套模拟数据，把一篇 SCI 级论文要用的每一个方法
> 从原理到实现走一遍。
>
> ⚠️ **本手册中所有数值结果均来自合成数据，仅用于验证流程与教学。
> 任何数字都不得作为科学结论引用。**

---

## 这个手册是什么

一句话：**它不是一个教程，而是一条已经被完整跑通的流水线。**

每一讲都是一个可独立执行的 Jupyter notebook，从加载数据到出图到落盘，
你打开就能跑，跑完就能看到结果，读完就能照着写论文的对应小节。

| 你会得到 | 位置 |
|---|---|
| 14 个可执行 notebook | `notebooks/S01` … `S14` |
| 31 张论文级图 | `figs/` |
| 13 张结果表（CSV）+ 6 个 RDS | `results/` |
| 一份 SCI 格式的结果汇总 | `results/RESULTS_TABLE.md` |
| 选题的独立评估 | `00_topic_evaluation.md` |

---

## 开始之前

### 环境

| 组件 | 版本 | 说明 |
|---|---|---|
| R | 4.6.1 | 4.3 以上即可 |
| IRkernel | — | `install.packages("IRkernel"); IRkernel::installspec()` |
| Jupyter | — | `jupyter lab` 或 `jupyter notebook` |

必需的 R 包：

```r
pkgs <- c("copula", "ggplot2", "dplyr", "tidyr", "zoo", "MASS",
          "scales", "viridisLite", "patchwork")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
```

### 目录结构

```
UGRP/
├── data/              原始与中间数据
├── method_manual/
│   ├── src/           讲义原文（Markdown，改这里）
│   ├── notebooks/     生成的 notebook（用 md2nb.py 转换）
│   ├── tools/md2nb.py Markdown → ipynb
│   ├── figs/          所有图
│   └── results/       所有表与 RDS
└── scripts/           早期原型脚本（可忽略）
```

### 一个必须知道的细节

讲义是**用 Markdown 写的**，`tools/md2nb.py` 负责转成 notebook。
这样写比手改 JSON 快一个数量级，而且 Git diff 可读。

```r
# 改完 src/S05_event_catalogue.md 之后：
#   python tools/md2nb.py S05
#   cd notebooks && jupyter nbconvert --to notebook --execute --inplace S05_event_catalogue.ipynb
```

> ⚠️ **不要直接改 `.ipynb`** —— 下次转换会被覆盖。永远改 `src/*.md`。

---

## 14 天路线

| Day | Notebook | 主题 | 对应论文章节 | 建议用时 |
|---|---|---|---|---|
| 1 | S01 | 数据审计、质控、湿球温度 | 3.1 | 4 h |
| 2 | S02 | 日历日滑动窗百分位阈值 | 3.2 | 4 h |
| 3 | S03 | 游程理论、去丛聚、事件属性 | 3.3 | 5 h |
| 4 | S04 | 三套定义（Tmax / Tmin / Twmax） | 3.4 | 3 h |
| 5 | S05 | 事件目录、样本量现实检查 | 4.1 | 4 h |
| 6 | S06 | 前期降雨状态 ARP / CDD | 3.5 | 5 h |
| 7 | S07 | RO2：相对风险、置换检验、FDR | 3.6 + 4.2 | 6 h |
| 8 | S08 | Copula 框架、Sklar、τ 的不变性 | 3.7.1 | 5 h |
| 9 | S09 | 边缘分布、截断、三种 PIT | 3.7.2 | 6 h |
| 10 | S10 | copula 拟合、AIC、尾部依赖 | 3.7.3 + 4.3 | 6 h |
| 11 | S11 | RO3：条件依赖与零结果 | 3.7.4 | 5 h |
| 12 | S12 | 联合重现期（AND / OR / Kendall） | 3.8 + 4.4 | 6 h |
| 13 | S13 | 定义敏感性、RDI、排序分歧 | 3.9 + 4.5 | 5 h |
| 14 | S14 | 非平稳性、敏感性、写作模板 | 3.10 + 4.6 + 5 | 6 h |

**总计约 70 小时。** 如果时间紧，优先级是 **S02 → S03 → S07 → S08 → S09 → S10 → S12**（核心方法链），
S04/S05/S13 可以略读，S14 必读。

---

## 怎么用（重要）

每一讲都按同一个节奏组织，**不要跳过任何一步**：

1. **先读「学习目标」** —— 三到五条，告诉你这一讲结束时你应该能回答什么
2. **跑代码，但先猜结果** —— 每个代码块跑之前，先在心里预测输出长什么样
3. **看 🔍 标记的解释** —— 那是「结果为什么是这样」，不是「代码做什么」
4. **看 🚨 标记的警告** —— 那是这一讲最容易出错、也是审稿人最常问的地方
5. **合上 notebook，回答「自检问题」** —— 答不上来就回去重读对应小节

> 🎯 **检验标准**：如果你能把某一讲的「常见错误」表默写出来，
> 并且能解释每一条*为什么*错，这一讲才算学会。

---

## 当前状态

下面这个格子会检查所有 notebook 是否已经执行过、结果文件是否齐全。

```r
for (p in c("00_setup.R", "notebooks/00_setup.R", "../00_setup.R")) {
  if (file.exists(p)) { source(p); break }
}
```

```r
nbs <- sort(list.files(file.path(PATHS$manual, "notebooks"),
                       pattern = "^S[0-9]{2}.*[.]ipynb$", full.names = TRUE))

ok <- vapply(nbs, function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "")
  grepl('"output_type"\\s*:\\s*"(stream|display_data|execute_result)"', txt) &&
    !grepl('"output_type"\\s*:\\s*"error"', txt)
}, logical(1))

status <- data.frame(
  notebook = sub("\\.ipynb$", "", basename(nbs)),
  executed = ifelse(ok, "OK", "未执行或有错"),
  size_kb  = round(file.size(nbs) / 1024, 1),
  row.names = NULL
)
print(status, row.names = FALSE)
cat("\n可执行:", sum(ok), "/", length(ok), "\n")
```

```r
res <- list.files(PATHS$mres, pattern = "\\.(csv|rds)$")
cat("结果文件:", length(res), "个\n")
print(sort(res))
```

---

## 换成真实数据

这是整套手册设计的核心：**换数据不改代码**。

### 需要改的只有一处

`notebooks/00_setup.R` 里的 `read_raw_station()` —— 它负责把原始文件的列名
映射到标准列名。

```r
# 标准列（其余代码全部依赖这 6 列）
# station, date, tmax, tmin, twmax, rain, rh
```

### 湿球温度

`twmax` 是核心变量，但原始站点数据通常**没有**它。两条路：

| 方式 | 说明 | 适用 |
|---|---|---|
| 由 T 和 RH 计算 | `stull_tw()` 已实现 Stull (2011) 公式 | 有逐时/日均 RH |
| 直接采购 | 部分再分析产品提供湿球温度 | 预算允许时 |

> ⚠️ **Stull 公式在 101.325 kPa 下标定**。
> Cameron Highlands（1545 m）等高原站会有系统性偏差 —— S01 会量化这一点，
> 论文里必须交代（或改用气压修正版本）。

### 换数据后的第一步

不要直接跑 S14。按顺序：

1. 跑 **S01** → 看质控报告，确认缺测率与异常值处理
2. 跑 **S02** → 看阈值曲线是否合理（应在年循环上有明显峰谷）
3. 跑 **S05** → **看样本量表**。如果大多数格子 < 30 个事件，
   那么 S11 的路线要按 S05 的三种方案重新选
4. 之后按顺序推进

---

## 每讲产出速查

| Notebook | 图 | 表 |
|---|---|---|
| S02 | 阈值年循环 | — |
| S03 | 游程示意 | — |
| S05 | 事件日历（Gantt） | `table1_event_summary.csv` |
| S06 | ARP–CDD 联合分布 | — |
| S07 | RR 森林图 | `table2_ro2_onset.csv`、`table2b_ro2_attributes.csv` |
| S08 | 六族形状对比 | — |
| S09 | PIT 诊断 | `S09_pit_pairs.rds`、`S09_marginals.rds` |
| S10 | 观测 vs 模拟 | `table3_copula_selection.csv`、`table3b_tail_dependence.csv` |
| S11 | Δτ 森林图、置换零分布 | `table4_conditional_tau.csv` |
| S12 | 重现期等值线 | `table5_return_periods.csv`、`table5b_ensemble.csv` |
| S13 | 排序连线图、RDI 零分布 | `table6_definition_ranks.csv` 等 |
| S14 | 非平稳性对比、敏感性点图 | — |

---

## 常见问题

<details>
<summary>跑某个 notebook 报错 `object 'xxx' not found`</summary>

99% 的原因是**跳着跑**。每个 notebook 的第一段代码是 `source("00_setup.R")`，
它定义了所有共享函数。从头跑（Kernel → Restart & Run All）。

</details>

<details>
<summary>`uniroot` 报 "f() values at end points not of same sign"</summary>

说明要求的重现期超出了该 copula 在这个依赖强度下能表达的范围。
S12 的 `solve_q()` 用 `tryCatch` 返回 `NA` —— **看到 NA 不要忽略**，
先看 `q` 是否已经接近 1。

</details>

<details>
<summary>`gofCopula` 报 "Sn not available for t copulas whose df are not fixed"</summary>

t copula 的 df 是拟合出来的连续值，`gofCopula` 与 `pCopula` 都要求整数 df。
解决办法是复制一份并取整：`tCopula(rho, df = round(df), df.fixed = TRUE)`。
S10 和 S12 都示范了这一点。

</details>

<details>
<summary>图里的中文变成方框</summary>

本手册所有图**只用英文标注**，就是为了避免字体问题。
如果你的系统连英文都显示异常，检查 `theme_paper()` 里的 `base_size` 与设备分辨率。

</details>

<details>
<summary>跑 S14 很慢</summary>

S14 要跑 6 个参数变体，每个约 15–25 秒，总计约 2 分钟。这是正常的。
其余 notebook 都在 30 秒内完成。

</details>

---

## 学完之后你应该能

- [ ] 不看笔记，写出日历日滑动窗百分位阈值的完整算法，并解释为什么要按日历日
- [ ] 解释去丛聚（declustering）解决的是什么问题，以及不做的后果
- [ ] 推导相对风险的分母为什么是「暴露日数」而不是「事件数」
- [ ] 说明 Kendall τ 为什么对单调变换不变，而 Pearson r 会变
- [ ] 写出截断离散边缘的似然，并解释为什么不能和不截断模型比 AIC
- [ ] 说清中点 PIT 与随机化 PIT 的区别，以及各自该用在哪一步
- [ ] 解释 ΔAIC < 2 意味着什么，以及论文里应该怎么写
- [ ] 区分 AND / OR / Kendall 三种重现期，并说明各自的适用场景
- [ ] 解释事件率换算为什么必须用记录长度而不是「有事件的年份数」
- [ ] 说明固定基准期为什么会制造虚假趋势，以及去趋势阈值的代价
- [ ] 面对一个「不显著」的结果，写出三段式解释（效应量 / 检验力 / 不能排除什么）

> 能做到以上十条，这篇论文的方法章就已经在本科水平之上了。
> 剩下的差距不在方法，在**数据质量**与**写作密度**。
