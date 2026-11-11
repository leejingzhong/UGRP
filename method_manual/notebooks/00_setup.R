# ============================================================================
# 00_setup.R  —  共享环境 / 常量 / 工具函数
#
# 每一个 notebook 的第一格都只做一件事：
#     source("00_setup.R")   （或自动搜索到它）
#
# 设计原则：
#   1. 所有可调参数集中在 CFG 里，正文代码里不出现魔法数字。
#   2. 所有函数只做一件事，并且有"为什么"注释。
#   3. 路径自动向上搜索项目根目录，所以在 UGRP/ 或 method_manual/notebooks/
#      下启动 Jupyter 都能跑。
# ============================================================================

# ---------------------------------------------------------------------------
# 0. 项目根目录探测
# ---------------------------------------------------------------------------
.find_root <- function(start = getwd(), max_up = 6) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_up)) {
    if (dir.exists(file.path(d, "02_raw")) && dir.exists(file.path(d, "scripts"))) return(d)
    p <- dirname(d)
    if (identical(p, d)) break
    d <- p
  }
  stop("找不到项目根目录（需含 02_raw/ 与 scripts/）。",
       "请把 Jupyter 的工作目录设在 UGRP/ 或其子目录下。")
}

ROOT <- .find_root()

PATHS <- list(
  root   = ROOT,
  raw    = file.path(ROOT, "02_raw"),
  clean  = file.path(ROOT, "03_clean"),
  events = file.path(ROOT, "04_events"),
  res    = file.path(ROOT, "05_results"),
  figs   = file.path(ROOT, "06_figures"),
  manual = file.path(ROOT, "method_manual"),
  mres   = file.path(ROOT, "method_manual", "results"),
  mfig   = file.path(ROOT, "method_manual", "figs")
)
invisible(lapply(PATHS[c("clean", "events", "res", "figs", "mres", "mfig")],
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# ---------------------------------------------------------------------------
# 1. 全局参数（所有可调项都在这里，正文不写魔法数字）
# ---------------------------------------------------------------------------
CFG <- list(
  # --- 阈值法 (S02) ---
  base_years  = 1985:2014,   # ETCCDI 推荐基准期
  halfwin     = 7,           # 日历日 ±7 天滑动窗 (共 15 天)
  prob_default= 0.90,        # 90 百分位
  min_obs     = 50,          # 滑动窗内最少有效观测数

  # --- 事件提取 (S03) ---
  min_len     = 3,           # 最短持续 3 天
  max_gap     = 2,           # ≤2 天中断合并为同一事件 (declustering)

  # --- 前期降雨 (S06) ---
  arp_window  = 7,           # 事件起始前 7 日累计降雨
  arp_dry     = 30,          # ARP < 30 → dry
  arp_wet     = 50,          # ARP > 50 → wet   (30–50 为中性带)

  # --- Copula (S08–S11) ---
  min_n_copula = 30,         # 单格最少事件数，低于此不拟合
  n_boot       = 1000,       # bootstrap / 置换次数
  boot_seed    = 11,

  # --- 复现性 ---
  seed = 20260101
)

DEFS <- c("Tmax", "Tmin", "Twmax")
# 定义名 → 数据列名。三套定义共用同一套方法，只是换了输入序列。
VAR_OF <- c(Tmax = "tmax", Tmin = "tmin", Twmax = "twmax")

# ---------------------------------------------------------------------------
# 2. 依赖包
# ---------------------------------------------------------------------------
.pkgs <- c("copula", "ggplot2", "dplyr", "tidyr", "lubridate",
           "zoo", "fitdistrplus", "scales")
for (p in .pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("正在安装缺失包: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(copula); library(ggplot2); library(dplyr)
  library(tidyr);  library(lubridate); library(zoo)
  library(fitdistrplus); library(scales)
})

# ---------------------------------------------------------------------------
# 3. 小工具
# ---------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

BANNER <- function(x) {
  cat("\n", strrep("=", 74), "\n ", x, "\n", strrep("=", 74), "\n", sep = "")
}

# 环形日历日距离：12 月 31 日与 1 月 1 日只差 1 天，不是 364 天
circ_dist <- function(a, b, period = 365) {
  d <- abs(a - b)
  pmin(d, period - d)
}

# ---------------------------------------------------------------------------
# 4. 湿球温度 (Stull 2011)
#    物理：湿球温度是"水蒸发冷却能达到的最低温度"，
#          因此它同时编码了气温与湿度 → 是人体热应激的单一指标。
#    适用范围：RH 5–99%，T −20~50 °C，误差 < 0.3 °C
# ---------------------------------------------------------------------------
stull_tw <- function(T, RH) {
  T * atan(0.151977 * sqrt(RH + 8.313659)) +
    atan(T + RH) -
    atan(RH - 1.676331) +
    0.00391838 * RH^1.5 * atan(0.023101 * RH) -
    4.686035
}

# ---------------------------------------------------------------------------
# 5. 日历日滑动窗百分位阈值 (S02)
#    为什么必须"按日历日"：如果直接用全年分位数，冬天的低温会把
#    夏天的阈值拉低，导致夏季虚假事件暴增。
#    为什么用 type = 8：ETCCDI / WMO 推荐的无偏分位数估计。
# ---------------------------------------------------------------------------
daily_threshold <- function(doy, values, years,
                            base_years = CFG$base_years,
                            prob       = CFG$prob_default,
                            halfwin    = CFG$halfwin,
                            min_obs    = CFG$min_obs,
                            period     = 365) {
  in_base <- years %in% base_years
  th <- rep(NA_real_, period)
  for (d in seq_len(period)) {
    sel <- in_base & circ_dist(doy, d, period) <= halfwin & !is.na(values)
    if (sum(sel) >= min_obs) {
      th[d] <- quantile(values[sel], prob, type = 8, na.rm = TRUE)
    }
  }
  idx <- which(!is.na(th))
  if (length(idx) > 1) {
    th <- approx(c(idx - period, idx, idx + period),
                 rep(th[idx], 3),
                 xout = seq_len(period), rule = 2)$y
  }
  th
}

# 把 365 天的阈值映射回每一行
map_threshold <- function(th, doy) th[doy]

# ---------------------------------------------------------------------------
# 6. 游程理论 (S03)
# ---------------------------------------------------------------------------
find_runs <- function(flag) {
  r <- rle(flag)
  data.frame(value  = r$values,
             start  = cumsum(r$lengths) - r$lengths + 1,
             end    = cumsum(r$lengths),
             length = r$lengths)
}

# declustering：把间隔 ≤ max_gap 天的两个超阈值片段合并成一个事件。
# 为什么：热浪常常"热三天、缓一天、再热三天"，若不合并会被算成两个事件，
#         导致事件数虚高、duration 系统性偏短。
decluster <- function(runs, min_len = CFG$min_len, max_gap = CFG$max_gap) {
  runs <- runs[runs$value == 1 & runs$length >= min_len, , drop = FALSE]
  if (nrow(runs) == 0) return(runs)
  out <- runs[1, , drop = FALSE]
  for (i in seq_len(nrow(runs))[-1]) {
    if (runs$start[i] - out$end[nrow(out)] - 1 <= max_gap) {
      out$end[nrow(out)]    <- runs$end[i]
      out$length[nrow(out)] <- out$end[nrow(out)] - out$start[nrow(out)] + 1
    } else {
      out <- rbind(out, runs[i, , drop = FALSE])
    }
  }
  out
}

# 缺测日必须硬性切断事件：不能让缺失数据把两个事件"缝"起来
split_at_gaps <- function(ev, bad_idx) {
  if (nrow(ev) == 0) return(ev)
  out <- list()
  for (i in seq_len(nrow(ev))) {
    s <- ev$start[i]; e <- ev$end[i]
    gaps <- bad_idx[bad_idx > s & bad_idx < e]
    cuts <- c(s - 1, gaps, e)
    for (j in seq_along(cuts)[-1]) {
      a <- cuts[j - 1] + 1; b <- cuts[j] - 1
      if (b >= a) out[[length(out) + 1]] <-
        data.frame(value = 1L, start = a, end = b, length = b - a + 1L)
    }
  }
  if (length(out) == 0) return(ev[0, , drop = FALSE])
  do.call(rbind, out)
}

# 单站单变量的完整事件提取（RO1 的核心）
extract_events <- function(df, var, defn = var,
                           prob = CFG$prob_default, cfg = CFG) {

  df <- df[format(df$date, "%m-%d") != "02-29", , drop = FALSE]
  df <- df[order(df$date), , drop = FALSE]

  doy  <- as.integer(format(df$date, "%j"))
  yr   <- as.integer(format(df$date, "%Y"))
  vals <- df[[var]]

  th  <- daily_threshold(doy, vals, yr, cfg$base_years, prob, cfg$halfwin, cfg$min_obs)
  thv <- map_threshold(th, doy)

  ex  <- as.integer(vals > thv)
  bad <- which(is.na(ex))
  ex[bad] <- 0L

  runs <- find_runs(ex)
  if (nrow(runs) == 0) return(NULL)
  runs <- runs[runs$value == 1, , drop = FALSE]
  if (nrow(runs) == 0) return(NULL)

  ev <- decluster(runs, cfg$min_len, cfg$max_gap)
  ev <- split_at_gaps(ev, bad)
  ev <- ev[ev$length >= cfg$min_len, , drop = FALSE]
  if (nrow(ev) == 0) return(NULL)

  n <- nrow(ev)
  duration <- ev$length
  peak_value <- intensity <- severity <- threshold_at_peak <- numeric(n)
  peak_date  <- as.Date(rep(NA_character_, n))
  n_missing  <- integer(n)
  n_exceed   <- integer(n)

  for (i in seq_len(n)) {
    idx <- ev$start[i]:ev$end[i]
    v   <- vals[idx]; t <- thv[idx]
    k   <- which.max(v)
    peak_value[i]         <- v[k]
    peak_date[i]          <- df$date[idx[k]]
    intensity[i]          <- v[k] - t[k]          # 强度 = 峰值超出量
    threshold_at_peak[i]  <- t[k]
    severity[i]           <- sum(pmax(v - t, 0), na.rm = TRUE)  # 累积超出量
    n_missing[i]          <- sum(is.na(v))
    n_exceed[i]           <- sum(v > t, na.rm = TRUE)
  }

  data.frame(
    station      = df$station[1],
    defn         = defn,
    start_date   = df$date[ev$start],
    end_date     = df$date[ev$end],
    start_year   = as.integer(format(df$date[ev$start], "%Y")),
    start_doy    = doy[ev$start],
    duration     = duration,
    peak_date    = as.Date(peak_date),
    peak_value   = peak_value,
    threshold_at_peak = threshold_at_peak,
    intensity    = intensity,
    severity     = severity,
    n_exceed     = n_exceed,
    n_missing    = n_missing,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 7. 前期降雨状态 (S06)
#    核心思想：用"相对该站该季节的历史位置（百分位）"而不是毫米数，
#              因为 5 mm 在东海岸和西海岸的意义完全不同。
# ---------------------------------------------------------------------------
antecedent_series <- function(d, cfg = CFG) {

  d <- d[format(d$date, "%m-%d") != "02-29", , drop = FALSE]
  d <- d[order(d$date), , drop = FALSE]

  doy <- as.integer(format(d$date, "%j"))
  yr  <- as.integer(format(d$date, "%Y"))
  W   <- cfg$arp_window

  # 截止到"前一天"的 W 日累计降雨（严格前期，不含当日）
  rw_now  <- zoo::rollsum(d$rain, k = W, align = "right", fill = NA)
  # round(): rollsum 走 cumsum，会残留 -3.8e-14 这类浮点噪声
  rw_prev <- round(c(NA_real_, rw_now[-length(rw_now)]), 6)

  # 连续干日（截至前一天）
  cdd <- numeric(nrow(d))
  for (i in 2:nrow(d)) {
    cdd[i] <- if (!is.na(d$rain[i - 1]) && d$rain[i - 1] < 1) cdd[i - 1] + 1 else 0
  }

  # ARP：对每个日历日，用基准期内 ±halfwin 天的历史分布做参照
  in_base <- yr %in% cfg$base_years
  arp <- rep(NA_real_, nrow(d))
  for (dd in seq_len(365)) {
    idx <- which(doy == dd)
    if (length(idx) == 0) next
    ref <- which(in_base & circ_dist(doy, dd) <= cfg$halfwin & !is.na(rw_prev))
    if (length(ref) < 30) next
    Fn <- ecdf(rw_prev[ref])
    arp[idx] <- 100 * Fn(rw_prev[idx])
  }

  d$rain_win_prev <- rw_prev
  d$cdd           <- cdd
  d$arp           <- arp
  d$state <- cut(arp,
                 breaks = c(-Inf, cfg$arp_dry, cfg$arp_wet, Inf),
                 labels = c("dry", "neutral", "wet"))
  d
}

# 把事件起始日的状态挂到事件表上（留一份"不含事件当天"的证据）
attach_antecedent <- function(ev, ant) {
  m <- match(ev$start_date, ant$date)
  ev$arp            <- ant$arp[m]
  ev$cdd            <- ant$cdd[m]
  ev$rain_win_prev  <- ant$rain_win_prev[m]
  ev$state          <- ant$state[m]
  ev
}

# ---------------------------------------------------------------------------
# 8. Copula 工具 (S08–S11)
# ---------------------------------------------------------------------------
# 离散变量（duration 是整数）的 PIT：在 F(x-) 与 F(x) 之间均匀随机取一点。
# 为什么：连续的 copula 理论要求 U ~ Unif(0,1)；整数变量直接取 ecdf 会
#         在格点上堆积质量，导致拟合出现伪装的完美拟合。
pit_discrete <- function(x, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Fn <- ecdf(x)
  u  <- runif(length(x))
  Fn(x - 1e-9) + u * (Fn(x) - Fn(x - 1e-9))
}

# 连续变量的 PIT（中位秩，避免 0 和 1）
pit_rank <- function(x) (rank(x) - 0.5) / length(x)

tau_of <- function(x, y) suppressWarnings(cor(x, y, method = "kendall"))

# 尾部依赖系数的解析式（copula 包不会自动给 Archimedean 族，必须自己写）
#   Clayton: lambda_L = 2^(-1/theta)
#   Gumbel : lambda_U = 2 - 2^(1/theta)
#   Joe    : lambda_U = 2 - 2^(1/theta)
#   t      : lambda = 2 * t_{nu}(-sqrt((nu+1)(1-rho)/(1+rho)))
#   Gaussian / Frank: 都为 0（这正是"它们拟合不了尾部依赖"的原因）
tail_dep <- function(fam, cop) {
  p  <- tryCatch(as.numeric(cop@parameters), error = function(e) NA_real_)
  th <- p[1]
  out <- c(lower = 0, upper = 0)
  if (is.na(th)) return(c(lower = NA_real_, upper = NA_real_))
  switch(fam,
    clayton = { out["lower"] <- 2^(-1 / th) },
    gumbel  = { out["upper"] <- 2 - 2^(1 / th) },
    joe     = { out["upper"] <- 2 - 2^(1 / th) },
    t       = {
      rho <- th
      df  <- if (length(p) > 1) p[2] else 4
      lam <- 2 * stats::pt(-sqrt((df + 1) * (1 - rho) / (1 + rho)), df)
      out[] <- lam
    },
    NULL)
  out
}

# 统一接口：给族名，返回一个未拟合的 copula 对象
make_copula <- function(fam, dim = 2) switch(fam,
  gaussian = normalCopula(dim = dim),
  t        = tCopula(dim = dim, df = 4, df.fixed = FALSE),
  clayton  = claytonCopula(dim = dim),
  gumbel   = gumbelCopula(dim = dim),
  frank    = frankCopula(dim = dim),
  joe      = joeCopula(dim = dim),
  stop("未知族: ", fam))

COPULA_FAMILIES <- c("gaussian", "t", "clayton", "gumbel", "frank", "joe")

fit_family <- function(u1, u2, fam) {
  tryCatch(
    fitCopula(make_copula(fam), cbind(u1, u2), method = "ml"),
    error = function(e) NULL
  )
}

# 一个分层（例如 dry）的完整 copula 分析
analyse_stratum <- function(dat, label, cfg = CFG, verbose = TRUE) {
  if (nrow(dat) < cfg$min_n_copula) {
    if (verbose) cat(sprintf("\n[%s] n = %d < %d，样本不足，跳过拟合\n",
                             label, nrow(dat), cfg$min_n_copula))
    return(NULL)
  }
  set.seed(cfg$boot_seed)
  u1 <- pit_discrete(dat$duration)
  u2 <- pit_rank(dat$intensity)

  fits <- lapply(COPULA_FAMILIES, function(f) fit_family(u1, u2, f))
  names(fits) <- COPULA_FAMILIES

  loglik <- vapply(fits, function(f) if (is.null(f)) NA_real_ else as.numeric(logLik(f)),
                   numeric(1))
  kpar   <- vapply(fits, function(f) if (is.null(f)) NA_real_ else length(f@estimate),
                   numeric(1))
  aic    <- -2 * loglik + 2 * kpar

  tab <- data.frame(family = COPULA_FAMILIES,
                    loglik = round(loglik, 2),
                    k      = kpar,
                    AIC    = round(aic, 2))
  tab <- tab[order(tab$AIC), ]

  best_name <- tab$family[1]
  best      <- fits[[best_name]]
  if (is.null(best)) return(NULL)

  tau <- tau_of(dat$duration, dat$intensity)
  td  <- tail_dep(best_name, best@copula)

  if (verbose) {
    cat(sprintf("\n[%s]  n = %d\n", label, nrow(dat)))
    print(tab, row.names = FALSE)
    cat(sprintf("  Kendall tau = %.3f   最优族 = %s\n", tau, best_name))
    cat(sprintf("  lambda_L = %.3f   lambda_U = %.3f\n", td["lower"], td["upper"]))
  }

  list(cop = best@copula, fit = best, family = best_name,
       tau = tau, lambda_L = unname(td["lower"]), lambda_U = unname(td["upper"]),
       n = nrow(dat), aic_table = tab, u = u1, v = u2)
}

# ---------------------------------------------------------------------------
# 9. 联合重现期 (S12)
# ---------------------------------------------------------------------------
T_AND <- function(u, v, cop) 1 / (1 - u - v + pCopula(cbind(u, v), cop))
T_OR  <- function(u, v, cop) 1 / (1 - pCopula(cbind(u, v), cop))

# Kendall 重现期：K_C(t) = P(C(U,V) <= t)，
# 用蒙特卡洛估计，对任何 copula 族都成立（不必逐族推解析式）
kendall_C_fun <- function(cop, n = 20000, seed = 20260101) {
  set.seed(seed)
  u <- rCopula(n, cop)
  w <- pCopula(u, cop)
  stats::ecdf(w)
}
T_KENDALL <- function(C, Kfun) 1 / (1 - Kfun(C))

# 把离散边缘 CDF 线性插值成连续函数，否则等值线会变成阶梯
smooth_cdf <- function(x) {
  xs <- sort(unique(x))
  Fs <- ecdf(x)(xs)
  Fs <- pmin(Fs, 1 - 1 / length(x))
  approxfun(xs, Fs, method = "linear", rule = 2, ties = "ordered")
}

# ---------------------------------------------------------------------------
# 10. 排序分歧指数 RDI (S13)
#    定义：sum|r_A(i) - r_B(i)| / (n^2 / 2)
#    含义：两套排序完全一致 → 0；完全反转 → 1
# ---------------------------------------------------------------------------
rdi <- function(rankA, rankB) {
  n <- length(rankA)
  sum(abs(rankA - rankB)) / (n^2 / 2)
}

# ---------------------------------------------------------------------------
# 11. 绘图规范
# ---------------------------------------------------------------------------
theme_paper <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor   = element_blank(),
          plot.title         = element_text(face = "bold", size = base_size + 1),
          plot.subtitle      = element_text(colour = "grey30"),
          strip.background   = element_rect(fill = "grey92", colour = NA),
          legend.position    = "bottom")
}

# 所有用模拟数据画的图都必须带这个水印，防止结果被误读为科学结论
watermark <- function(label = "SYNTHETIC DATA - pipeline validation only") {
  annotate("text", x = Inf, y = -Inf, label = label,
           hjust = 1.02, vjust = -0.6, size = 3, colour = "grey55", fontface = "italic")
}

save_fig <- function(p, name, w = 6.5, h = 4.5, dpi = 200) {
  f <- file.path(PATHS$mfig, paste0(name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = dpi, bg = "white")
  invisible(f)
}

# ---------------------------------------------------------------------------
# 12. 数据准备：QC → 湿球 → 落盘
#
#     真数据到位后，唯一需要改的是 read_raw_station() 里的列名映射。
#     其余全部保持不动 —— 这是"换数据不改代码"的关键。
# ---------------------------------------------------------------------------
MISSING_CODES <- c(-999, -9999, 9999, 99999)

# 单站质控。每一步都记录被改动的数量，方便在 S01 里做质控报告。
qc_one_station <- function(df, cfg = CFG) {
  num <- intersect(c("tmax", "tmin", "rain", "rh"), names(df))
  n_before <- nrow(df)

  # (1) 缺测编码统一成 NA
  n_miss_code <- 0
  for (v in num) {
    hit <- df[[v]] %in% MISSING_CODES
    n_miss_code <- n_miss_code + sum(hit)
    df[[v]][hit] <- NA
  }

  # (2) 日期：去闰日、排序、去重
  df$date <- as.Date(df$date)
  df <- df[!is.na(df$date), , drop = FALSE]
  df <- df[format(df$date, "%m-%d") != "02-29", , drop = FALSE]
  df <- df[order(df$date), , drop = FALSE]
  n_dup <- sum(duplicated(df$date))
  df <- df[!duplicated(df$date), , drop = FALSE]

  # (3) 物理一致性
  bad_range <- rep(FALSE, nrow(df))
  if ("tmax" %in% num) bad_range <- bad_range | (!is.na(df$tmax) & (df$tmax < -20 | df$tmax > 50))
  if ("tmin" %in% num) bad_range <- bad_range | (!is.na(df$tmin) & (df$tmin < -20 | df$tmin > 50))
  if ("rh"   %in% num) bad_range <- bad_range | (!is.na(df$rh)   & (df$rh   <   0 | df$rh   > 100))
  for (v in intersect(c("tmax", "tmin", "rh"), num)) df[[v]][bad_range] <- NA

  bad_neg <- if ("rain" %in% num) !is.na(df$rain) & df$rain < 0 else rep(FALSE, nrow(df))
  df$rain[bad_neg] <- NA

  # 顺序颠倒时两侧都置 NA：不猜哪个是真值
  bad_ord <- if (all(c("tmax", "tmin") %in% num))
    !is.na(df$tmax) & !is.na(df$tmin) & df$tmax < df$tmin else rep(FALSE, nrow(df))
  df$tmax[bad_ord] <- NA; df$tmin[bad_ord] <- NA

  # (4) 湿球温度
  #     注意：这里用 stull_tw(tmax, rh) 是"日最高气温 + 日均相对湿度"的
  #     近似。真实研究里 Twmax 应直接向 MetMalaysia 购买（见 S04）。
  if (all(c("tmax", "rh") %in% names(df))) {
    df$twmax <- stull_tw(df$tmax, df$rh)
  }

  attr(df, "qc") <- c(rows_in = n_before, miss_code = n_miss_code, dup = n_dup,
                      range = sum(bad_range), neg_rain = sum(bad_neg),
                      order_bad = sum(bad_ord))
  df
}

# 读一个站的原始文件。真数据到手后改这个函数即可。
read_raw_station <- function(f) {
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

load_stations_meta <- function() {
  f <- file.path(PATHS$raw, "stations_meta.csv")
  if (!file.exists(f)) return(data.frame(station = character(0)))
  read.csv(f, stringsAsFactors = FALSE)
}

# 构建清洗数据；已存在且完整时直接读缓存
build_clean_daily <- function(force = FALSE) {
  out  <- file.path(PATHS$clean, "clean_daily.csv")
  need <- c("station", "date", "tmax", "tmin", "rain", "rh", "twmax")
  n_st <- nrow(load_stations_meta())

  if (!force && file.exists(out)) {
    d <- read.csv(out, stringsAsFactors = FALSE)
    if (all(need %in% names(d)) && length(unique(d$station)) >= max(n_st, 1)) {
      d$date <- as.Date(d$date)
      return(d)
    }
  }

  files <- list.files(PATHS$raw, pattern = "\\.csv$", full.names = TRUE)
  files <- files[!grepl("stations_meta", files)]
  raw   <- do.call(rbind, lapply(files, read_raw_station))

  parts <- lapply(split(raw, raw$station), qc_one_station)
  qc_tab <- do.call(rbind, lapply(names(parts), function(s) {
    q <- attr(parts[[s]], "qc")
    data.frame(station = s, t(as.data.frame(q)), row.names = NULL)
  }))

  d <- do.call(rbind, parts)
  write.csv(d, out, row.names = FALSE)
  write.csv(qc_tab, file.path(PATHS$mres, "qc_summary.csv"), row.names = FALSE)
  d
}

load_clean_daily <- function(force = FALSE) build_clean_daily(force)

# 全站 × 三定义的完整事件目录（RO1 的最终产物），带缓存
build_event_catalogue <- function(force = FALSE, cfg = CFG) {
  f <- file.path(PATHS$mres, "events_with_antecedent.rds")
  if (!force && file.exists(f)) return(readRDS(f))

  d <- load_clean_daily()
  out <- list()
  for (st in unique(d$station)) {
    s   <- subset(d, station == st)
    ant <- antecedent_series(s, cfg)
    evs <- lapply(DEFS, function(dn) extract_events(s, VAR_OF[[dn]], defn = dn, cfg = cfg))
    evs <- do.call(rbind, Filter(Negate(is.null), evs))
    if (!is.null(evs) && nrow(evs) > 0) out[[st]] <- attach_antecedent(evs, ant)
  }
  ev <- do.call(rbind, out)
  ev$id <- seq_len(nrow(ev))
  saveRDS(ev, f)
  write.csv(ev, file.path(PATHS$mres, "events_with_antecedent.csv"), row.names = FALSE)
  ev
}

load_catalogue <- function(force = FALSE) build_event_catalogue(force)

cat("\n[00_setup.R] 已加载\n")
cat("  项目根目录 : ", PATHS$root, "\n", sep = "")
cat("  基准期     : ", min(CFG$base_years), "-", max(CFG$base_years), "\n", sep = "")
cat("  阈值       : ", CFG$prob_default * 100, "百分位, ±", CFG$halfwin, "天滑窗\n", sep = "")
cat("  事件规则   : 最短 ", CFG$min_len, " 天, 合并间隔 ≤ ", CFG$max_gap, " 天\n", sep = "")
cat("  干/湿切点  : ARP < ", CFG$arp_dry, " / > ", CFG$arp_wet, "\n", sep = "")
