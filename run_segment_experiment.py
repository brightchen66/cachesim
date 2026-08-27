# -*- coding: utf-8 -*-
"""分段 LocalRatioCaching 实验：有界离线视野下 4-近似算法的退化规律。

将整条 trace 均分为 k 段（k=10/20/50/100），段内独立运行修改后的
LocalRatioCaching（死页按“下次请求为该段第 n_seg+1 次”建虚拟实例），
把各段局部解直接拼接成合并解，与整条 trace 上的全局解比较。

模型口径（分段 = 有界离线视野 / batch 化）：
- 每段独立、从空缓存出发：算法只能看到段内请求，段间不共享缓存状态与调度信息；
- 合并 = 拼接：段首请求若对象在上一段末尾仍在缓存，也按冷启动计 miss——
  这是分段信息损失的主要来源（段边界的 compulsory miss）；
- 段内“请求一次就不再请求”的死页按段末建虚拟实例（全局运行按 trace 末尾建）。

冷启动下界（任何分段方案的固有损失，与段内算法好坏无关）：
  miss >= Σ_段 段内不同对象数；bit 模型 cost >= Σ_段 段内不同对象字节和。
用于区分“分段固有损失”与“算法在段内的容量压力损失”。

用法:
    python run_segment_experiment.py                          # 默认 twitter29_test10.csv
    python run_segment_experiment.py twitter45_test.csv wiki2018_test10.csv
    python run_segment_experiment.py --pcts 1 5 --segments 10 100
"""

import argparse
from pathlib import Path

from cache_sim.algorithms.local_ratio_caching import LocalRatioCaching
from cache_sim.algorithms.registry import get_algorithm
from cache_sim.engine.simulator import ContentSimulator
from cache_sim.traceparser.loader import check_trace, iter_requests

TRACE_DIR = Path("traces")
STAT_KEYS = ("hits", "misses", "byte_total", "byte_hit", "evictions",
             "cost", "rounds", "vinst")


def split_segments(items, k):
    """将 items 均分为 k 段：前 n%k 段各多 1 条请求。"""
    n = len(items)
    base, rem = n // k, n % k
    out, start = [], 0
    for i in range(k):
        size = base + (1 if i < rem else 0)
        if size > 0:
            out.append(items[start:start + size])
        start += size
    return out


def seg_floor(segs):
    """各段冷启动下界：(Σ段内不同对象数, Σ段内不同对象字节和)。"""
    floor_miss = 0
    floor_cost = 0.0
    for seg in segs:
        first = {}
        for _t, oid, sz, _e in seg:
            first.setdefault(oid, sz)
        floor_miss += len(first)
        floor_cost += float(sum(first.values()))
    return floor_miss, floor_cost


def run_one(trace, capacity, cost_model):
    """单次运行（全局或单段），返回汇总指标。"""
    res = LocalRatioCaching().run(trace, capacity=capacity, cost_model=cost_model)
    return {
        "hits": res.hits,
        "misses": res.misses,
        "byte_total": res.byte_total,
        "byte_hit": res.byte_hit,
        "evictions": res.evictions,
        "cost": res.extra["cost"],
        "rounds": res.extra["num_rounds"],
        "vinst": res.extra["num_virtual_instances"],
    }


def run_segmented(trace, capacity, k, cost_model):
    """k 段独立求解后拼接：各段指标求和 + 段信息 + 冷启动下界。"""
    segs = split_segments(trace, k)
    tot = {key: 0 for key in STAT_KEYS}
    for seg in segs:
        r = run_one(seg, capacity, cost_model)
        for key in STAT_KEYS:
            tot[key] += r[key]
    tot["num_segments"] = len(segs)
    tot["seg_len_min"] = min(len(s) for s in segs)
    tot["seg_len_max"] = max(len(s) for s in segs)
    tot["floor_miss"], tot["floor_cost"] = seg_floor(segs)
    return tot


def _seg_len_txt(s):
    if s["seg_len_min"] == s["seg_len_max"]:
        return str(s["seg_len_min"])
    return f"{s['seg_len_min']}~{s['seg_len_max']}"


def run_dataset(path, pcts, segments, cost_model):
    stats = check_trace(path)
    trace = [(t, oid, sz, e) for t, oid, sz, e in iter_requests(path)]
    n = len(trace)
    print("=" * 108)
    print(f"=== {path.name}: n={n}, 唯一对象={stats['num_unique_objects']}, "
          f"总字节={stats['total_bytes']}, cost_model={cost_model} ===")

    # 全局冷启动下界（compulsory）
    gfirst = {}
    for _t, oid, sz, _e in trace:
        gfirst.setdefault(oid, sz)
    gfloor_miss = len(gfirst)
    gfloor_cost = float(sum(gfirst.values()))

    grid_cost = {}   # (pct, k) -> cost/全局
    grid_miss = {}   # (pct, k) -> miss/全局
    for pct in pcts:
        cap = int(round(stats["total_bytes"] * pct / 100))
        g = run_one(trace, cap, cost_model)
        belady = get_algorithm("belady")
        bres = ContentSimulator(belady, capacity=cap,
                                cost_model=cost_model).run(trace)
        print(f"\n── [{pct}%] S={cap} （参照：Belady 全局 miss={bres.misses}；"
              f"全局冷启动下界 miss={gfloor_miss} / cost={gfloor_cost:.0f}）──")
        print(f"{'方案':<10}{'段长':>9}{'miss':>8}{'cost(字节)':>13}"
              f"{'evict':>8}{'rounds':>8}{'cost/全局':>10}{'miss/全局':>10}"
              f"{'miss/下界':>10}")
        print(f"{'全局':<10}{'-':>9}{g['misses']:>8}{g['cost']:>13.0f}"
              f"{g['evictions']:>8}{g['rounds']:>8}{1.0:>10.4f}{1.0:>10.4f}"
              f"{g['misses'] / gfloor_miss:>10.4f}")
        for k in segments:
            s = run_segmented(trace, cap, k, cost_model)
            grid_cost[(pct, k)] = s["cost"] / g["cost"] if g["cost"] else 1.0
            grid_miss[(pct, k)] = s["misses"] / g["misses"] if g["misses"] else 1.0
            print(f"{'分段k=' + str(k):<10}{_seg_len_txt(s):>9}"
                  f"{s['misses']:>8}{s['cost']:>13.0f}{s['evictions']:>8}"
                  f"{s['rounds']:>8}{grid_cost[(pct, k)]:>10.4f}"
                  f"{grid_miss[(pct, k)]:>10.4f}"
                  f"{s['misses'] / s['floor_miss']:>10.4f}")

    # 汇总网格（Markdown，便于粘贴报告）
    for title, grid in (("cost/全局", grid_cost), ("miss/全局", grid_miss)):
        print(f"\n── 汇总网格：{title} ──")
        head = "| 缓存比例 | " + " | ".join(f"k={k}" for k in segments) + " |"
        print(head)
        print("|---|" + "---|" * len(segments))
        for pct in pcts:
            cells = " | ".join(f"{grid[(pct, k)]:.4f}" for k in segments)
            print(f"| {pct}% | {cells} |")


def main():
    ap = argparse.ArgumentParser(
        description="分段 LocalRatioCaching（有界离线视野）vs 全局解对比实验")
    ap.add_argument("datasets", nargs="*",
                    help="trace 文件名或路径（默认 twitter29_test10.csv）")
    ap.add_argument("--pcts", type=int, nargs="+", default=[1, 2, 5, 10],
                    help="缓存容量占总字节百分比（默认 1 2 5 10）")
    ap.add_argument("--segments", type=int, nargs="+", default=[10, 20, 50, 100],
                    help="分段数 k（默认 10 20 50 100）")
    ap.add_argument("--cost-model", default="bit",
                    choices=["bit", "fault", "general"], help="成本模型（默认 bit）")
    args = ap.parse_args()

    names = args.datasets or ["twitter29_test10.csv"]
    for name in names:
        path = Path(name)
        if not path.exists():
            alt = TRACE_DIR / name
            if alt.exists():
                path = alt
            else:
                print(f"跳过 {name}: 文件不存在")
                continue
        run_dataset(path, args.pcts, args.segments, args.cost_model)


if __name__ == "__main__":
    main()
