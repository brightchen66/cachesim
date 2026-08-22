r"""General Caching 4-approximation -- Bar-Noy et al. (STOC 2000) §4 + §4.1.

论文《A Unified Approach to Approximating Resource Allocation and Scheduling》将
**general caching**（缓存大小固定、页面有任意 size 与任意 reload cost）归约为
**loss minimization** 问题，并用 local-ratio 技术给出 **4-近似** 离线算法。

================================================================================
问题（§4.1）
================================================================================
- 缓存大小 S>0；页面 j 有 size 0<s(j)<S 与 reload cost c(j)>0。
- 给定请求序列 r(1),...,r(n)（时间 i 的请求为 r(i)）。
- 替换调度：每个时刻 i 页面 r(i) 必在缓存中，且缓存内页面 size 之和 ≤ S。
- 代价 = Σ c(r(i))，对所有 i 使得 r(i) 在时刻 i−1 不在缓存中（即重新载入的代价）。
- 目标：最小化总 reload cost（离线，已知全部请求）。

================================================================================
归约到 loss minimization（§4.1，复述 Albers et al. [1]）
================================================================================
关键观察：若某页面在两次连续请求之间被驱逐过，不妨“请求后立即驱逐、下次请求时再取回”。
故只需考虑两类调度：两次连续请求之间，页面要么**全程在缓存**，要么**全程不在**。

对页面 p 在时刻 j<i 的两次连续请求，构造一个**活动实例**：
  - 时间区间 [j+1, i-1]（j+1 到 i-1 之间的中间时刻，p 可能被保留）；
  - width = s(p)（p 占用的缓存空间）；
  - penalty = c(p)（若不保留 p，则在 i 重新载入需付的代价）。
“调度（选中）该实例” = 在 (j,i) 期间保留 p，则 i 时刻命中、省下 c(p)；
“未调度” = 驱逐并重新载入，付 penalty c(p)。故**最小化未调度 penalty = 最小化 reload cost**。

死页扩展（超出论文归约的部分）：论文只为两次真实请求间的驻留建实例，末次请求后无后续
请求的页面不建实例（其驻留不影响代价）。本实现按用户指定的模型，假设死页的“下次请求”
为 n+1，为末次请求 j（j≤n-1）建立**虚拟实例** [j+1,n]：未调度 <=> 末次请求后驱逐；
被调度 <=> 保留至 trace 结尾。虚拟实例参与 local-ratio 优化（与真实实例竞争空间与
penalty），但其 penalty 仅存在于模型（n+1 不在真实请求序列内），**不计入真实 reload
cost**；代价是算法的 4-近似保证只对“真实+虚拟”的扩展目标成立，对真实 cost 不再有
直接的 4-近似保证（真实 cost ≤ 扩展代价 ≤ 4·OPT(扩展)）。

资源宽度 Width(t) = S − s(r(t))（时刻 t 请求页 r(t) 必在缓存，留给“保留页”的空间）。
仅当 i−j≥2（存在至少一个中间时刻）才创建实例；i−j==1（连续请求）页面必然保留（命中）；
首次请求 / 不可缓存页（s>S）为强制未命中（固定代价，不参与优化）。

================================================================================
Loss minimization local-ratio 算法（§4，4-近似）
================================================================================
每个活动仅单实例。维护各实例当前 penalty π(I) 与“是否存活”。迭代：

1. 找 t* = argmax Δ(t)，其中 Δ(t) = Σ_{存活 I 覆盖 t} width(I) − Width(t)（过载量）。
   若 max Δ ≤ 0：当前所有存活实例**共同可行**（处处 W_live ≤ Width），全部调度，终止。
2. Δ* = Δ(t*)；Z(t*) = 覆盖 t* 的存活实例集合；
   p = min_{I∈Z} π(I) / min(Δ*, width(I))  （使某个 π 恰好降到 0 的最大标量）。
3. 对 I∈Z：π(I) -= p·min(Δ*, width(I))。penalty 降到 ≤0 的实例**删除**（入栈），
   从存活集中移除（Δ 在其区间上下降 width(I)）。其余存活实例 penalty 降低但仍在。
阶段 2（贪心补回）：存活到终止的实例全部调度（先占用 slack）；再按栈 LIFO 弹出，
   若加入后仍可行（区间内最小 slack ≥ width）则调度。

**4-近似**：任何**极大**可行调度对 (M, p1) 都是 4-近似（Lemma 4.2：不稳定点 t 的
未调度 penalty ≤ 2pΔ*；t* 两侧最近不稳定点合计 ≤ 4pΔ* = 最优 p1-penalty 下界）。
local-ratio 定理把该界传递到原 penalty。

================================================================================
实现说明
================================================================================
- 离线、自带顶层 run() 循环（计算完整调度），不兼容 select_victim 接口；
  经 LocalRatioCachingSimulator 路由（schedule_based 标志）。
- 离线但**非最优**（4-近似），故不置 offline=True（否则被误设竞争比=1.0），
  而是与 Belady 基线对比计算两个竞争比。
- 数据结构：
  * 值域线段树（range-add / range-min / range-max+argmax）：阶段 1 维护 Δ(t) 取全局最大，
    阶段 2 维护 slack(t)=Width(t)−已调度宽度 取区间最小判可行性。
  * 区间线段树（按区间规范分解存实例 id 集合）：O(log n + |Z|) 查询覆盖 t* 的存活实例。
- cost 模型：bit/general -> c(j)=size(j)（优化字节代价）；fault -> c(j)=1（优化未命中次数）。
- evictions 口径（调度->驱逐的等价转换，两类驱逐均计数）：
  * 真实实例未调度：页面在上次请求后立即驱逐、下次请求重载（真实 reload 事件）；
  * 虚拟实例（死页）未调度：末次请求后驱逐（模型事件）；被调度则保留至结尾、不驱逐
    （末次请求恰在 n 的页面与零大小页同样保留至结尾）。
  即 evictions = 未调度真实实例数 + 未调度虚拟实例数。
"""

from typing import Dict, Iterable, List, Optional, Tuple

from cache_sim.algorithms.registry import register
from cache_sim.core.models import SimulationResult


# =====================================================================================
# 线段树：range-add + range-min/range-max(+argmax)，1-indexed 坐标 [1, n]
# =====================================================================================
class _SegTree:
    """支持区间加、区间最小、区间最大（带 argmax）的线段树（lazy propagation）。

    每个节点维护 mn / mx（该区间真值，已含本节点 lazy）与 lazy（待下传给子节点的增量）。
    """

    __slots__ = ("n", "mn", "mx", "lazy")

    def __init__(self, n: int, init_vals: List[float]):
        self.n = n
        size = 4 * n + 4
        self.mn = [0.0] * size
        self.mx = [0.0] * size
        self.lazy = [0.0] * size
        if n > 0:
            self._build(1, 1, n, init_vals)

    def _build(self, node: int, nl: int, nr: int, init_vals: List[float]) -> None:
        if nl == nr:
            v = init_vals[nl]
            self.mn[node] = v
            self.mx[node] = v
            return
        mid = (nl + nr) >> 1
        self._build(node << 1, nl, mid, init_vals)
        self._build((node << 1) | 1, mid + 1, nr, init_vals)
        self._pull(node)

    def _pull(self, node: int) -> None:
        lc, rc = node << 1, (node << 1) | 1
        self.mn[node] = self.mn[lc] if self.mn[lc] < self.mn[rc] else self.mn[rc]
        self.mx[node] = self.mx[lc] if self.mx[lc] > self.mx[rc] else self.mx[rc]

    def _apply(self, node: int, v: float) -> None:
        self.mn[node] += v
        self.mx[node] += v
        self.lazy[node] += v

    def _push(self, node: int) -> None:
        lz = self.lazy[node]
        if lz != 0.0:
            self._apply(node << 1, lz)
            self._apply((node << 1) | 1, lz)
            self.lazy[node] = 0.0

    def range_add(self, l: int, r: int, v: float) -> None:
        if l > r or self.n == 0:
            return
        self._add(1, 1, self.n, l, r, v)

    def _add(self, node: int, nl: int, nr: int, l: int, r: int, v: float) -> None:
        if r < nl or nr < l:
            return
        if l <= nl and nr <= r:
            self._apply(node, v)
            return
        self._push(node)
        mid = (nl + nr) >> 1
        self._add(node << 1, nl, mid, l, r, v)
        self._add((node << 1) | 1, mid + 1, nr, l, r, v)
        self._pull(node)

    def range_min(self, l: int, r: int) -> float:
        if l > r or self.n == 0:
            return float("inf")
        return self._qmin(1, 1, self.n, l, r)

    def _qmin(self, node: int, nl: int, nr: int, l: int, r: int) -> float:
        if r < nl or nr < l:
            return float("inf")
        if l <= nl and nr <= r:
            return self.mn[node]
        self._push(node)
        mid = (nl + nr) >> 1
        a = self._qmin(node << 1, nl, mid, l, r)
        b = self._qmin((node << 1) | 1, mid + 1, nr, l, r)
        return a if a < b else b

    def range_max_argmax(self, l: int, r: int) -> Tuple[float, int]:
        """返回 [l,r] 内 (最大值, 任意一个取到最大值的下标)。"""
        if l > r or self.n == 0:
            return (float("-inf"), -1)
        return self._qmax(1, 1, self.n, l, r)

    def _qmax(self, node: int, nl: int, nr: int, l: int, r: int) -> Tuple[float, int]:
        if r < nl or nr < l:
            return (float("-inf"), -1)
        if l <= nl and nr <= r:
            return (self.mx[node], self._argmax_leaf(node, nl, nr))
        self._push(node)
        mid = (nl + nr) >> 1
        a = self._qmax(node << 1, nl, mid, l, r)
        b = self._qmax((node << 1) | 1, mid + 1, nr, l, r)
        return a if a[0] >= b[0] else b

    def _argmax_leaf(self, node: int, nl: int, nr: int) -> int:
        """下钻到该节点区间内一个取到 mx[node] 的叶子下标。"""
        while nl != nr:
            self._push(node)
            mid = (nl + nr) >> 1
            lc = node << 1
            if self.mx[lc] >= self.mx[(node << 1) | 1]:
                node, nr = lc, mid
            else:
                node, nl = (node << 1) | 1, mid + 1
        return nl


# =====================================================================================
# 区间线段树：按区间规范分解存实例 id，点查询返回覆盖该点的全部实例
# =====================================================================================
class _IntervalTree:
    """每个实例按 [a,b] 的规范分解存入 O(log n) 个节点；点查询收集根→叶路径上的全部实例。

    用 dict 按需建节点集合，避免 4n 个空集合的内存开销（n 大时关键）。
    """

    __slots__ = ("n", "nodes")

    def __init__(self, n: int):
        self.n = n
        self.nodes: Dict[int, set] = {}

    def insert(self, iid: int, l: int, r: int) -> None:
        self._ins(1, 1, self.n, l, r, iid)

    def _ins(self, node: int, nl: int, nr: int, l: int, r: int, iid: int) -> None:
        if l <= nl and nr <= r:
            self.nodes.setdefault(node, set()).add(iid)
            return
        mid = (nl + nr) >> 1
        if l <= mid:
            self._ins(node << 1, nl, mid, l, r, iid)
        if r > mid:
            self._ins((node << 1) | 1, mid + 1, nr, l, r, iid)

    def remove(self, iid: int, l: int, r: int) -> None:
        self._rem(1, 1, self.n, l, r, iid)

    def _rem(self, node: int, nl: int, nr: int, l: int, r: int, iid: int) -> None:
        s = self.nodes.get(node)
        if s is not None:
            s.discard(iid)
            if not s:
                del self.nodes[node]
        if l <= nl and nr <= r:
            return
        mid = (nl + nr) >> 1
        if l <= mid:
            self._rem(node << 1, nl, mid, l, r, iid)
        if r > mid:
            self._rem((node << 1) | 1, mid + 1, nr, l, r, iid)

    def query(self, pt: int) -> List[int]:
        """覆盖点 pt 的全部实例 id（沿根→叶路径收集，无重复）。"""
        res: List[int] = []
        node, nl, nr = 1, 1, self.n
        while True:
            s = self.nodes.get(node)
            if s:
                res.extend(s)
            if nl == nr:
                break
            mid = (nl + nr) >> 1
            if pt <= mid:
                node, nr = node << 1, mid
            else:
                node, nl = (node << 1) | 1, mid + 1
        return res


# =====================================================================================
# General Caching 4-近似算法
# =====================================================================================
class LocalRatioCaching:
    """论文 §4.1 + §4 的 general caching 离线 4-近似算法。

    通过 :meth:`run` 消费 (time, id, size[, extra]) trace，返回 :class:`SimulationResult`。
    """

    # 标记：基于完整调度的算法，不兼容 select_victim 接口，走专用模拟器
    schedule_based = True
    # 理论近似比上界
    APPROX_RATIO = 4.0

    def __init__(self, eps: float = 1e-9, verbose: bool = False):
        self.eps = float(eps)
        # verbose=True 时逐轮打印 local-ratio 递归的最拥塞时间点 t* 等信息
        self.verbose = bool(verbose)

    def solve(self, trace: Iterable[Tuple], capacity: int,
              cost_model: str = "bit") -> Dict:
        """求解完整替换调度（供测试与复用）。

        返回字典含：
          - scheduled: 被调度（保留）的实例 id 集合；
          - iterations: 每轮递归记录（t* / Δ* / |Z| / p / 删除数与删除对象 size 范围）；
          - inst_a/inst_b/inst_w/inst_pen/inst_req/inst_virtual: 实例元数据；
          - inst_for_req: 请求 i -> 真实实例 iid；inst_for_final: 末次请求 i -> 虚拟实例 iid；
          - prev/nxt/req_size/page_size/width_at/cost_model/S/n: 调度上下文。
        """
        if capacity <= 0:
            raise ValueError(f"capacity 必须为正，得到 {capacity}")
        S = int(capacity)
        eps = self.eps
        items = list(trace)  # 离线：物化整条 trace
        n = len(items)

        # ---- 每页规范 size（首次出现），每请求的实际 size/cost ----
        page_size: Dict[int, int] = {}
        req_size = [0] * (n + 1)       # req_size[i] = 请求 i (1-indexed) 的 size
        for i in range(1, n + 1):
            sz = int(items[i - 1][2])
            req_size[i] = sz
            oid = items[i - 1][1]
            if oid not in page_size:
                page_size[oid] = sz

        def cost_of(oid: int) -> float:
            if cost_model == "fault":
                return 1.0
            return float(page_size.get(oid, 0))

        def cacheable(oid: int) -> bool:
            return page_size.get(oid, 0) <= S

        # ---- prev[i] = 页 r(i) 上一次请求的序号（1-indexed），0 表示无 ----
        prev = [0] * (n + 1)
        last: Dict[int, int] = {}
        for i in range(1, n + 1):
            oid = items[i - 1][1]
            prev[i] = last.get(oid, 0)
            last[oid] = i

        # ---- next[i] = 页 r(i) 下一次请求的序号（1-indexed），0 表示无（该请求后为死页）----
        nxt = [0] * (n + 1)
        last.clear()
        for i in range(n, 0, -1):
            oid = items[i - 1][1]
            nxt[i] = last.get(oid, 0)
            last[oid] = i

        # ---- Width(t) for t=1..n：r(t) 可缓存则 S-s(r(t))，否则 S（不缓存，全空间留给保留页）----
        width_at = [0] * (n + 1)
        for t in range(1, n + 1):
            oid_t = items[t - 1][1]
            st = page_size.get(oid_t, 0)
            width_at[t] = S if st > S else S - st

        # ---- 构造活动实例 ----
        # 实例 iid：区间 [a,b]、width w、penalty π、对应请求序号（见 inst_req）。
        # 两类实例：
        #   * 真实实例：页面在两次真实连续请求 (j,i) 之间（i-j>=2），区间 [j+1,i-1]；
        #     未调度 <=> 请求 j 后驱逐、i 时重载（penalty = 真实 reload cost）。
        #   * 虚拟实例（死页扩展）：末次请求 j 后无后续请求（nxt[j]==0），假设
        #     “下次请求”为 n+1，区间 [j+1,n]（要求 j<=n-1，j==n 时区间空、必然保留）；
        #     未调度 <=> 末次请求后驱逐。其 penalty 只存在于模型（n+1 不在真实请求
        #     序列内），**不计入真实 reload cost**。
        inst_a: List[int] = []
        inst_b: List[int] = []
        inst_w: List[float] = []
        inst_pen: List[float] = []
        inst_req: List[int] = []        # 该实例对应的请求序号 i
        inst_virtual: List[bool] = []   # True = 虚拟实例（死页，假设下次请求 n+1）
        inst_for_req: Dict[int, int] = {}   # 请求 i -> 真实实例 iid（i-j>=2 且可缓存时）
        inst_for_final: Dict[int, int] = {} # 末次请求 i（nxt=0）-> 虚拟实例 iid
        for i in range(1, n + 1):
            oid = items[i - 1][1]
            j = prev[i]
            if j == 0 or not cacheable(oid):
                continue
            w = float(page_size[oid])
            if w <= 0:
                continue  # 零大小页：不占空间、必然可保留（真实请求命中、死页保留至结尾）
            if i - j >= 2:
                # 真实实例：区间内保留 <=> 请求 i 命中
                iid = len(inst_a)
                inst_a.append(j + 1)
                inst_b.append(i - 1)
                inst_w.append(w)
                inst_pen.append(cost_of(oid))   # penalty = reload cost c(p)
                inst_req.append(i)
                inst_virtual.append(False)
                inst_for_req[i] = iid
            if nxt[i] == 0 and i <= n - 1:
                # 虚拟实例：末次请求后保留至结尾（“下次请求”为虚拟的 n+1）
                iid = len(inst_a)
                inst_a.append(i + 1)
                inst_b.append(n)
                inst_w.append(w)
                inst_pen.append(cost_of(oid))   # 仅模型 penalty，不计入真实 cost
                inst_req.append(i)              # 对应末次请求（虚拟下次请求为 n+1）
                inst_virtual.append(True)
                inst_for_final[i] = iid

        m = len(inst_a)
        scheduled: set = set()
        iterations: List[Dict] = []
        if m > 0:
            scheduled, iterations = self._solve_schedule(
                n, S, eps, width_at,
                inst_a, inst_b, inst_w, inst_pen)

        return {
            "scheduled": scheduled,
            "iterations": iterations,
            "inst_a": inst_a, "inst_b": inst_b, "inst_w": inst_w,
            "inst_pen": inst_pen, "inst_req": inst_req,
            "inst_virtual": inst_virtual,
            "inst_for_req": inst_for_req, "inst_for_final": inst_for_final,
            "prev": prev, "nxt": nxt, "req_size": req_size, "page_size": page_size,
            "width_at": width_at, "cost_model": cost_model, "S": S, "n": n,
            "items": items, "cacheable": cacheable, "cost_of": cost_of,
        }

    def run(self, trace: Iterable[Tuple], capacity: int,
            cost_model: str = "bit", dataset_name: str = "") -> SimulationResult:
        """运行离线 general caching 4-近似。

        Args:
            trace: (time, id, size[, extra]) 记录可迭代对象。
            capacity: 缓存容量 S（字节）。
            cost_model: "bit"/"general" -> c(j)=size；"fault" -> c(j)=1。
            dataset_name: 数据集名（用于结果）。
        """
        sol = self.solve(trace, capacity, cost_model)
        S = sol["S"]
        n = sol["n"]
        items = sol["items"]
        scheduled = sol["scheduled"]
        iterations = sol["iterations"]
        inst_for_req = sol["inst_for_req"]
        inst_for_final = sol["inst_for_final"]
        inst_virtual = sol["inst_virtual"]
        inst_pen = sol["inst_pen"]
        prev = sol["prev"]
        nxt = sol["nxt"]
        page_size = sol["page_size"]
        req_size = sol["req_size"]
        cost_of = sol["cost_of"]
        cacheable = sol["cacheable"]

        hits = 0
        misses = 0
        byte_total = 0
        byte_hit = 0
        evictions = 0
        cost = 0.0

        # ---- 由调度统计命中/未命中/代价/驱逐 ----
        # 驱逐口径（与调度的等价转换一致，两类事件均计数）：
        #   (a) 真实实例未调度：页面在上次请求后立即被驱逐、本次重载 -- 一次驱逐；
        #   (b) 虚拟实例（死页，假设下次请求 n+1）未调度：末次请求后驱逐 -- 一次驱逐；
        #       被调度 / 末次请求恰在 n（区间空，必保留）/ 零大小页：保留至结尾，不驱逐。
        # 故 evictions = 未调度真实实例数 + 未调度虚拟实例数；
        # cost 仍为**真实** reload 代价（虚拟 penalty 不计入，n+1 不在请求序列内）。
        virtual_penalty = sum(
            inst_pen[iid] for iid in range(len(inst_virtual))
            if inst_virtual[iid] and iid not in scheduled)
        for i in range(1, n + 1):
            oid = items[i - 1][1]
            sz = req_size[i]
            byte_total += sz
            j = prev[i]
            c = cost_of(oid)
            cacheable_i = cacheable(oid)
            if j == 0:
                # 首次请求：强制载入（不可缓存页不入缓存，仅计 miss）
                misses += 1
                cost += c
            elif not cacheable_i:
                # 不可缓存页：无法驻留，每次请求都未命中
                misses += 1
                cost += c
            elif i - j < 2:
                # 连续请求：页面自上次请求仍保留，命中
                hits += 1
                byte_hit += sz
            elif page_size.get(oid, 0) <= 0:
                # 零大小页：不占空间、必然可保留（solve 跳过建实例），命中
                hits += 1
                byte_hit += sz
            else:
                iid = inst_for_req.get(i)
                if iid is not None and iid in scheduled:
                    hits += 1
                    byte_hit += sz
                else:
                    misses += 1
                    cost += c
                    evictions += 1   # (a) 实例未调度：请求 j 后被驱逐，本次重载
            if nxt[i] == 0 and cacheable_i:
                # (b) 死页：末次请求后是否驱逐由虚拟实例（假设下次请求 n+1）决定
                viid = inst_for_final.get(i)
                if viid is not None and viid not in scheduled:
                    evictions += 1       # 虚拟实例未调度：末次请求后驱逐
                # 否则（虚拟实例被调度 / 末次请求恰在 n / 零大小页）：
                # 保留至 trace 结尾，不驱逐

        return SimulationResult(
            cache_type="content",
            algorithm="LocalRatioCaching",
            dataset=dataset_name,
            config={
                "capacity": S,
                "cost_model": cost_model,
                "offline": True,        # 算法本身是离线的
                "approx_ratio": self.APPROX_RATIO,
            },
            total_requests=n,
            hits=hits,
            misses=misses,
            byte_total=byte_total,
            byte_hit=byte_hit,
            evictions=evictions,
            competitive_ratio=None,     # 由调用方对比 Belady 回填
            byte_competitive_ratio=None,
            extra={
                "cost": cost,            # 算法优化的总 reload cost
                "num_instances": len(inst_for_req),
                "num_scheduled": len(scheduled),
                "approx_ratio_bound": self.APPROX_RATIO,
                "num_rounds": len(iterations),          # local-ratio 递归轮数
                # 每轮递归记录：round / t_star / delta_star / num_z / p /
                # num_deleted / del_size_min/max/sum（删除实例的对象 size 范围）
                "iterations": iterations,
                # 虚拟实例（死页，假设下次请求 n+1）统计；其 penalty 为模型代价，不计入 cost
                "num_virtual_instances": sum(inst_virtual),
                "virtual_penalty": virtual_penalty,
            },
        )

    # ---------------------------------------------------------------------------------
    # local-ratio 调度求解：返回 (被调度实例 id 集合, 每轮递归记录列表)
    # ---------------------------------------------------------------------------------
    def _solve_schedule(self, n: int, S: int, eps: float, width_at: List[int],
                        inst_a: List[int], inst_b: List[int],
                        inst_w: List[float], inst_pen: List[float]) -> Tuple[set, List[Dict]]:
        m = len(inst_a)
        # 值域线段树：Δ(t) = W_live(t) - Width(t)。初始 = -Width(t)，每实例 +w 于其区间。
        init_delta = [0.0] * (n + 1)
        for t in range(1, n + 1):
            init_delta[t] = -float(width_at[t])
        val_tree = _SegTree(n, init_delta)
        for iid in range(m):
            val_tree.range_add(inst_a[iid], inst_b[iid], inst_w[iid])

        # 区间线段树：存存活实例 id，O(log n + |Z|) 查询覆盖 t* 的实例
        itree = _IntervalTree(n)
        for iid in range(m):
            itree.insert(iid, inst_a[iid], inst_b[iid])

        pen = list(inst_pen)        # 当前 penalty（可变副本）
        alive = [True] * m
        stack: List[int] = []       # 删除顺序（阶段 2 LIFO 弹出）
        # 每轮递归记录：round / t_star（最拥塞时间点）/ delta_star（过载量）/
        # num_z（覆盖 t* 的存活实例数）/ p（局部 penalty 比例因子）/
        # num_deleted（本轮删除数）/ del_size_min/max/sum（删除实例的对象 size 范围）
        iterations: List[Dict] = []

        # ---- 阶段 1：local-ratio 迭代 ----
        # 在 [2, n] 寻找最大过载点（真实实例区间 ⊆ [2,n-1]，虚拟实例延伸到 n）。
        lo_q, hi_q = 2, n
        if lo_q > hi_q:
            # 无中间时刻：所有实例“存活到终止”，全部调度（若可行）
            return (self._finalize(m, n, eps, width_at, inst_a, inst_b, inst_w,
                                   alive, stack), iterations)

        while True:
            val, t_star = val_tree.range_max_argmax(lo_q, hi_q)
            if val <= eps:
                break  # max Δ ≤ 0：所有存活实例共同可行，终止
            delta_star = val
            Z = itree.query(t_star)
            # 计算 p = min_{I∈Z} π(I) / min(Δ*, w(I))，跳过 min(Δ*,w)=0 的实例
            best_p = float("inf")
            for iid in Z:
                if not alive[iid]:
                    continue
                mw = delta_star if delta_star < inst_w[iid] else inst_w[iid]
                if mw <= 0:
                    continue
                r = pen[iid] / mw
                if r < best_p:
                    best_p = r
            if best_p == float("inf") or best_p <= 0:
                # 无可降低的实例（理论上不应发生：Δ*>0 必有正宽实例）：安全终止
                break
            # 降低 Z 中存活实例的 penalty，删除降到 ≤0 的
            to_delete: List[int] = []
            for iid in Z:
                if not alive[iid]:
                    continue
                mw = delta_star if delta_star < inst_w[iid] else inst_w[iid]
                if mw <= 0:
                    continue
                pen[iid] -= best_p * mw
                if pen[iid] <= eps:
                    to_delete.append(iid)
            round_deleted = list(to_delete)
            for iid in to_delete:
                alive[iid] = False
                val_tree.range_add(inst_a[iid], inst_b[iid], -inst_w[iid])
                itree.remove(iid, inst_a[iid], inst_b[iid])
                stack.append(iid)
            if not to_delete:
                # penalty 全降但无一到 0（浮点）：强制删除当前 penalty 最小者避免死循环
                mn, mn_id = float("inf"), -1
                for iid in Z:
                    if alive[iid] and pen[iid] < mn:
                        mn, mn_id = pen[iid], iid
                if mn_id < 0:
                    break
                alive[mn_id] = False
                val_tree.range_add(inst_a[mn_id], inst_b[mn_id], -inst_w[mn_id])
                itree.remove(mn_id, inst_a[mn_id], inst_b[mn_id])
                stack.append(mn_id)
                round_deleted.append(mn_id)
            # ---- 记录本轮递归（论文中每次递归调用对应此处一轮迭代）----
            d_sizes = [inst_w[iid] for iid in round_deleted]
            iterations.append({
                "round": len(iterations) + 1,
                "t_star": t_star,          # 最拥塞时间点：argmax Δ(t)
                "delta_star": delta_star,  # 过载量 A*
                "num_z": len(Z),           # 覆盖 t* 的存活实例数 |Z(t*)|
                "p": best_p,               # 使某 penalty 恰好降为 0 的比例因子
                "num_deleted": len(round_deleted),
                # 本轮删除实例对应对象的 size 范围与总和（width = s(p)）
                "del_size_min": min(d_sizes) if d_sizes else None,
                "del_size_max": max(d_sizes) if d_sizes else None,
                "del_size_sum": sum(d_sizes) if d_sizes else 0.0,
            })
            if self.verbose:
                rng_txt = (f"（对象大小 {min(d_sizes):g}–{max(d_sizes):g} 字节，"
                           f"合计 {sum(d_sizes):g}）") if d_sizes else ""
                print(f"[local-ratio 第{len(iterations)}轮] 最拥塞时间点 t*={t_star}, "
                      f"Δ*={delta_star:g}, |Z(t*)|={len(Z)}, p={best_p:g}, "
                      f"删除 {len(round_deleted)} 个实例{rng_txt}")

        return (self._finalize(m, n, eps, width_at, inst_a, inst_b, inst_w,
                               alive, stack), iterations)

    def _finalize(self, m: int, n: int, eps: float, width_at: List[int],
                  inst_a: List[int], inst_b: List[int], inst_w: List[float],
                  alive: List[bool], stack: List[int]) -> set:
        """阶段 2：存活到终止的实例全部调度，再 LIFO 贪心补回栈中实例。"""
        # slack(t) = Width(t) - 已调度宽度。初始 = Width(t)。
        init_slack = [0.0] * (n + 1)
        for t in range(1, n + 1):
            init_slack[t] = float(width_at[t])
        slack_tree = _SegTree(n, init_slack)
        scheduled: set = set()

        # 存活到终止的实例全部调度（max Δ≤0 保证共同可行）
        for iid in range(m):
            if alive[iid]:
                scheduled.add(iid)
                slack_tree.range_add(inst_a[iid], inst_b[iid], -inst_w[iid])

        # LIFO 弹出，若区间内最小 slack ≥ w 则调度
        for idx in range(len(stack) - 1, -1, -1):
            iid = stack[idx]
            mn = slack_tree.range_min(inst_a[iid], inst_b[iid])
            if mn >= inst_w[iid] - eps:
                scheduled.add(iid)
                slack_tree.range_add(inst_a[iid], inst_b[iid], -inst_w[iid])
        return scheduled


# =====================================================================================
# 暴力最优（仅用于小规模验证 4-近似）：DP over (时刻, 缓存状态)
# =====================================================================================
def optimal_caching_bruteforce(trace: Iterable[Tuple], capacity: int,
                               cost_model: str = "bit") -> float:
    """小规模 general caching 的精确最优总 reload cost（验证用，状态空间 2^页面数）。

    DP：dp[t][state] = 服务前 t 个请求、时刻 t 缓存状态为 state（含 r(t)、size 和≤S）的最小代价。
    仅适用于页面数很少（≤ ~14）的实例；用于单元测试验证 4-近似。
    """
    items = list(trace)
    n = len(items)
    if n == 0:
        return 0.0
    S = int(capacity)
    page_size: Dict[int, int] = {}
    for it in items:
        oid = it[1]
        if oid not in page_size:
            page_size[oid] = int(it[2])

    def cost_of(oid: int) -> float:
        return 1.0 if cost_model == "fault" else float(page_size.get(oid, 0))

    def size_of(oid: int) -> int:
        return page_size.get(oid, 0)

    # 所有可能入缓存的页面（size ≤ S）
    cacheable_pages = [p for p in page_size if page_size[p] <= S]
    # 枚举所有合法缓存状态（子集，size 和 ≤ S）—— 仅页面数少时可行
    all_subsets: List[frozenset] = []
    pidx = {p: i for i, p in enumerate(cacheable_pages)}
    P = len(cacheable_pages)

    def enum(i: int, cur: frozenset, cur_sz: int):
        if i == P:
            all_subsets.append(cur)
            return
        p = cacheable_pages[i]
        enum(i + 1, cur, cur_sz)
        if cur_sz + page_size[p] <= S:
            enum(i + 1, cur | {p}, cur_sz + page_size[p])

    if P <= 16:
        enum(0, frozenset(), 0)
    else:
        raise ValueError("页面数过多，无法暴力求最优")

    INF = float("inf")
    # 时刻 0 缓存为空
    dp: Dict[frozenset, float] = {frozenset(): 0.0}
    for t in range(1, n + 1):
        rt = items[t - 1][1]
        srt = size_of(rt)
        ndp: Dict[frozenset, float] = {}
        for prev_state, prev_cost in dp.items():
            hit = rt in prev_state
            add_cost = 0.0 if hit else cost_of(rt)
            base = prev_cost + add_cost
            # 若 rt 不可缓存，则 t 时刻状态 = prev_state（不接纳 rt）
            if srt > S:
                ns = prev_state
                if ns not in ndp or base < ndp[ns]:
                    ndp[ns] = base
                continue
            # 枚举下一状态：prev_state ∪ {rt} 的子集且含 rt、size 和 ≤ S
            cand = set(prev_state)
            cand.add(rt)
            cand_sz = sum(size_of(p) for p in cand)
            # 枚举 cand 的合法子集（含 rt）
            cand_list = list(cand)
            Q = len(cand_list)
            for mask in range(1 << Q):
                ns_set = set()
                ns_sz = 0
                ok = False
                for k in range(Q):
                    if mask & (1 << k):
                        p = cand_list[k]
                        ns_set.add(p)
                        ns_sz += size_of(p)
                        if p == rt:
                            ok = True
                if not ok or ns_sz > S:
                    continue
                ns = frozenset(ns_set)
                if ns not in ndp or base < ndp[ns]:
                    ndp[ns] = base
        dp = ndp
    return min(dp.values()) if dp else 0.0


@register("local_ratio_caching")
def _make_local_ratio_caching(eps: float = 1e-9,
                              verbose: bool = False) -> LocalRatioCaching:
    return LocalRatioCaching(eps=eps, verbose=verbose)
