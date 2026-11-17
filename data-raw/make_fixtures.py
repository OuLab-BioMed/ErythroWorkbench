"""make_fixtures.py — 从库内真实数据集生成跨语言对拍 fixture。

为什么用真实矩阵的子集而不是构造数据：秩是在**全基因背景**下算的，
随手造 47 个 marker 的矩阵会让 rank 全落在 1..47，与真实情形（HBB 排在
两万个基因的最前面）完全不同，对拍通过也说明不了什么。

做法：保留全部 47 个谱系 marker，再随机抽 N 个背景基因（固定 seed），
两边都从这一份矩阵跑到最终分数。Python 侧的期望输出一并落盘。

用法（在包根目录执行）：
    python data-raw/make_fixtures.py
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parent.parent
REPO = PKG.parent.parent                     # EryWB/
sys.path.insert(0, str(REPO / "app"))
sys.path.insert(0, str(REPO))

from backend.engines import toolkit_core as tc          # noqa: E402
from backend.engines import sample_types as st          # noqa: E402

DATA = REPO / "app" / "data"
BACKEND_DATA = REPO / "app" / "backend" / "data"
OUT = PKG / "tests" / "testthat" / "fixtures"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20261117
N_BACKGROUND = 3000


def load_matrix(dataset: str) -> pd.DataFrame:
    # 读库内矩阵，兼容三种首列写法：probe_id(+gene_symbol) / probe_id / metabolite。
    # 样本列名常带 GEO 导出的三引号（GSM 名被三个双引号包住），不 strip 掉会让
    # 注释对齐失败，且报错信息完全看不出原因。
    p = DATA / "matrices" / dataset / "matrix.tsv"
    m = pd.read_csv(p, sep="\t")
    if "gene_symbol" in m.columns:
        ids = m["gene_symbol"].astype(str).str.strip().str.strip('"')
        m = m.drop(columns=["gene_symbol"])
    else:
        ids = m[m.columns[0]].astype(str).str.strip().str.strip('"')
    m = m.drop(columns=[m.columns[0]])
    m.columns = [str(c).strip().strip('"') for c in m.columns]
    m.index = ids
    x = m.apply(pd.to_numeric, errors="coerce")
    return x.groupby(level=0).mean()


def write_subset(dataset: str, markers, name: str) -> pd.DataFrame:
    x = load_matrix(dataset)
    keep = [g for g in markers if g in x.index]
    rng = np.random.default_rng(SEED)
    bg = [g for g in x.index if g not in set(keep)]
    bg = list(rng.choice(bg, size=min(N_BACKGROUND, len(bg)), replace=False))
    sub = x.loc[keep + bg]
    out = sub.copy()
    out.insert(0, "gene_symbol", out.index)
    out.to_csv(OUT / name, sep="\t", index=False, compression="gzip",
               float_format="%.6g")
    print(f"  {name}: {sub.shape[0]} features x {sub.shape[1]} samples "
          f"({len(keep)} markers + {len(bg)} background)")
    return sub


def purity_expected(x: pd.DataFrame, layer: str = "转录组") -> pd.DataFrame:
    """按引擎口径算一遍，作为 R 侧的对拍目标。"""
    refs = json.load(open(BACKEND_DATA / "layer_references.json", encoding="utf-8"))
    ery, leu = refs["markers"]["ery"], refs["markers"]["leu"]
    ## 与引擎 run_purity 同一条入口：先归一到标准符号空间再算
    names, X = tc.to_space(x, "gene")
    L, _ = tc.as_log2(X)
    names = list(names)
    idx = {n: i for i, n in enumerate(names)}
    R = np.apply_along_axis(
        lambda v: tc.rank_profile(v) / max(1.0, (np.isfinite(v).sum() - 1)), 0, L)
    trefs = json.load(open(BACKEND_DATA / "toolkit_refs.json", encoding="utf-8"))["purity"]
    k = float(trefs["identity"]["k"])
    c = float(trefs["identity"]["center"])
    cp = trefs.get("cell_prior") or {}
    auc_max = float(cp.get("auc_leu_endpoint_max") or 9)
    rows = []
    for j, s in enumerate(x.columns):
        e_full = np.array([R[idx[g], j] if g in idx else np.nan for g in ery], float)
        l_full = np.array([R[idx[g], j] if g in idx else np.nan for g in leu], float)
        cons = st.layer_consistency(np.concatenate([e_full, l_full]), layer)
        er, lr = np.nanmean(e_full), np.nanmean(l_full)
        d = er - lr
        a = e_full[np.isfinite(e_full)]
        b = l_full[np.isfinite(l_full)]
        auc = float(np.mean(a[:, None] > b[None, :]) + 0.5 * np.mean(a[:, None] == b[None, :]))
        strong = bool(len(b) >= 8 and auc <= auc_max + 0.005)
        h = st.headline_score(cons.get("score"), layer, "high", True, False)
        rows.append({
            "sample": str(s), "ery_rank": round(er, 3), "leu_rank": round(lr, 3),
            "delta": round(d, 3), "auc": round(auc, 3),
            "rna_identity_score": round(float(1 / (1 + np.exp(-np.clip(k * (d - c), -30, 30)))), 3),
            "ery_consistency": cons.get("score"),
            "consistency_rho": cons.get("rho"),
            "n_markers": cons.get("n_markers"),
            "conflict": int(bool(auc < 0.5) and not strong),
            "headline": h["score"], "prior": h["prior"], "ceiling": h["ceiling"],
        })
    return pd.DataFrame(rows)


def projection_expected(x: pd.DataFrame) -> pd.DataFrame:
    """投影标签转移的引擎输出，作为 R 侧对拍目标。"""
    from backend.engines import reference_projection as rp

    a = rp.atlas()
    samples = a["samples"]
    ref = np.array([i for i, s in enumerate(samples) if s["labeled"]])
    gate = float(a["ood_gate"]["rank_p5"])
    mgate = float(json.load(open(BACKEND_DATA / "reference_projection_validation.json",
                                 encoding="utf-8"))["spaces"]["rank"]["margin_gate_p10"])
    markers = a["markers"]
    names, X = tc.to_space(x, "gene")
    L, _ = tc.as_log2(X)
    names = list(names)
    idx = {n: i for i, n in enumerate(names)}
    R = np.apply_along_axis(
        lambda v: tc.rank_profile(v) / max(1.0, (np.isfinite(v).sum() - 1)), 0, L)
    rows = []
    for j, col in enumerate(x.columns):
        q = np.array([R[idx[g], j] if g in idx else np.nan for g in markers], float)
        sims = np.full(len(ref), np.nan, float)
        for k, i in enumerate(ref):
            v = rp._sim(q, rp._vec(samples[int(i)], "rank"), "rank")
            sims[k] = np.nan if v is None else v
        full = np.full(len(samples), np.nan, float)
        full[ref] = sims
        lab, prob, best, ent, why = rp._vote(full, ref, samples, gate, mgate)
        rows.append({"sample": str(col), "label": lab or "", "probability": prob,
                     "best_rho": best, "entropy": ent, "abstain": why or "",
                     "n_markers": int(np.isfinite(q).sum())})
    return pd.DataFrame(rows)


def clock_expected(dataset: str, space: str) -> pd.DataFrame:
    """时钟模块的引擎输出（预训练投影路径）。"""
    from backend.engines import toolkit_tools as tt

    x = load_matrix(dataset)
    names, X = tc.to_space(x, space)          # 与引擎同一套名空间标准化
    finite = np.isfinite(X)
    # _prepare 不做 as_log2，且用**行中位数**填补缺失；少这一步会把"整行缺失"
    # 的样本当成"测到了、值等于中位数"，给出一个看起来正常却没有依据的天数
    X = np.where(finite, X, np.nanmedian(X, axis=1, keepdims=True))
    prep = {"names": names, "X": X, "space": space,
            "mat": x, "finite": finite, "ann": pd.DataFrame()}
    out = tt._clock_reference(prep)
    if not out.get("ok"):
        print("  clock reference failed:", out.get("reason"))
        return pd.DataFrame()
    preds = out["predictions"]
    return pd.DataFrame([{
        "sample": r["sample"], "day": r.get("day"), "lo": r.get("lo"),
        "hi": r.get("hi"), "coverage": r.get("coverage"),
        "raw": r.get("raw"), "clipped": r.get("clipped"),
    } for r in preds])


def module_clock_expected(dataset: str, space: str = "metabolite") -> pd.DataFrame:
    """跨层外推（模块桥）路径的引擎输出。"""
    from backend.engines import toolkit_tools as tt

    x = load_matrix(dataset)
    names, X = tc.to_space(x, space)
    X = np.where(np.isfinite(X), X, np.nanmedian(X, axis=1, keepdims=True))
    prep = {"names": names, "X": X, "space": space, "mat": x, "ann": pd.DataFrame()}
    out = tt._clock_module(prep)
    if not out.get("ok"):
        print("  module clock failed:", out.get("reason"))
        return pd.DataFrame()
    return pd.DataFrame([{
        "sample": r["sample"], "day": r.get("day"), "lo": r.get("lo"),
        "hi": r.get("hi"), "raw": r.get("raw"), "clipped": r.get("clipped"),
    } for r in out["predictions"]])


def pathway_expected(x: pd.DataFrame, mode: str = "within") -> pd.DataFrame:
    """通路指纹的引擎输出。"""
    names, X = tc.to_space(x, "gene")
    L, _ = tc.as_log2(X)
    Z = tc.within_sample_z(L) if mode == "within" else tc.zscore_rows(L)
    Z = np.where(np.isfinite(Z), Z, 0.0)
    members = tc.pathway_members("gene", scoring_only=True)
    pw = tc.pathway_kb()["pathways"]
    pos = {str(n).upper(): i for i, n in enumerate(names)}
    rows = []
    for pid, ms in members.items():
        gi = [pos[str(m).upper()] for m in ms if str(m).upper() in pos]
        if len(gi) < 2:
            continue
        sub = Z[gi]
        nm = pw[pw["pathway_id"] == pid]
        frac = [max((sub[:, j] > 0).sum(), (sub[:, j] < 0).sum()) / sub.shape[0]
                for j in range(sub.shape[1])]
        row = {"pathway_id": pid, "n_members": len(ms), "n_measured": len(gi),
               "consistency": round(float(np.median(frac)), 3)}
        for j, s in enumerate(x.columns):
            row[str(s)] = round(float(np.mean(sub[:, j])), 3)
        rows.append(row)
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    mxs = df[x.columns.astype(str).tolist()].abs().max(axis=1)
    df["_key"] = mxs * np.sqrt(df["n_measured"])
    df = df.sort_values("_key", ascending=False).drop(columns=["_key"])
    return df.reset_index(drop=True)


def attribution_expected(x: pd.DataFrame, n_perm: int = 200,
                         parity_samples=None, parity_n: int = 15):
    """条件谱匹配的引擎输出，外加可复用的置换矩阵。

    与 run_attribution() 同一条路径，只差 upload 那一层：身份用
    si.scope(..., 无注释)，所以 identity 全按样本名判定。

    置换序列是唯一无法跨语言对齐的东西（numpy PCG64 vs R 的 MT）。所以这里
    把 Python 实际用到的置换导出来，R 侧喂同一组，p 值才能逐位比。只对
    parity_samples 指定的样本导出，否则 fixture 体积会到几十 MB。
    """
    from backend.engines import toolkit_core as tc
    from backend.engines import toolkit_tools as tt
    from backend.engines import sample_identity as si

    cent_all = json.load(open(BACKEND_DATA / "toolkit_refs.json",
                              encoding="utf-8"))["attribution"]["centroids"]
    cent = [c for c in cent_all if c["space"] == "gene"]
    ref_cent = [c for c in cent if c.get("role") == "celltype_reference"]
    ery_cent = [c for c in cent if c.get("role", "erythroid") == "erythroid"]

    names, X = tc.to_space(x, "gene")
    names = [str(n) for n in names]
    X = np.where(np.isfinite(X), X, np.nanmedian(X, axis=1, keepdims=True))
    Q = np.apply_along_axis(
        lambda v: tc.rank_profile(v) / max(1.0, (np.isfinite(v).sum() - 1)), 0, X)
    qpos = {n: i for i, n in enumerate(names)}
    n_feat = len(names)

    def aligned(pool):
        out = []
        for c in pool:
            shared = [(i, qpos[n]) for i, n in enumerate(c["features"]) if n in qpos]
            if len(shared) >= tc.MIN_SHARED:
                out.append({
                    "dataset_id": c["dataset_id"], "group": c["group"],
                    "disease": c["disease"], "layer": c["layer"], "n": c["n"],
                    "role": c.get("role", "erythroid"), "shared": len(shared),
                    "_qi": np.array([q for _, q in shared], int),
                    "_cv": np.array([c["profile"][i] for i, _ in shared], float),
                })
        return out

    al = aligned(ery_cent)
    al_ref = aligned(ref_cent)
    sc = si.scope(list(x.columns), pd.DataFrame(), {})
    idmap = {str(s["sample"]): s for s in sc["samples"]}

    def null_with_perms(q, members, perms):
        refs_rb = []
        for m in members:
            cv = np.asarray(m["_cv"], float)
            rb = tc.rank_profile(cv)
            rb = rb - rb.mean()
            nb = np.sqrt(float((rb * rb).sum()))
            refs_rb.append((rb, nb))
        fat = []
        for m, (rb, nb) in zip(members, refs_rb):
            A = q[perms[:, m["_qi"]]]
            ra = tt._rank_rows(A)
            ra = ra - ra.mean(axis=1, keepdims=True)
            na = np.sqrt((ra * ra).sum(axis=1))
            with np.errstate(divide="ignore", invalid="ignore"):
                rho = (ra @ rb) / (na * nb)
            rho[~np.isfinite(rho)] = -np.inf
            fat.append(rho)
        best = np.vstack(fat).max(axis=0)
        return best[np.isfinite(best)]

    rng = np.random.default_rng(20261110)
    rows, perms_out = [], {}
    for j, sample in enumerate(x.columns):
        q = Q[:, j]
        hits = []
        for a in al:
            rho = tc.spearman(q[a["_qi"]], a["_cv"])
            if np.isfinite(rho):
                hits.append({**{k: v for k, v in a.items() if not k.startswith("_")},
                             "rho": round(float(rho), 3), "_raw": float(rho), "_a": a})
        if not hits:
            rows.append({"sample": str(sample), "ok": 0, "top_label": "",
                         "best_rho": np.nan, "mean_rho": np.nan,
                         "background_median_rho": np.nan, "n_compared": 0,
                         "runner_label": "", "runner_best_rho": np.nan,
                         "consistency": np.nan, "p_value": np.nan,
                         "null_q95": np.nan, "margin": np.nan,
                         "identity": "", "n_conditions": 0})
            continue
        ref_best = None
        for a in al_ref:
            rho = tc.spearman(q[a["_qi"]], a["_cv"])
            if np.isfinite(rho) and (ref_best is None or rho > ref_best["rho"]):
                ref_best = {"rho": round(float(rho), 3), "role": a["role"]}

        by_cond = {}
        for h in hits:
            by_cond.setdefault(h["disease"], []).append(h)
        conds = []
        for label, hs in by_cond.items():
            rhos = [h["rho"] for h in hs]
            best = max(hs, key=lambda h: h["rho"])
            conds.append({"label": label, "n_cohorts": len(hs),
                          "mean_rho": round(float(np.mean(rhos)), 3),
                          "best_rho": round(float(best["rho"]), 3),
                          "best_dataset": best["dataset_id"]})
        conds.sort(key=lambda c: (-c["best_rho"], -c["mean_rho"]))
        topc = conds[0]
        members = [h["_a"] for h in by_cond[topc["label"]]]

        # 只有指定样本导出置换；其余样本用 Python 自己的 RNG，p 只作参考
        if parity_samples is not None and j in parity_samples:
            perms = np.empty((parity_n, n_feat), dtype=np.int64)
            for t in range(parity_n):
                perms[t] = rng.permutation(n_feat)
            null = null_with_perms(q, members, perms)
            perms_out[j] = perms
        else:
            perms = None
            for _ in range(n_perm):
                rng.permutation(n_feat)      # 保持与引擎相同的调用次数
            null = _null_reference(q, members, n_perm, rng)

        p = (float((1 + np.sum(null >= topc["best_rho"])) / (1 + len(null)))
             if len(null) else np.nan)
        bg = float(np.median([h["rho"] for h in hits]))
        ident = idmap.get(str(sample)) or {}
        rows.append({
            "sample": str(sample), "ok": 1,
            "top_label": topc["label"], "best_rho": topc["best_rho"],
            "mean_rho": topc["mean_rho"], "background_median_rho": round(bg, 3),
            # 未舍入：舍入值在 .0005 边界上 R 与 numpy 会差一个末位，用这个列做
            # 1e-9 级别的严格对拍，才分得清「浮点末位」和「算法真的不一样」
            "best_rho_raw": float(max(h["_raw"] for h in by_cond[topc["label"]])),
            "mean_rho_raw": float(np.mean([h["_raw"] for h in by_cond[topc["label"]]])),
            "bg_raw": bg,
            "n_compared": len(hits),
            "runner_label": conds[1]["label"] if len(conds) > 1 else "",
            "runner_best_rho": conds[1]["best_rho"] if len(conds) > 1 else np.nan,
            "consistency": (round(topc["mean_rho"] / topc["best_rho"], 3)
                            if topc["best_rho"] else np.nan),
            "p_value": round(p, 4) if np.isfinite(p) else np.nan,
            "null_q95": round(float(np.quantile(null, 0.95)), 3) if len(null) else np.nan,
            "margin": round(topc["best_rho"] - bg, 3),
            "identity": ident.get("compartment", ""),
            "n_conditions": len(conds),
        })
    return pd.DataFrame(rows), perms_out


def _null_reference(q, members, n_perm, rng):
    """引擎逐次实现的等价版，用于没导出置换的样本。"""
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(len(q))
        qq = q[perm]
        vals = [tc.spearman(qq[m["_qi"]], m["_cv"]) for m in members]
        vals = [v for v in vals if np.isfinite(v)]
        if vals:
            null.append(float(np.max(vals)))
    return np.asarray(null)


def age_expected(n_stage_rows: int = 5):
    """红系衰老：把参考面板自己投回去，五个阶段应各自落在 0/1/2/3/4。"""
    from backend.engines import age_axis as aa

    mat = pd.read_csv(aa.REF_PATH, sep="\t", index_col=0).apply(
        pd.to_numeric, errors="coerce")
    mat = mat[mat.notna().all(axis=1)]
    mat.insert(0, "lipid", mat.index)
    mat.to_csv(OUT / "age_lipid_input.tsv.gz", sep="\t", index=False,
               compression="gzip", float_format="%.6g")

    res = aa.project_matrix(mat.set_index("lipid"))
    v = aa.validation()
    rows = []
    for s in res["samples"]:
        rows.append({
            "sample": s["sample"], "ok": int(bool(s.get("ok"))),
            "position": s.get("position", np.nan), "rho": s.get("rho", np.nan),
            "p_perm": s.get("p_perm", np.nan),
            "band_lo": (s.get("band") or [np.nan, np.nan])[0],
            "band_hi": (s.get("band") or [np.nan, np.nan])[1],
            "n_matched": s.get("n_matched", np.nan),
            "n_within": s.get("n_within", np.nan),
            "p_binom": s.get("p_value", np.nan),
            "zone_key": (s.get("zone") or {}).get("key", ""),
            "tier_A": (s.get("tiers") or {}).get("A", np.nan),
            "tier_B": (s.get("tiers") or {}).get("B", np.nan),
            "reason": s.get("reason", ""),
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "age_expected.csv", index=False)

    # 形状检验的置换矩阵：_perm_p 每次调用都重新 default_rng(seed)，所以五个
    # 样本（keep 数都是 45）共用同一组排列。导出来才能让 p_perm 逐位对拍。
    n_keep = int(df["n_matched"].iloc[0])
    rng = np.random.default_rng(20261112)
    perms = np.empty((200, n_keep), dtype=np.int64)
    for t in range(200):
        perms[t] = rng.permutation(n_keep)
    with gzip.open(OUT / "age_perms.csv.gz", "wt") as fh:
        np.savetxt(fh, perms + 1, fmt="%d", delimiter=",")
    print(f"  age perms: {perms.shape}")
    cov = res["coverage"]
    meta = {
        "loso_mae": v["loso_mae"], "loso_max": v["loso_max"],
        "loso_mae_pct": v["loso_mae_pct"], "loso_errors": v["loso_errors"],
        "n_matched": cov["n_matched"], "n_anchor": cov["n_anchor"],
        "tier_A": cov["tier_A"], "tier_B": cov["tier_B"],
        "n_ok": res["summary"]["n_ok"], "n_refused": res["summary"]["n_refused"],
        "median_position": res["summary"]["median_position"],
        "median_zone": res["summary"]["median_zone"]["key"],
        "n_young": res["summary"]["n_young"], "n_old": res["summary"]["n_old"],
        "n_between": res["summary"]["n_between"],
    }
    json.dump(meta, open(OUT / "age_model_expected.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    # 一条弃权路径要有对拍：只留 6 个脂质 → insufficient_overlap
    short = mat.set_index("lipid").iloc[:6]
    r2 = aa.project_matrix(short)
    json.dump({"ok": bool(r2["ok"]), "reason": r2.get("reason", ""),
               "n_samples": len(r2["samples"]),
               "first_reason": (r2["samples"][0].get("reason") or "")},
              open(OUT / "age_refusal_expected.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return df


def age_prior_expected(x: pd.DataFrame, tag: str = ""):
    """脂质调控先验外推（C 档）：基因层输入 → 批内相对秩位置。

    tag 用来区分两条路径的期望文件：marker 子集走弃权，84 个先验基因走成功。
    """
    from backend.engines import age_prior as ap

    res = ap.project_prior(x)
    out = {
        "ok": bool(res.get("ok")), "unverifiable": bool(res.get("unverifiable")),
        "reason_code": res.get("reason_code", ""),
        "n_hit": res.get("n_hit"), "n_directional": res.get("n_directional"),
        "n_context_only": res.get("n_context_only"),
        "coverage": res.get("coverage"),
        "direction_agreement": res.get("direction_agreement"),
        "n_tier_a": res.get("n_tier_a"), "n_tier_b": res.get("n_tier_b"),
        "feature_loo_r2": res.get("feature_loo_r2"),
        "confidence": res.get("confidence"), "grade": res.get("grade"),
        "confidence_parts": res.get("confidence_parts"),
        "diagnostics": {
            "r2": (res.get("diagnostics") or {}).get("r2"),
            "loo_r2": (res.get("diagnostics") or {}).get("loo_r2"),
            "stage_recovery_rho": (res.get("diagnostics") or {}).get("stage_recovery_rho"),
            "beta": (res.get("diagnostics") or {}).get("beta"),
            "directions": (res.get("diagnostics") or {}).get("directions"),
            "permutation_p": (res.get("diagnostics") or {}).get("permutation_p"),
        },
    }
    name = ("age_prior_expected.json" if not tag
            else f"age_prior_{tag}_expected.json")
    json.dump(out, open(OUT / name, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if res.get("ok"):
        pd.DataFrame(res["samples"]).to_csv(OUT / "age_prior_samples_expected.csv",
                                            index=False)
    return out


def lipid_parse_expected():
    """脂质名解析的期望值：面板全部 151 行 + 文献里常见的变体写法。

    解析器错一个字符，某些分子种就永远对不上参考轨迹，而症状只是"命中数偏低"，
    看不出是解析错了。所以逐名对拍。
    """
    from backend.engines import age_axis as aa

    names = list(pd.read_csv(aa.REF_PATH, sep="\t", usecols=[0])["lipid"].astype(str))
    names += [
        "PC(34:1)", "PC(34:1)e", "PC(O-34:1)", "PC(P-36:2)", "PE(38:4)",
        "SM(t42:3)", "SM(d42:1)", "SM(42:1)", "TG(18:2_17:1_18:2)",
        "TG(16:0/18:1/18:2)", "DG(16:0-18:1)", "Cer(d18:1/16:0)",
        "SPHINGOMYELIN 34:1", "PHOSPHATIDYLCHOLINE(36:2)",
        "LPC(18:0)", "LPE(20:4)", "Cholesterol", "CE(18:2)",
        "PI(38:4)", "PS(36:1)", "PC:0:0", "PE:40:7 ()", "  PC:34:1 (16:0-18:1w9)  ",
        "PC:34", "PC", "34:1", "PC:x:y", "PC(16:0_18:1)", "PE(d18:1-16:0)",
    ]
    rows = []
    for nm in names:
        p = aa.parse_lipid(nm)
        rows.append({
            "name": nm,
            "cls": "" if not p else p["cls"],
            "c": np.nan if not p else p["c"],
            "db": np.nan if not p else p["db"],
            "ether": "" if not p else ("TRUE" if p["ether"] else "FALSE"),
            "base": "" if not p else p["base"],
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "lipid_parse_expected.csv", index=False)
    print(f"  lipid names: {len(df)} (parsed {int((df['cls'] != '').sum())})")
    return df


def age_feature_perms():
    """特征模型的置换矩阵（1000 × 45），让 permutation_p 可逐位对拍。

    feature_model() 里这个序列由 numpy 的 PCG64 产生，R 没法复现；把实际用到的
    1000 个排列导出来喂给 R，才能证明"检验统计量一致"而不是"两边都算出一个小数"。
    """
    from backend.engines import age_axis as aa
    from backend.engines import age_prior as ap

    n = len(aa.model()["keys"])
    rng = np.random.default_rng(ap.PERM_SEED)
    perms = np.empty((ap.PERM_N, n), dtype=np.int64)
    for t in range(ap.PERM_N):
        perms[t] = rng.permutation(n)
    with gzip.open(OUT / "age_feature_perms.csv.gz", "wt") as fh:
        np.savetxt(fh, perms + 1, fmt="%d", delimiter=",")
    print(f"  feature perms: {perms.shape}")
    return perms


def main():
    refs = json.load(open(BACKEND_DATA / "layer_references.json", encoding="utf-8"))
    markers = refs["markers"]["all"]
    print(f"markers: {len(markers)}")

    # 转录组：纯 RBC / PBMC / 储存 8 天三档，端点判别力最强的一个数据集
    sub = write_subset("GSE63703", markers, "gse63703_subset.tsv.gz")
    exp = purity_expected(sub)
    exp.to_csv(OUT / "purity_expected.csv", index=False)
    print(exp.to_string(index=False))

    pexp = projection_expected(sub)
    pexp.to_csv(OUT / "projection_expected.csv", index=False)
    print(pexp.to_string(index=False))

    # 通路指纹（within / cohort 两种归一）
    for md in ("within", "cohort"):
        pex = pathway_expected(sub, md)
        pex.to_csv(OUT / f"pathway_{md}_expected.csv", index=False)
        print(f"  pathway {md}: {pex.shape[0]} pathways")
        print(pex.head(4).to_string(index=False))

    # 条件谱匹配（RBC_1 与 PBMC_1 导出置换，p 可逐位比）
    aexp, aperms = attribution_expected(sub, parity_samples=[0, 3], parity_n=15)
    aexp.to_csv(OUT / "attribution_expected.csv", index=False)
    print(aexp.to_string(index=False))
    for j, perms in aperms.items():
        with gzip.open(OUT / f"attribution_perms_{j}.csv.gz", "wt") as fh:
            np.savetxt(fh, perms + 1, fmt="%d", delimiter=",")
        print(f"  perms for sample {j}: {perms.shape}")

    # 特征模型置换（让 permutation_p 可比）
    age_feature_perms()

    # 脂质名解析（逐名对拍）
    lipid_parse_expected()

    # 红系衰老（参考面板自投影：五阶段应落在 0..4）
    aexp2 = age_expected()
    print(aexp2.to_string(index=False))

    # 脂质调控先验外推（C 档）。marker 子集只命中 7 条带方向先验（阈值 8）→
    # 走弃权分支；为了也覆盖**成功路径**，另用全矩阵里命中的 84 个先验基因。
    apx = age_prior_expected(sub, tag="markers")
    print("  age prior (markers, expect abstain):",
          {k: apx[k] for k in ("ok", "reason_code", "n_directional")})

    from backend.engines import age_prior as _ap
    full = load_matrix("GSE63703")
    pri = {str(r["gene"]).upper() for r in _ap.priors()}
    sel = [i for i in full.index if str(i).upper() in pri]
    sub_prior = full.loc[sel]
    sub_prior.insert(0, "gene", sub_prior.index)
    sub_prior.to_csv(OUT / "age_prior_input.tsv.gz", sep="\t", index=False,
                     compression="gzip", float_format="%.6g")
    apx2 = age_prior_expected(sub_prior.set_index("gene"))
    print("  age prior (84 priors, expect ok):",
          {k: apx2[k] for k in ("ok", "grade", "coverage", "direction_agreement",
                                "n_directional", "confidence")})

    # miRNA 时钟（同层模型）
    mx = load_matrix("GSE114990")
    mx.insert(0, "feature", mx.index)
    mx.to_csv(OUT / "gse114990_mirna.tsv.gz", sep="\t", index=False,
              compression="gzip", float_format="%.6g")
    cexp = clock_expected("GSE114990", "mirna")
    if len(cexp):
        cexp.to_csv(OUT / "clock_mirna_expected.csv", index=False)
        print(f"  mirna clock: {cexp.shape[0]} samples")
        print(cexp.head(8).to_string(index=False))

    # 代谢组时钟（同层 + 模块桥，两条路径都测）
    sx = load_matrix("ST004148_metabolomics")
    sx.insert(0, "feature", sx.index)
    sx.to_csv(OUT / "st004148_metabolome.tsv.gz", sep="\t", index=False,
              compression="gzip", float_format="%.6g")
    sexp = clock_expected("ST004148_metabolomics", "metabolite")
    if len(sexp):
        sexp.to_csv(OUT / "clock_metabolite_expected.csv", index=False)
        print(f"  metabolite clock: {sexp.shape[0]} samples")
        print(sexp.head(8).to_string(index=False))
    mexp = module_clock_expected("ST004148_metabolomics", "metabolite")
    if len(mexp):
        mexp.to_csv(OUT / "clock_module_expected.csv", index=False)
        print(f"  module clock: {mexp.shape[0]} samples")
        print(mexp.head(8).to_string(index=False))

    meta = {
        "seed": SEED,
        "n_background": N_BACKGROUND,
        "engine_version": open(REPO / "VERSION.txt", encoding="utf-8").readline().strip(),
        "datasets": {"gse63703_subset.tsv.gz": "GSE63703 (mature RBC vs PBMC vs 8-day stored)"},
        "generated_by": "data-raw/make_fixtures.py",
    }
    json.dump(meta, open(OUT / "fixtures.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    total = sum(f.stat().st_size for f in OUT.iterdir())
    print(f"\nfixtures: {len(list(OUT.iterdir()))} files / {total/1024:.0f} KB")


if __name__ == "__main__":
    main()
