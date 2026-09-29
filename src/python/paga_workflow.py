"""Scanpy PAGA 独立轨迹入口：h5ad → 群间拓扑；可选 root → DPT pseudotime。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# 功能：把当前 Matplotlib 图同时保存为 PDF 和 600 dpi PNG。
def save_current_figure(stem: Path, width: float = 8, height: float = 7) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(figure)


# 功能：计算 PAGA 连通拓扑，并在指定生物学 root 时计算 DPT pseudotime。
def run_paga(
    input_h5ad: str,
    groupby: str,
    root_group: str | None,
    output_dir: str,
    object_out: str | None,
    threshold: float,
    n_pcs: int,
    n_neighbors: int,
) -> None:
    output = Path(output_dir).expanduser().resolve()
    object_path = (
        Path(object_out).expanduser().resolve()
        if object_out
        else output / "paga.h5ad"
    )
    output.mkdir(parents=True, exist_ok=True)
    object_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_h5ad)
    if groupby not in adata.obs:
        raise KeyError(f"groupby column not found: {groupby}")
    adata.obs[groupby] = adata.obs[groupby].astype("category")

    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.paga(adata, groups=groupby)

    sc.pl.paga(adata, color=groupby, threshold=threshold, show=False)
    save_current_figure(output / "paga_graph")
    sc.pl.paga_compare(
        adata,
        color=groupby,
        threshold=threshold,
        node_size_scale=2,
        edge_width_scale=0.4,
        show=False,
    )
    save_current_figure(output / "paga_compare", width=12, height=6)

    categories = adata.obs[groupby].cat.categories
    connectivities = adata.uns["paga"]["connectivities"].toarray()
    pd.DataFrame(connectivities, index=categories, columns=categories).to_csv(
        output / "paga_connectivities.csv",
        index=True,
        index_label=groupby,
    )

    # PAGA 只估计群间拓扑；有明确起始群时才计算具有方向含义的 DPT。
    if root_group is not None:
        root_cells = np.flatnonzero(
            adata.obs[groupby].astype(str).to_numpy() == root_group
        )
        if len(root_cells) == 0:
            raise ValueError(f"root_group not found in {groupby}: {root_group}")
        adata.uns["iroot"] = int(root_cells[0])
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
        adata.obs[[groupby, "dpt_pseudotime"]].reset_index(names="cell").to_csv(
            output / "paga_dpt_pseudotime.csv",
            index=False,
        )
        if "X_umap" in adata.obsm:
            sc.pl.umap(adata, color="dpt_pseudotime", show=False)
            save_current_figure(output / "paga_dpt_pseudotime")

    adata.write_h5ad(object_path)
    (output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# PAGA 轨迹分析结果",
                "",
                f"- 输入对象：{Path(input_h5ad).expanduser().resolve()}",
                f"- 分组列：{groupby}",
                f"- root group：{root_group if root_group is not None else '未指定'}",
                f"- edge threshold：{threshold}",
                f"- PAGA 对象：{object_path}",
                "- PAGA 图表示群间连通拓扑；未指定 root 时不能解释为有方向的发育顺序。",
                "- 指定 root 后的 DPT 是图上的相对顺序，不等于真实时间。",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


# 功能：定义 h5ad、分群列、可选 root 和图参数。
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Scanpy PAGA and optional DPT.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--groupby", required=True)
    parser.add_argument("--root-group", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--object-out", default=None)
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_paga(
        input_h5ad=args.input_h5ad,
        groupby=args.groupby,
        root_group=args.root_group,
        output_dir=args.output_dir,
        object_out=args.object_out,
        threshold=args.threshold,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
    )
