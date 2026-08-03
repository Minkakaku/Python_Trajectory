"""Portable Scanpy PAGA trajectory workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def save_current_figure(stem: Path, width: float = 8, height: float = 7) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_paga(
    input_h5ad: str,
    groupby: str,
    output_dir: str,
    object_out: str,
    threshold: float,
    n_pcs: int,
    n_neighbors: int,
) -> None:
    output = Path(output_dir).expanduser().resolve()
    object_path = Path(object_out).expanduser().resolve()
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

    adata.write_h5ad(object_path)
    (output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# PAGA 轨迹分析结果",
                "",
                f"- 输入对象：{Path(input_h5ad).expanduser().resolve()}",
                f"- 分组列：{groupby}",
                f"- edge threshold：{threshold}",
                f"- PAGA 对象：{object_path}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Scanpy PAGA trajectory analysis.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--groupby", required=True)
    parser.add_argument("--output-dir", default="results/04_trajectory/paga")
    parser.add_argument("--object-out", default="data/processed/paga.h5ad")
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_paga(
        input_h5ad=args.input_h5ad,
        groupby=args.groupby,
        output_dir=args.output_dir,
        object_out=args.object_out,
        threshold=args.threshold,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
    )