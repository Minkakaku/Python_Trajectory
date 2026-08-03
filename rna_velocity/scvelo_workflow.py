"""Portable scVelo RNA-velocity workflow for h5ad files with spliced/unspliced layers."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scvelo as scv


def save_current_figure(stem: Path, width: float = 9, height: float = 7) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_velocity(
    input_h5ad: str,
    output_dir: str,
    object_out: str,
    groupby: str,
    basis: str,
    mode: str,
    n_top_genes: int,
    n_pcs: int,
    n_neighbors: int,
) -> None:
    output = Path(output_dir).expanduser().resolve()
    object_path = Path(object_out).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    object_path.parent.mkdir(parents=True, exist_ok=True)

    adata = scv.read(input_h5ad, cache=False)
    missing_layers = {"spliced", "unspliced"}.difference(adata.layers)
    if missing_layers:
        raise KeyError(f"Missing velocity layers: {', '.join(sorted(missing_layers))}")
    if groupby not in adata.obs:
        raise KeyError(f"groupby column not found: {groupby}")

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=n_top_genes)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    if mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.latent_time(adata)
    else:
        scv.tl.velocity(adata, mode=mode)

    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)

    scv.pl.velocity_embedding_stream(adata, basis=basis, color=groupby, show=False)
    save_current_figure(output / "velocity_stream")

    scv.pl.velocity_embedding(adata, basis=basis, color=groupby, arrow_length=3, show=False)
    save_current_figure(output / "velocity_arrows")

    columns = [groupby, "velocity_pseudotime"]
    if "latent_time" in adata.obs:
        columns.append("latent_time")
    adata.obs[columns].reset_index(names="cell").to_csv(
        output / "velocity_pseudotime.csv",
        index=False,
    )

    adata.write_h5ad(object_path)
    (output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# scVelo RNA velocity 结果",
                "",
                f"- 输入对象：{Path(input_h5ad).expanduser().resolve()}",
                f"- 模式：{mode}",
                f"- 分组列：{groupby}",
                f"- embedding：{basis}",
                f"- scVelo 对象：{object_path}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run scVelo RNA-velocity analysis.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--groupby", required=True)
    parser.add_argument("--output-dir", default="results/04_trajectory/rna_velocity")
    parser.add_argument("--object-out", default="data/processed/scvelo.h5ad")
    parser.add_argument("--basis", default="umap")
    parser.add_argument(
        "--mode",
        choices=["deterministic", "stochastic", "dynamical"],
        default="stochastic",
    )
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_velocity(
        input_h5ad=args.input_h5ad,
        output_dir=args.output_dir,
        object_out=args.object_out,
        groupby=args.groupby,
        basis=args.basis,
        mode=args.mode,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
    )