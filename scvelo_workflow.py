"""scVelo 独立 RNA velocity 入口，要求 h5ad 保留 spliced/unspliced layers。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import scvelo as scv


# 功能：把当前 Matplotlib 图同时保存为 PDF 和 300 dpi PNG。
def save_current_figure(stem: Path, width: float = 9, height: float = 7) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


# 功能：运行 scVelo moments、velocity graph、velocity pseudotime 和可选 latent time。
def run_velocity(
    input_h5ad: str,
    output_dir: str,
    object_out: str | None,
    groupby: str,
    basis: str,
    mode: str,
    min_shared_counts: int,
    n_top_genes: int,
    n_pcs: int,
    n_neighbors: int,
) -> None:
    output = Path(output_dir).expanduser().resolve()
    object_path = (
        Path(object_out).expanduser().resolve()
        if object_out
        else output / "scvelo.h5ad"
    )
    output.mkdir(parents=True, exist_ok=True)
    object_path.parent.mkdir(parents=True, exist_ok=True)

    adata = scv.read(input_h5ad, cache=False)
    missing_layers = {"spliced", "unspliced"}.difference(adata.layers)
    if missing_layers:
        raise KeyError(
            "h5ad 缺少 velocity layers："
            + ", ".join(sorted(missing_layers))
            + "。普通表达矩阵不能补造 RNA velocity。"
        )
    if groupby not in adata.obs:
        raise KeyError(f"groupby column not found: {groupby}")

    # scVelo 按 spliced/unspliced 共同计数过滤并计算一、二阶 moments。
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=min_shared_counts,
        n_top_genes=n_top_genes,
    )
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    if mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    if mode == "dynamical":
        # latent time 使用动力学参数和 velocity graph，必须在建图后计算。
        scv.tl.latent_time(adata)

    scv.pl.velocity_embedding_stream(
        adata, basis=basis, color=groupby, show=False
    )
    save_current_figure(output / "velocity_stream")
    scv.pl.velocity_embedding(
        adata,
        basis=basis,
        color=groupby,
        arrow_length=3,
        show=False,
    )
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
                f"- 分组列/embedding：{groupby}/{basis}",
                f"- min_shared_counts：{min_shared_counts}；top genes：{n_top_genes}",
                f"- scVelo 对象：{object_path}",
                "- 输入必须保留由 velocyto/kallisto|bustools 等流程得到的 spliced/unspliced layers。",
                "- 流线方向是转录动力学模型推断；需结合已知生物学方向、采样设计和其他轨迹方法复核。",
                "- velocity pseudotime/latent time 是相对动力学顺序，不等于真实时间。",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


# 功能：定义 velocity 输入层要求、模型和邻居参数。
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run scVelo from h5ad with spliced/unspliced layers.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--groupby", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--object-out", default=None)
    parser.add_argument("--basis", default="umap")
    parser.add_argument(
        "--mode",
        choices=["deterministic", "stochastic", "dynamical"],
        default="stochastic",
    )
    parser.add_argument("--min-shared-counts", type=int, default=20)
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
        min_shared_counts=args.min_shared_counts,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
    )