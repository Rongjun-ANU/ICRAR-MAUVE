import matplotlib.pyplot as plt
import numpy as np

def plot_rho_on_ax(ax, result, title):
    observed = result["observed"]
    centers = observed["centers"]
    ax.fill_between(
        centers,
        result["null"]["p16_rho"],
        result["null"]["p84_rho"],
        color="tab:orange",
        alpha=0.25,
        label="Null 16-84%",
    )
    ax.plot(
        centers,
        result["null"]["median_rho"],
        color="tab:orange",
        linewidth=2.0,
        label="Null median",
    )
    ax.plot(
        centers,
        observed["rho"],
        color="black",
        linewidth=2.2,
        marker="o",
        label="Observed (error-qualified subset)",
    )
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_xlabel(r"$\log\,\Sigma_*$ bin center", fontsize=11)
    ax.set_ylabel(r"Spearman $\rho$", fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=10)
