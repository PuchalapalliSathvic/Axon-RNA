"""Scientific figures sourced exclusively from real data and out-of-fold predictions."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT/".cache/matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay, ConfusionMatrixDisplay

plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                     "figure.dpi":120, "savefig.dpi":220})
DEST = ROOT/"results/figures"
COLORS = ["#446b9e", "#d89135", "#238878", "#a45a89"]
NAMES = {"sequence_only":"Sequence", "sequence_structure":"Sequence + structure",
         "sequence_knowledge":"Sequence + knowledge", "sequence_structure_knowledge":"Sequence + structure + knowledge"}


def save(fig, name):
    fig.savefig(DEST/f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(ROOT/"results/strong_candidates.csv").head(20)
    studies = pd.read_csv(ROOT/"data/real/processed/all_real_studies.csv")
    comparison = pd.read_csv(ROOT/"results/model_comparison.csv")
    pred = pd.read_csv(ROOT/"results/model_predictions.csv")
    data = pd.read_csv(ROOT/"data/real/modeling/real_model_dataset.csv")
    data = data[data.included_in_models]
    fig, ax = plt.subplots(figsize=(8,6), layout="constrained")
    view = candidates.iloc[::-1]
    ax.barh(view.transcript_id, view.median_log2fc, color=np.where(view.median_log2fc<0,COLORS[0],COLORS[1]))
    ax.axvline(0,color="0.3",lw=.7)
    ax.set(title="Real AD cohorts: top 20 strong consensus candidates", xlabel="Median study effect (difference in mean log₂(CPM + 1))")
    save(fig,"top_consensus_genes")

    heat = studies.pivot(index="transcript_id",columns="study_id",values="log2fc").reindex(candidates.transcript_id)
    fig,ax = plt.subplots(figsize=(6,7),layout="constrained")
    limit = np.abs(heat.to_numpy()).max()
    im=ax.imshow(heat,cmap="RdBu_r",vmin=-limit,vmax=limit,aspect="auto")
    ax.set(xticks=range(len(heat.columns)),xticklabels=heat.columns,yticks=range(len(heat)),yticklabels=heat.index,
           title="Real AD cohorts: direction and effect by study")
    fig.colorbar(im,ax=ax,label="Difference in mean log₂(CPM + 1)")
    save(fig,"consensus_direction_heatmap")

    fig,ax=plt.subplots(figsize=(9,4.5),layout="constrained")
    x=np.arange(len(comparison))
    ax.bar(x-.17,comparison.roc_auc,.34,label="ROC-AUC",color=COLORS[0])
    ax.bar(x+.17,comparison.average_precision,.34,label="Average precision",color=COLORS[2])
    for i,r in comparison.iterrows():
        ax.text(i-.17,r.roc_auc+.012,f"{r.roc_auc:.3f}",ha="center",fontsize=9)
        ax.text(i+.17,r.average_precision+.012,f"{r.average_precision:.3f}",ha="center",fontsize=9)
    ax.set(xticks=x,xticklabels=[NAMES[n].replace(" + ","\n+ ") for n in comparison.model],ylim=(0,1.07),
           ylabel="Pooled held-out score",title=f"Real gene-level labels • {len(pred)} genes • identical five-fold splits")
    ax.legend(loc="lower right")
    save(fig,"model_comparison")

    for kind,display in [("roc",RocCurveDisplay),("precision_recall",PrecisionRecallDisplay)]:
        fig,ax=plt.subplots(figsize=(7,5.5),layout="constrained")
        for (name,label),color in zip(NAMES.items(),COLORS):
            display.from_predictions(pred.label,pred[name],name=label,color=color,ax=ax)
        if kind=="roc":
            ax.plot([0,1],[0,1],"--",color="0.6",lw=1)
        else:
            ax.axhline(pred.label.mean(),ls="--",color="0.6",label="Reference prevalence")
        ax.set(title="Real gene labels: out-of-fold predictions")
        ax.legend(fontsize=8,loc="lower right" if kind=="roc" else "lower left")
        save(fig,f"{kind}_curve")

    fig,axes=plt.subplots(2,2,figsize=(9,8),layout="constrained")
    for ax,(name,label) in zip(axes.flat,NAMES.items()):
        ConfusionMatrixDisplay.from_predictions(pred.label,pred[name]>=.5,display_labels=["Reference","Strong"],
                                                colorbar=False,cmap="Blues",ax=ax)
        ax.set_title(label,fontsize=10)
    fig.suptitle("Real gene labels: held-out classifications at threshold 0.5")
    save(fig,"confusion_matrix")

    fig,axes=plt.subplots(1,3,figsize=(10,4),layout="constrained")
    for ax,col,label in zip(axes,["normalized_mfe","paired_fraction","window_gc_content"],
                           ["MFE / folded nt (kcal/mol/nt)","Paired-base fraction","Window GC fraction"]):
        groups=[data.loc[data.label==i,col].dropna().to_numpy() for i in [0,1]]
        ax.boxplot(groups,tick_labels=["Reference","Strong"],showfliers=False)
        rng=np.random.default_rng(42)
        for i,g in enumerate(groups,1):
            ax.scatter(rng.normal(i,.045,len(g)),g,s=8,alpha=.35,color=COLORS[i-1])
        ax.set_ylabel(label)
    fig.suptitle("Real representative transcripts: 5′ windows up to 1,000 nt")
    save(fig,"structure_feature_comparison")

    imp=pd.read_csv(ROOT/"results/knowledge_feature_importance.csv")
    agg=imp.groupby("feature").auc_decrease.agg(["mean","std"]).sort_values("mean")
    fig,ax=plt.subplots(figsize=(8,4.5),layout="constrained")
    ax.barh(agg.index,agg["mean"],xerr=agg["std"],color=COLORS[2],capsize=3)
    ax.axvline(0,color="0.3",lw=.8)
    ax.set(xlabel="Held-out ROC-AUC decrease after permutation (mean ± fold SD)",
           title="Generic knowledge importance • sequence + knowledge model")
    save(fig,"knowledge_feature_importance")

    consensus=pd.read_csv(ROOT/"results/consensus_genes.csv")
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    axes[0].axis("off")
    summary=(f"AXON-RNA | real-data analysis\n\n"
             f"3 human AD brain cohorts\n{len(consensus):,} genes in the union\n"
             f"{int((consensus.n_studies==3).sum()):,} measured in all three\n"
             f"{int(consensus.ad_dysregulated.sum()):,} consensus-positive (screening)\n"
             f"65 strong candidates\n{len(pred)} modeled genes; {int(pred.label.sum())} strong\n\n"
             "Gene-level labels; representative transcript inputs\n"
             "Five-fold held-out comparison, seed 42\n"
             "Exploratory; no external-cohort validation")
    axes[0].text(.02,.96,summary,va="top",fontsize=12,linespacing=1.6)
    axes[1].barh([NAMES[n] for n in comparison.model],comparison.roc_auc,color=COLORS)
    axes[1].set(xlim=(0,1),xlabel="Pooled held-out ROC-AUC")
    axes[1].axvline(.5,ls="--",color="0.5")
    for i,v in enumerate(comparison.roc_auc):
        axes[1].text(v+.01,i,f"{v:.3f}",va="center")
    save(fig,"project_summary")


if __name__ == "__main__":
    main()
