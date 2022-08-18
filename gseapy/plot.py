# -*- coding: utf-8 -*-
import operator
import sys
from typing import Dict, Iterable, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def zscore(data2d, axis=0):
    """Standardize the mean and variance of the data axis Parameters.

    :param data2d: DataFrame to normalize.
    :param axis: int, Which axis to normalize across. If 0, normalize across rows,
                  if 1, normalize across columns. If None, don't change data

    :Returns: Normalized DataFrame. Normalized data with a mean of 0 and variance of 1
              across the specified axis.

    """
    if axis is None:
        # normalized to mean and std using entire matrix
        # z_scored = (data2d - data2d.values.mean()) / data2d.values.std(ddof=1)
        return data2d
    assert axis in [0, 1]
    z_scored = data2d.apply(
        lambda x: (x - x.mean()) / x.std(ddof=1), axis=operator.xor(1, axis)
    )
    return z_scored


def _skip_ticks(labels, tickevery):
    """Return ticks and labels at evenly spaced intervals."""
    n = len(labels)
    if tickevery == 0:
        ticks, labels = [], []
    elif tickevery == 1:
        ticks, labels = np.arange(n) + 0.5, labels
    else:
        start, end, step = 0, n, tickevery
        ticks = np.arange(start, end, step) + 0.5
        labels = labels[start:end:step]
    return ticks, labels


def _auto_ticks(ax, labels, axis):
    transform = ax.figure.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(transform)
    size = [bbox.width, bbox.height][axis]
    axis = [ax.xaxis, ax.yaxis][axis]
    (tick,) = ax.xaxis.set_ticks([0])
    fontsize = tick.label1.get_size()
    max_ticks = int(size // (fontsize / 72))
    if max_ticks < 1:
        tickevery = 1
    else:
        tickevery = len(labels) // max_ticks + 1
    return tickevery


def heatmap(
    df,
    z_score=None,
    title="",
    figsize=(5, 5),
    cmap="RdBu_r",
    xticklabels=True,
    yticklabels=True,
    ofname=None,
    **kwargs,
):
    """Visualize the dataframe.

    :param df: DataFrame from expression table.
    :param z_score: z_score axis{0, 1}. If None, don't normalize data.
    :param title: gene set name.
    :param outdir: path to save heatmap.
    :param figsize: heatmap figsize.
    :param cmap: matplotlib colormap.
    :param ofname: output file name. If None, don't save figure

    """
    df = zscore(df, axis=z_score)
    df = df.iloc[::-1]
    # If working on commandline, don't show figure
    if hasattr(sys, "ps1") and (ofname is None):
        fig = plt.figure(figsize=figsize)
    else:
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    vmin = np.percentile(df.min(), 2)
    vmax = np.percentile(df.max(), 98)
    matrix = ax.pcolormesh(df.values, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    xstep = _auto_ticks(ax, df.columns.values, 0)
    ystep = _auto_ticks(ax, df.index.values, 1)
    xticks, xlabels = _skip_ticks(df.columns.values, tickevery=xstep)
    yticks, ylabels = _skip_ticks(df.index.values, tickevery=ystep)
    ax.set_ylim([0, len(df)])
    ax.set(xticks=xticks, yticks=yticks)
    ax.set_xticklabels(xlabels if xticklabels else "", fontsize=14, rotation=90)
    ax.set_yticklabels(ylabels if yticklabels else "", fontsize=14)
    ax.set_title("%s\nHeatmap of the Analyzed Geneset" % title, fontsize=20)
    ax.tick_params(
        axis="both", which="both", bottom=False, top=False, right=False, left=False
    )
    # cax=fig.add_axes([0.93,0.25,0.05,0.20])
    cbar = fig.colorbar(matrix, shrink=0.3, aspect=10)
    cbar_title = "z-score" if z_score is not None else ""
    cbar.ax.set_title(cbar_title, loc="left", fontweight="bold")
    for key, spine in cbar.ax.spines.items():
        spine.set_visible(False)
    # cbar = colorbar(matrix)

    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(False)
        # cbar.ax.spines[side].set_visible(False)
    if ofname is not None:
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches="tight", dpi=300)
    return


class GSEAPlot(object):
    def __init__(
        self,
        rank_metric,
        term,
        hits,
        nes,
        pval,
        fdr,
        RES,
        pheno_pos="",
        pheno_neg="",
        figsize=(6, 5.5),
        cmap="seismic",
        ofname=None,
        **kwargs,
    ):
        # center color map at midpoint = 0
        self._norm = MidpointNormalize(midpoint=0)
        # dataFrame of ranked matrix scores
        self._x = np.arange(len(rank_metric))
        self.rankings = np.asarray(rank_metric)
        self.RES = np.asarray(RES)
        self._im_matrix = np.tile(self.rankings, (2, 1))

        self.figsize = figsize
        self.term = term
        self.cmap = cmap
        self.ofname = ofname

        self._pos_label = pheno_pos
        self._neg_label = pheno_neg
        self._zero_score_ind = np.abs(self.rankings).argmin()
        self._z_score_label = "Zero score at " + str(self._zero_score_ind)
        self._hit_indices = hits
        self.module = "tmp" if ofname is None else ofname.split(".")[-2]
        if self.module == "ssgsea":
            self._nes_label = "ES: " + "{:.3f}".format(float(nes))
            self._pval_label = "Pval: invliad for ssgsea"
            self._fdr_label = "FDR: invalid for ssgsea"
        else:
            self._nes_label = "NES: " + "{:.3f}".format(float(nes))
            self._pval_label = "Pval: " + "{:.3e}".format(float(pval))
            self._fdr_label = "FDR: " + "{:.3e}".format(float(fdr))

        # output truetype
        plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
        # in most case, we will have many plots, so do not display plots
        # It's also usefull to run this script on command line.

        # GSEA Plots
        if hasattr(sys, "ps1") and (self.ofname is None):
            # working inside python console, show figure
            self.fig = plt.figure(figsize=self.figsize)
        else:
            # If working on command line, don't show figure
            self.fig = Figure(figsize=self.figsize)
            self._canvas = FigureCanvas(self.fig)

        self.fig.suptitle(self.term, fontsize=16, wrap=True, fontweight="bold")

    def axes_rank(self, rect):
        """
        rect : sequence of float
               The dimensions [left, bottom, width, height] of the new axes. All
               quantities are in fractions of figure width and height.
        """
        # Ranked Metric Scores Plot
        ax1 = self.fig.add_axes(rect, sharex=self.ax)
        if self.module == "ssgsea":
            ax1.fill_between(self._x, y1=np.log(self.rankings), y2=0, color="#C9D3DB")
            ax1.set_ylabel("log ranked metric", fontsize=16, fontweight="bold")
        else:
            ax1.fill_between(self._x, y1=self.rankings, y2=0, color="#C9D3DB")
            ax1.set_ylabel("Ranked list metric", fontsize=16, fontweight="bold")

        ax1.text(
            0.05,
            0.9,
            self._pos_label,
            color="red",
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax1.transAxes,
        )
        ax1.text(
            0.95,
            0.05,
            self._neg_label,
            color="Blue",
            horizontalalignment="right",
            verticalalignment="bottom",
            transform=ax1.transAxes,
        )
        # the x coords of this transformation are data, and the y coord are axes
        trans1 = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
        ax1.vlines(
            self._zero_score_ind,
            0,
            1,
            linewidth=0.5,
            transform=trans1,
            linestyles="--",
            color="grey",
        )

        hap = self._zero_score_ind / max(self._x)
        if hap < 0.25:
            ha = "left"
        elif hap > 0.75:
            ha = "right"
        else:
            ha = "center"
        ax1.text(
            hap,
            0.5,
            self._z_score_label,
            horizontalalignment=ha,
            verticalalignment="center",
            transform=ax1.transAxes,
            fontsize=14,
        )
        ax1.set_xlabel("Rank in Ordered Dataset", fontsize=16, fontweight="bold")
        ax1.spines["top"].set_visible(False)
        ax1.tick_params(
            axis="both", which="both", top=False, right=False, left=False, labelsize=14
        )
        ax1.locator_params(axis="y", nbins=5)
        ax1.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda tick_loc, tick_num: "{:.1f}".format(tick_loc))
        )

    def axes_hits(self, rect):
        """
        rect : sequence of float
               The dimensions [left, bottom, width, height] of the new axes. All
               quantities are in fractions of figure width and height.
        """
        # gene hits
        ax2 = self.fig.add_axes(rect, sharex=self.ax)
        # the x coords of this transformation are data, and the y coord are axes
        trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
        ax2.vlines(self._hit_indices, 0, 1, linewidth=0.5, transform=trans2)
        ax2.spines["bottom"].set_visible(False)
        ax2.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            right=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )

    def axes_cmap(self, rect):
        """
        rect : sequence of float
               The dimensions [left, bottom, width, height] of the new axes. All
               quantities are in fractions of figure width and height.
        """
        # colormap
        ax3 = self.fig.add_axes(rect, sharex=self.ax)
        ax3.imshow(
            self._im_matrix,
            aspect="auto",
            norm=self._norm,
            cmap=self.cmap,
            interpolation="none",
        )  # cm.coolwarm
        ax3.spines["bottom"].set_visible(False)
        ax3.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            right=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )

    def axes_stat(self, rect):
        """
        rect : sequence of float
               The dimensions [left, bottom, width, height] of the new axes. All
               quantities are in fractions of figure width and height.
        """
        # Enrichment score plot

        ax4 = self.fig.add_axes(rect)
        ax4.plot(self._x, self.RES, linewidth=4, color="#88C544")
        ax4.text(0.1, 0.1, self._fdr_label, transform=ax4.transAxes, fontsize=14)
        ax4.text(0.1, 0.2, self._pval_label, transform=ax4.transAxes, fontsize=14)
        ax4.text(0.1, 0.3, self._nes_label, transform=ax4.transAxes, fontsize=14)

        # the y coords of this transformation are data, and the x coord are axes
        trans4 = transforms.blended_transform_factory(ax4.transAxes, ax4.transData)
        ax4.hlines(0, 0, 1, linewidth=0.5, transform=trans4, color="grey")
        ax4.set_ylabel("Enrichment Score", fontsize=16, fontweight="bold")
        # ax4.set_xlim(min(self._x), max(self._x))
        ax4.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            right=False,
            labelbottom=False,
            labelsize=14,
        )
        ax4.locator_params(axis="y", nbins=5)
        # FuncFormatter need two argument, I don't know why. this lambda function used to format yaxis tick labels.
        ax4.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda tick_loc, tick_num: "{:.1f}".format(tick_loc))
        )

        self.ax = ax4

    def add_axes(self):
        """
        Please check matplotlib docs about how to `add_axes` to figure.

        Here is a more flexible way to create a new gseaplot.
        For example, don't show ranking and merge hits and colormap together
        just used:

            self.axes_stat([0.1,0.2,0.8,0.8]) # axes_stat should be called first
            self.axes_cmap([0.1,0.1,0.8,0.1])
            self.axes_hits([0.1,0.1,0.8,0.1])

        """
        self.axes_stat([0.1, 0.5, 0.8, 0.4])
        self.axes_hits([0.1, 0.45, 0.8, 0.05])
        self.axes_cmap([0.1, 0.40, 0.8, 0.05])
        self.axes_rank([0.1, 0.1, 0.8, 0.3])
        # self.fig.subplots_adjust(hspace=0)
        # self.fig.tight_layout()

    def savefig(self, bbox_inches="tight", dpi=300):

        # if self.ofname is not None:
        if hasattr(sys, "ps1") and (self.ofname is not None):
            self.fig.savefig(self.ofname, bbox_inches=bbox_inches, dpi=dpi)
        elif self.ofname is None:
            return
        else:
            self._canvas.print_figure(self.ofname, bbox_inches=bbox_inches, dpi=300)
        return


def gseaplot(
    rank_metric: Iterable,
    term: str,
    hits: List[int],
    nes: float,
    pval: float,
    fdr: float,
    RES: float,
    pheno_pos: str = "",
    pheno_neg: str = "",
    figsize: Tuple[float] = (6, 5.5),
    cmap: str = "seismic",
    ofname: Optional[str] = None,
    **kwargs,
):
    """This is the main function for reproducing the gsea plot.

    :param rank_metric: pd.Series for rankings, rank_metric.values.
    :param term: gene_set name
    :param hits: hits indices of rank_metric.index presented in gene set S.
    :param nes: Normalized enrichment scores.
    :param pval: nominal p-value.
    :param fdr: false discovery rate.
    :param RES: running enrichment scores.
    :param pheno_pos: phenotype label, positive correlated.
    :param pheno_neg: phenotype label, negative correlated.
    :param figsize: matplotlib figsize.
    :param ofname: output file name. If None, don't save figure

    """
    g = GSEAPlot(
        rank_metric,
        term,
        hits,
        nes,
        pval,
        fdr,
        RES,
        pheno_pos,
        pheno_neg,
        figsize,
        cmap,
        ofname,
    )
    g.add_axes()
    g.savefig()


def isfloat(x):
    try:
        float(x)
    except:
        return False
    else:
        return True


def dotplot(
    df,
    column: str = "Adjusted P-value",
    title: str = "",
    cutoff: float = 0.05,
    top_term: int = 10,
    size: float = 10,
    figsize: Tuple[float] = (6, 5.5),
    cmap: str = "viridis_r",
    ofname: Optional[str] = None,
    **kwargs,
):
    """Visualize enrichr results.

    :param df: GSEApy Enrichr DataFrame results.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: terms with 'column' value < cut-off are shown.
    :param top_term: number of enriched terms to show.
    :param ascending: bool, the order of y axis.
    :param size: float, scale the scatter size to get proper visualization.
    :param norm: maplotlib.colors.Normalize object.
    :param legend: bool, whether to show legend.
    :param figsize: tuple, figure size.
    :param cmap: matplotlib colormap
    :param ofname: output file name. If None, don't save figure

    """

    colname = column

    # check if any values in `df[colname]` can't be coerced to floats
    can_be_coerced = df[colname].map(isfloat)
    if np.sum(~can_be_coerced) > 0:
        raise ValueError("some value in %s could not be typecast to `float`" % colname)
    # subset
    df = df[df[colname] <= cutoff]
    if len(df) < 1:
        msg = "Warning: No enrich terms when cutoff = %s" % cutoff
        return msg
    # sorting the dataframe for better visualization
    if colname in ["Adjusted P-value", "P-value"]:
        # df.loc[:, colname] = df[colname].map(float)
        # df = df.assign(logAP=lambda x: -x[colname].apply(np.log10))
        df = df.assign(p_inv=np.log(1 / df[colname]))
        colname = "p_inv"
        cbar_title = r"$Log \frac{1}{P val}$"

    # get top_terms
    df = df.sort_values(by=colname).iloc[-top_term:, :]
    # get scatter area
    temp = df["Overlap"].str.split("/", expand=True).astype(int)
    df = df.assign(Hits_ratio=temp.iloc[:, 0] / temp.iloc[:, 1])
    # make area bigger to better visualization
    # area = df["Hits_ratio"] * plt.rcParams["lines.linewidth"] * 100
    area = np.pi * (df["Hits_ratio"] * size * plt.rcParams["lines.linewidth"]).pow(2)

    # set xaxis values
    xlabel = "Combined Score"
    if "Combined Score" in df.columns:
        x = df["Combined Score"].values
    elif "Odds Ratio" in df.columns:
        x = df["Odds Ratio"].values
        xlabel = "Odds Ratio"
    else:
        # revert back to p_inv
        x = df.loc[colname].values
        xlabel = cbar_title

    # y axis index and values
    # y = [i for i in range(0, len(df))]
    ylabels = df["Term"].values

    # create scatter plot
    if hasattr(sys, "ps1") and (ofname is None):
        # working inside python console, show figure
        fig, ax = plt.subplots(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)

    # scatter colormap range
    colmap = df[colname].round().astype("int")
    vmin = np.percentile(colmap.min(), 2)
    vmax = np.percentile(colmap.max(), 98)
    sc = ax.scatter(
        x=x,
        y=ylabels,
        s=area,
        edgecolors="face",
        c=colmap,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    ax.set_xlabel(xlabel, fontsize=14, fontweight="bold")
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.set_axisbelow(True)  # set grid blew other element
    ax.grid(axis="y")  # zorder=-1.0
    ax.margins(x=0.25)

    # We change the fontsize of minor ticks label
    # ax.tick_params(axis='y', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=14)

    # scatter size legend
    # we use the *func* argument to supply the inverse of the function
    # used to calculate the sizes from above. The *fmt* ensures to string you want
    handles, labels = sc.legend_elements(
        prop="sizes",
        num=3,  # fmt="$ {x:.2f}",
        color="gray",
        func=lambda s: np.sqrt(s / np.pi) / plt.rcParams["lines.linewidth"] / size,
    )
    ax.legend(
        handles,
        labels,
        title="% Path\nis DEG",
        bbox_to_anchor=(1.02, 0.9),
        loc="upper left",
        frameon=False,
    )
    # colorbar
    # cax = fig.add_axes([1.0, 0.20, 0.03, 0.22])
    cbar = fig.colorbar(
        sc,
        shrink=0.2,
        aspect=10,
        anchor=(0.0, 0.2),  # (0.0, 0.2),
        location="right"
        # cax=cax,
    )
    # cbar.ax.tick_params(right=True)
    cbar.ax.set_title(cbar_title, loc="left", fontweight="bold")
    for key, spine in cbar.ax.spines.items():
        spine.set_visible(False)

    ax.set_title(title, fontsize=20, fontweight="bold")
    if ofname is not None:
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches="tight", dpi=300)
        return
    return ax


def ringplot(
    df,
    x: Optional[str] = None,
    column: str = "Adjusted P-value",
    title: str = "",
    cutoff: float = 0.05,
    top_term: int = 10,
    size: float = 10,
    figsize: Tuple[float] = (6, 5.5),
    cmap: str = "viridis_r",
    ofname: Optional[str] = None,
    **kwargs,
):
    """
    when multiple samples exist in the input dataframe, use ringplot instead of heatmap

    :param df: GSEApy Enrichr DataFrame results.
    :param x: column name of x axis data from df.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: terms with 'column' value < cut-off are shown.
    :param top_term: number of enriched terms to show.
    :param sizes: scatter size. Not functional for now
    :param figsize: tuple, figure size.
    :param cmap: matplotlib colormap
    :param ofname: output file name. If None, don't save figure

    """

    colname = column

    # check if any values in `df[colname]` can't be coerced to floats
    can_be_coerced = df[colname].map(isfloat)
    if np.sum(~can_be_coerced) > 0:
        raise ValueError("some value in %s could not be typecast to `float`" % colname)
    # subset
    df = df[df[colname] <= cutoff]
    if len(df) < 1:
        msg = "Warning: No enrich terms when cutoff = %s" % cutoff
        return msg
    # sorting the dataframe for better visualization
    if colname in ["Adjusted P-value", "P-value"]:
        df = df.assign(p_inv=np.log(1 / df[colname]))
        colname = "p_inv"
        cbar_title = r"$Log \frac{1}{P val}$"

    if (x is not None) and (x in df.columns):
        # get top term of each group
        df = (
            df.groupby(x)
            .apply(lambda x: x.sort_values(by=colname).tail(top_term))
            .reset_index(drop=True)
        )
    else:
        df = df.sort_values(by=colname).tail(top_term)

    xlabel = ""
    # set xaxis values, so you could get dotplot
    if (x is not None) and (x in df.columns):
        xlabel = x
    elif "Combined Score" in df.columns:
        xlabel = "Combined Score"
        x = xlabel
    elif "Odds Ratio" in df.columns:
        xlabel = "Odds Ratio"
        x = xlabel
    else:
        # revert back to p_inv
        x = colname
        xlabel = cbar_title

    # get scatter area
    temp = df["Overlap"].str.split("/", expand=True).astype(int)
    df = df.assign(Hits_ratio=temp.iloc[:, 0] / temp.iloc[:, 1])
    # Because the hits_ratio is much too small when being provided as size for ``s``,
    # we normalize it to some useful point sizes, s=0.3*(raito*3)**2
    df = df.assign(
        area=np.pi * (df["Hits_ratio"] * size * plt.rcParams["lines.linewidth"]).pow(2)
    )

    # create scatter plot
    if hasattr(sys, "ps1") and (ofname is None):
        # working inside python console, show figure
        fig, ax = plt.subplots(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)

    # scatter colormap range
    colmap = df[colname].round().astype("int")
    vmin = np.percentile(colmap.min(), 2)
    vmax = np.percentile(colmap.max(), 98)
    # outer ring
    ax.scatter(
        x=x,
        y="Term",
        s=df["area"].max() * 1.5,
        edgecolors="gray",
        c="white",  # colmap,
        data=df,
    )
    sc = ax.scatter(
        x=x,
        y="Term",
        s="area",
        edgecolors="face",
        c=colname,  # colmap,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        data=df,
    )

    # ax.set_xlabel(xlabel, fontsize=14, fontweight="bold")
    # ax.yaxis.set_major_locator(plt.FixedLocator(y))
    # ax.yaxis.set_major_formatter(plt.FixedFormatter(ylabels))
    ax.xaxis.set_tick_params(
        labelsize=14,
        labelrotation=90,
    )
    ax.yaxis.set_tick_params(labelsize=16)
    ax.set_axisbelow(True)  # set grid blew other element
    ax.grid(axis="both")  # zorder=-1.0
    ax.margins(x=0.25)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

    # We change the fontsize of minor ticks label
    # ax.tick_params(axis='y', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=14)

    # scatter size legend
    # we use the *func* argument to supply the inverse of the function
    # used to calculate the sizes from above. The *fmt* ensures to string you want
    handles, labels = sc.legend_elements(
        prop="sizes",
        num=4,  # fmt="$ {x:.2f}",
        color="gray",
        func=lambda s: np.sqrt(s / np.pi) / plt.rcParams["lines.linewidth"] / size,
    )
    ax.legend(
        handles,
        labels,
        title="% Path\nis DEG",
        bbox_to_anchor=(1.02, 0.9),
        loc="upper left",
        frameon=False,
    )
    # colorbar
    # cax = fig.add_axes([1.0, 0.20, 0.03, 0.22])
    cbar = fig.colorbar(
        sc,
        shrink=0.2,
        aspect=10,
        anchor=(0.0, 0.2),  # (0.0, 0.2),
        location="right"
        # cax=cax,
    )
    # cbar.ax.tick_params(right=True)
    cbar.ax.set_title(cbar_title, loc="left", fontweight="bold")
    for key, spine in cbar.ax.spines.items():
        spine.set_visible(False)

    ax.set_title(title, fontsize=20, fontweight="bold")
    if ofname is not None:
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches="tight", dpi=300)
        return
    return ax


def enrichmap(
    df,
):
    """
    Node (inner circle) size corresponds to the number of genes in dataset 1 within the geneset
    Colour of the node (inner circle) corresponds to the significance of the geneset for dataset 1.
    Edge size corresponds to the number of genes that overlap between the two connected genesets.
    Green edges correspond to both datasets when it is the only colour edge.
    When there are two different edge colours, green corresponds to dataset 1 and blue corresponds to dataset 2.
    """
    pass


def barplot(
    df,
    column="Adjusted P-value",
    title="",
    cutoff=0.05,
    top_term=10,
    figsize=(6.5, 6),
    color="salmon",
    ofname=None,
    **kwargs,
):
    """Visualize enrichr results.

    :param df: GSEApy Enrichr DataFrame results.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title.
    :param cutoff: terms with 'column' value < cut-off are shown.
    :param top_term: number of top enriched terms to show.
    :param figsize: tuple, matplotlib figsize.
    :param color: color for bars.
    :param ofname: output file name. If None, don't save figure

    """
    colname = column
    # check if any values in `df[colname]` can't be coerced to floats
    can_be_coerced = df[colname].map(isfloat)
    if np.sum(~can_be_coerced) > 0:
        raise ValueError("some value in %s could not be typecast to `float`" % colname)
    df.loc[:, colname] = df[colname].map(float)
    df = df[df[colname] <= cutoff]
    if len(df) < 1:
        msg = "Warning: No enrich terms using library %s when cutoff = %s" % (
            title,
            cutoff,
        )
        return msg
    if colname in ["Adjusted P-value", "P-value"]:
        df = df.assign(logAP=lambda x: -x[colname].apply(np.log10))
        colname = "logAP"

    dd = df.sort_values(by=colname).iloc[-top_term:, :]
    # create bar plot
    if hasattr(sys, "ps1") and (ofname is None):
        # working inside python console, show (True) figure
        fig = plt.figure(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    bar = dd.plot.barh(x="Term", y=colname, color=color, alpha=0.75, fontsize=16, ax=ax)

    if column in ["Adjusted P-value", "P-value"]:
        xlabel = "-log$_{10}$(%s)" % column
    else:
        xlabel = column
    bar.set_xlabel(xlabel, fontsize=16, fontweight="bold")
    bar.set_ylabel("")
    bar.set_title(title, fontsize=24, fontweight="bold")
    bar.xaxis.set_major_locator(MaxNLocator(integer=True))
    bar.legend_.remove()
    adjust_spines(ax, spines=["left", "bottom"])

    if hasattr(sys, "ps1") and (ofname is not None):
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches="tight", dpi=300)
        return
    elif ofname is None:
        return
    elif not hasattr(sys, "ps1") and (ofname is not None):
        canvas.print_figure(ofname, bbox_inches="tight", dpi=300)
        return
    return ax


def adjust_spines(ax, spines):
    """function for removing spines and ticks.

    :param ax: axes object
    :param spines: a list of spines names to keep. e.g [left, right, top, bottom]
                    if spines = []. remove all spines and ticks.

    """
    for loc, spine in ax.spines.items():
        if loc in spines:
            # spine.set_position(('outward', 10))  # outward by 10 points
            # spine.set_smart_bounds(True)
            continue
        else:
            spine.set_color("none")  # don't draw spine

    # turn off ticks where there is no spine
    if "left" in spines:
        ax.yaxis.set_ticks_position("left")
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if "bottom" in spines:
        ax.xaxis.set_ticks_position("bottom")
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
