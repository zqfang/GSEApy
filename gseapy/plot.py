# -*- coding: utf-8 -*-
import operator
import sys
import warnings
from typing import Iterable, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from gseapy.scipalette import SciPalette


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # Note also that we must extrapolate beyond vmin/vmax
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y, left=-np.inf, right=np.inf))

    def inverse(self, value):
        y, x = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.interp(value, x, y, left=-np.inf, right=np.inf)


def zscore(data2d: pd.DataFrame, axis: Optional[int] = 0):
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


class Heatmap(object):
    def __init__(
        self,
        df: pd.DataFrame,
        z_score: Optional[int] = None,
        title: Optional[str] = None,
        figsize: Tuple[float] = (5, 5),
        cmap: Optional[str] = None,
        xticklabels: bool = True,
        yticklabels: bool = True,
        ofname: Optional[str] = None,
        **kwargs,
    ):
        self.title = "Heatmap of the Analyzed Geneset" if title is None else title
        self.figsize = figsize
        self.xticklabels = xticklabels
        self.yticklabels = yticklabels
        self.ofname = ofname

        # scale dataframe
        df = df.astype(float)
        df = zscore(df, axis=z_score)
        df = df.iloc[::-1]
        self._df = df
        self.cbar_title = "Scaled Exp" if z_score is None else "Z-Score"
        self.cmap = cmap
        if cmap is None:
            self.cmap = SciPalette.create_colormap()  # navyblue2darkred
        self._zscore = z_score

    def _skip_ticks(self, labels, tickevery):
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

    def _auto_ticks(self, ax, labels, axis):
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

    def get_ax(self):
        if hasattr(sys, "ps1") and (self.ofname is None):
            fig = plt.figure(figsize=self.figsize)
        else:
            fig = Figure(figsize=self.figsize)
            canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        self.fig = fig
        return ax

    def draw(self):
        df = self._df
        ax = self.get_ax()
        vmin = np.percentile(df, 2)
        vmax = np.percentile(df, 98)
        if self._zscore is None:
            norm = Normalize(vmin=vmin, vmax=vmax)
            cbar_locator = MaxNLocator(nbins=5, integer=True)
        else:
            norm = MidpointNormalize(vmin=vmin, vmax=vmax, vcenter=0)
            cbar_locator = MaxNLocator(nbins=3, symmetric=True)  # symmetric=True
        matrix = ax.pcolormesh(
            df.values,
            cmap=self.cmap,
            norm=norm,
            rasterized=True,
        )
        xstep = self._auto_ticks(ax, df.columns.values, 0)
        ystep = self._auto_ticks(ax, df.index.values, 1)
        xticks, xlabels = self._skip_ticks(df.columns.values, tickevery=xstep)
        yticks, ylabels = self._skip_ticks(df.index.values, tickevery=ystep)
        ax.set_ylim([0, len(df)])
        ax.set(xticks=xticks, yticks=yticks)
        ax.set_xticklabels(
            xlabels if self.xticklabels else "", fontsize=14, rotation=90
        )
        ax.set_yticklabels(ylabels if self.yticklabels else "", fontsize=14)
        ax.set_title(self.title, fontsize=20, fontweight="bold")
        ax.tick_params(
            axis="both", which="both", bottom=False, top=False, right=False, left=False
        )
        # cax=fig.add_axes([0.93,0.25,0.05,0.20])
        cbar = self.fig.colorbar(matrix, shrink=0.3, aspect=10)  #  ticks=[-1, 0, 1]
        cbar.ax.yaxis.set_tick_params(
            color="white", direction="in", left=True, right=True
        )
        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        cbar.locator = cbar_locator  # LinearLocator(3)
        cbar.update_ticks()
        cbar.ax.set_title(self.cbar_title, loc="left", fontweight="bold")
        for key, spine in cbar.ax.spines.items():
            spine.set_visible(False)
        # cbar = colorbar(matrix)

        for side in ["top", "right", "left", "bottom"]:
            ax.spines[side].set_visible(False)
            # cbar.ax.spines[side].set_visible(False)
        return ax


def heatmap(
    df: pd.DataFrame,
    z_score: Optional[int] = None,
    title: str = "",
    figsize: Tuple[float] = (5, 5),
    cmap: Optional[str] = None,
    xticklabels: bool = True,
    yticklabels: bool = True,
    ofname: Optional[str] = None,
    **kwargs,
):
    """Visualize the dataframe.

    :param df: DataFrame from expression table.
    :param z_score: 0, 1, or None. z_score axis{0, 1}. If None, not scale.
    :param title: figure title.
    :param figsize: heatmap figsize.
    :param cmap: matplotlib colormap. e.g. "RdBu_r".
    :param xticklabels: bool, whether to show xticklabels.
    :param xticklabels: bool, whether to show xticklabels.
    :param ofname: output file name. If None, don't save figure

    """

    ht = Heatmap(df, z_score, title, figsize, cmap, xticklabels, yticklabels, ofname)
    ax = ht.draw()
    if ofname is None:
        return ax
    # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
    ht.fig.savefig(ofname, bbox_inches="tight", dpi=300)


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

        # dataFrame of ranked matrix scores
        self._x = np.arange(len(rank_metric))
        self.rankings = np.asarray(rank_metric)
        self.RES = np.asarray(RES)

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
        # center color map at midpoint = 0
        vmin = np.percentile(self.rankings.min(), 2)
        vmax = np.percentile(self.rankings.max(), 98)
        midnorm = MidpointNormalize(vmin=vmin, vcenter=0, vmax=vmax)
        # colormap
        ax3 = self.fig.add_axes(rect, sharex=self.ax)
        ax3.pcolormesh(
            self.rankings[np.newaxis, :],
            rasterized=True,
            norm=midnorm,
            cmap=self.cmap,
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
        ax4.hlines(0, 0, 1, linewidth=1, transform=trans4, color="grey")
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


class DotPlot(object):
    def __init__(
        self,
        df: pd.DataFrame,
        x: Optional[str] = None,
        y: str = "Term",
        hue: str = "Adjusted P-value",
        size_scale: float = 5.0,
        thresh: float = 0.05,
        n_terms: int = 10,
        title: str = "",
        figsize: Tuple[float] = (6, 5.5),
        cmap: str = "viridis_r",
        ofname: Optional[str] = None,
        **kwargs,
    ):
        self.marker = "o"
        if "marker" in kwargs:
            self.marker = kwargs["marker"]
        self.y = y
        self.x = x
        self.hue = str(hue)
        self.colname = str(hue)
        self.figsize = figsize
        self.cmap = cmap
        self.ofname = ofname
        self.scale = size_scale
        self.title = title
        self.n_terms = n_terms
        self.thresh = thresh
        self._df = self.process(df)
        plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})

    def isfloat(self, x):
        try:
            float(x)
        except:
            return False
        else:
            return True

    def process(self, df):
        # check if any values in `df[colname]` can't be coerced to floats
        can_be_coerced = df[self.colname].map(self.isfloat).sum()
        if can_be_coerced < len(df):
            msg = "some value in %s could not be typecast to `float`" % self.colname
            raise ValueError(msg)
        # subset
        mask = df[self.colname] <= self.thresh
        if self.colname in ["Combined Score", "NES", "ES", "Odds Ratio"]:
            mask.loc[:] = True

        df = df.loc[mask]
        if len(df) < 1:
            msg = "Warning: No enrich terms when cutoff = %s" % self.thresh
            raise ValueError(msg)
        self.cbar_title = self.colname
        # clip GSEA lower bounds
        # if self.colname in ["NOM p-val", "FDR q-val"]:
        #     df[self.colname].clip(1e-5, 1.0, inplace=True)
        # sorting the dataframe for better visualization
        if self.colname in ["Adjusted P-value", "P-value", "NOM p-val", "FDR q-val"]:
            # get top_terms
            df = df.sort_values(by=self.colname)
            df[self.colname].replace(
                0, method="bfill", inplace=True
            )  ## asending order, use bfill
            df = df.assign(p_inv=np.log10(1 / df[self.colname].astype(float)))
            self.colname = "p_inv"
            self.cbar_title = r"$\log_{10} \frac{1}{P val}$"

        # get top terms; sort ascending
        if (self.x is not None) and (self.x in df.columns):
            # get top term of each group
            df = (
                df.groupby(self.x)
                .apply(lambda _x: _x.sort_values(by=self.colname).tail(self.n_terms))
                .reset_index(drop=True)
            )
        else:
            df = df.sort_values(by=self.colname).tail(self.n_terms)  # acending
        # get scatter area
        ol = df.columns[df.columns.isin(["Overlap", "Tag %"])]
        temp = (
            df[ol].squeeze(axis=1).str.split("/", expand=True).astype(int)
        )  # axis=1, in case you have only 1 row
        df = df.assign(Hits_ratio=temp.iloc[:, 0] / temp.iloc[:, 1])
        return df

    def get_ax(self):
        """
        setup figure axes
        """
        # create fig
        if hasattr(sys, "ps1") and (self.ofname is None):
            # working inside python console, show figure
            fig = plt.figure(figsize=self.figsize)
        else:
            # If working on commandline, don't show figure
            fig = Figure(figsize=self.figsize)
            _canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        self.fig = fig
        return ax

    def set_x(self):
        """
        set x-axis's value
        """
        x = self.x
        xlabel = ""
        # set xaxis values, so you could get dotplot
        if (x is not None) and (x in self._df.columns):
            xlabel = x
        elif "Combined Score" in self._df.columns:
            xlabel = "Combined Score"
            x = xlabel
        elif "Odds Ratio" in self._df.columns:
            xlabel = "Odds Ratio"
            x = xlabel
        elif "NES" in self._df.columns:
            xlabel = "NES"
            x = xlabel
        else:
            # revert back to p_inv
            x = self.colname
            xlabel = self.cbar_title

        return x, xlabel

    def scatter(self, outer_ring=False):
        """
        build scatter
        """
        # scatter colormap range
        # df = df.assign(colmap=self._df[self.colname].round().astype("int"))
        # make area bigger to better visualization
        # area = df["Hits_ratio"] * plt.rcParams["lines.linewidth"] * 100
        df = self._df.assign(
            area=(
                self._df["Hits_ratio"] * self.scale * plt.rcParams["lines.markersize"]
            ).pow(2)
        )
        colmap = df[self.colname].astype(int)
        vmin = np.percentile(colmap.min(), 2)
        vmax = np.percentile(colmap.max(), 98)
        # vmin = np.percentile(df.colmap.min(), 2)
        # vmax = np.percentile(df.colmap.max(), 98)
        ax = self.get_ax()
        # if self.x is None:
        x, xlabel = self.set_x()
        y = self.y
        # outer ring
        if outer_ring:
            smax = df["area"].max()
            # TODO:
            # Matplotlib BUG: when setting edge colors,
            # there's the center of scatter could not aligned.
            # Instead, I have to add more dots in the plot to get the ring
            blk_sc = ax.scatter(
                x=x,
                y=y,
                s=smax * 1.6,
                edgecolors="none",
                c="black",
                data=df,
                marker=self.marker,
                zorder=0,
            )
            wht_sc = ax.scatter(
                x=x,
                y=y,
                s=smax * 1.3,
                edgecolors="none",
                c="white",
                data=df,
                marker=self.marker,
                zorder=1,
            )
            # data = np.array(rg.get_offsets()) # get data coordinates
        # inner circle
        sc = ax.scatter(
            x=x,
            y=y,
            data=df,
            s="area",
            edgecolors="none",
            c=self.colname,
            cmap=self.cmap,
            vmin=vmin,
            vmax=vmax,
            marker=self.marker,
            zorder=2,
        )
        ax.set_xlabel(xlabel, fontsize=14, fontweight="bold")
        ax.xaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_tick_params(labelsize=16)
        ax.set_axisbelow(True)  # set grid blew other element
        ax.grid(axis="y", zorder=-1)  # zorder=-1.0
        ax.margins(x=0.25)

        # We change the fontsize of minor ticks label
        # ax.tick_params(axis='y', which='major', labelsize=16)
        # ax.tick_params(axis='both', which='minor', labelsize=14)

        # scatter size legend
        # we use the *func* argument to supply the inverse of the function
        # used to calculate the sizes from above. The *fmt* ensures to string you want
        handles, labels = sc.legend_elements(
            prop="sizes",
            num=3,  #
            fmt="{x:.2f}",
            color="gray",
            func=lambda s: np.sqrt(s) / plt.rcParams["lines.markersize"] / self.scale,
        )
        ax.legend(
            handles,
            labels,
            title="% Genes\nin set",
            bbox_to_anchor=(1.02, 0.9),
            loc="upper left",
            frameon=False,
            labelspacing=1.0,
        )
        ax.set_title(self.title, fontsize=20, fontweight="bold")
        self.add_colorbar(sc)
        return ax

    def add_colorbar(self, sc):
        """
        :param sc: matplotlib.Scatter
        """
        # colorbar
        # cax = fig.add_axes([1.0, 0.20, 0.03, 0.22])
        cbar = self.fig.colorbar(
            sc,
            shrink=0.25,
            aspect=10,
            anchor=(0.0, 0.2),  # (0.0, 0.2),
            location="right"
            # cax=cax,
        )
        # cbar.ax.tick_params(direction='in')
        cbar.ax.yaxis.set_tick_params(
            color="white", direction="in", left=True, right=True
        )
        cbar.ax.set_title(self.cbar_title, loc="left", fontweight="bold")
        for key, spine in cbar.ax.spines.items():
            spine.set_visible(False)

    def barh(self, color=None, group=None, ax=None):
        """
        Barplot
        """
        if ax is None:
            ax = self.get_ax()
        x, xlabel = self.set_x()
        bar = self._df.plot.barh(
            x=self.y, y=self.colname, alpha=0.75, fontsize=16, ax=ax
        )
        if self.hue in ["Adjusted P-value", "P-value", "FDR q-val", "NOM p-val"]:
            xlabel = "$- \log_{10}$ (%s)" % self.hue
        else:
            xlabel = self.hue
        bar.set_xlabel(xlabel, fontsize=16, fontweight="bold")
        bar.set_ylabel("")
        bar.set_title(self.title, fontsize=24, fontweight="bold")
        bar.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

        # get default color cycle
        if (not isinstance(color, str)) and isinstance(color, Iterable):
            _colors = list(color)
        else:
            prop_cycle = plt.rcParams["axes.prop_cycle"]
            _colors = prop_cycle.by_key()["color"]
        colors = _colors
        # remove old legend first
        bar.legend_.remove()
        if (group is not None) and (group in self._df.columns):
            num_grp = self._df[group].value_counts(sort=False)
            # set colors for each bar (groupby hue)
            colors = []
            legend_elements = []
            for i, n in enumerate(num_grp):
                # cycle _colors if num_grp > len(_colors)
                c = _colors[i % len(_colors)]
                colors += [c] * n
                ele = Line2D(
                    xdata=[0],
                    ydata=[0],
                    marker="o",
                    color="w",
                    label=num_grp.index[i],
                    markerfacecolor=c,
                    markersize=8,
                )
                legend_elements.append(ele)
            # add custom legend
            ax.legend(
                handles=legend_elements,
                loc="upper left",
                title=group,
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
            )
        # update color of bars
        for j, b in enumerate(ax.patches):
            c = colors[j % len(colors)]
            b.set_facecolor(c)

        # self.adjust_spines(ax, spines=["left", "bottom"])
        for side in ["right", "top"]:
            ax.spines[side].set_visible(False)
        # set ticks
        ax.tick_params(axis="both", which="both", top=False, right=False)
        return ax

    def to_edgelist(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        return two dataframe of nodes, and edges
        """
        num_nodes = len(self._df)
        # build graph
        # G = nx.Graph()
        group_loc = None
        if self.x is not None:
            group_loc = self._df.columns.get_loc(self.x)
        term_loc = self._df.columns.get_loc(self.y)  # "Terms"
        if "Genes" in self._df.columns:
            gene_loc = self._df.columns.get_loc("Genes")
        elif "Lead_genes":
            gene_loc = self._df.columns.get_loc("Lead_genes")
        else:
            raise KeyError("Sorry, could not locate enriched gene list")
        # build graph
        genes = self._df.iloc[:, gene_loc].str.split(";")
        ns_loc = self._df.columns.get_loc("Hits_ratio")
        edge_list = []
        nodes = []
        for i in range(num_nodes):
            nodes.append([i, self._df.iloc[i, term_loc], self._df.iloc[i, ns_loc]])
            if group_loc is not None:
                nodes[-1].append(self._df.iloc[i, group_loc])
            for j in range(i + 1, num_nodes):
                set_i = set(genes.iloc[i])
                set_j = set(genes.iloc[j])
                ov = set_i.intersection(set_j)
                if len(ov) < 1:
                    continue
                jaccard_coefficient = len(ov) / len(set_i.union(set_j))
                overlap_coefficient = len(ov) / min(len(set_i), len(set_j))
                edge = [
                    i,
                    j,
                    self._df.iloc[i, term_loc],
                    self._df.iloc[j, term_loc],
                    jaccard_coefficient,
                    overlap_coefficient,
                    ",".join(ov),
                ]
                edge_list.append(edge)
                # G.add_edge(src,
                # targ,
                # jaccard= jaccard_coefficient,
                # overlap = overlap_coefficient,
                # genes = list(ov))
        edges = pd.DataFrame(
            edge_list,
            columns=[
                "src_idx",
                "targ_idx",
                "src_name",
                "targ_name",
                "jaccard_coef",
                "overlap_coef",
                "overlap_genes",
            ],
        )
        node_c = ["node_idx", "node_name", "node_size"]
        if group_loc is not None:
            node_c += ["node_group"]
        nodes = pd.DataFrame(nodes, columns=node_c)
        return nodes, edges


def dotplot(
    df: pd.DataFrame,
    column: str = "Adjusted P-value",
    group: Optional[str] = None,
    title: str = "",
    cutoff: float = 0.05,
    top_term: int = 10,
    size: float = 5,
    figsize: Tuple[float] = (4, 6),
    cmap: str = "viridis_r",
    ofname: Optional[str] = None,
    xticklabels_rot: Optional[float] = None,
    yticklabels_rot: Optional[float] = None,
    marker: str = "o",
    show_ring: bool = False,
    **kwargs,
):
    """Visualize GSEApy Results.
    When multiple datasets exist in the input dataframe, the `group` argument is your friend.

    :param df: GSEApy DataFrame results.
    :param column: column name in `df` to map the dot colors. Default: Adjusted P-value
    :param group: group by the variable in `df` that will produce categorical scatterplot.
    :param title: figure title
    :param cutoff: terms with `column` value < cut-off are shown. Work only for
                   ("Adjusted P-value", "P-value", "NOM p-val", "FDR q-val")
    :param top_term: number of enriched terms to show.
    :param size: float, scale the dot size to get proper visualization.
    :param figsize: tuple, matplotlib figure size.
    :param cmap: matplotlib colormap for mapping the `column` semantic.
    :param ofname: output file name. If None, don't save figure
    :param marker: the matplotlib.markers. See https://matplotlib.org/stable/api/markers_api.html
    :param show_ring bool: whether to show outer ring.

    :return: matplotlib.Axes. return None if given ofname.
             Only terms with `column` <= `cut-off` are plotted.
    """

    dot = DotPlot(
        df=df,
        x=group,
        y="Term",
        hue=column,
        title=title,
        thresh=cutoff,
        n_terms=int(top_term),
        size_scale=size,
        figsize=figsize,
        cmap=cmap,
        ofname=ofname,
        marker=marker,
        **kwargs,
    )
    ax = dot.scatter(outer_ring=show_ring)

    if xticklabels_rot:
        for label in ax.get_xticklabels():
            label.set_ha("right")
            label.set_rotation(xticklabels_rot)

    if yticklabels_rot:
        for label in ax.get_yticklabels():
            label.set_ha("right")
            label.set_rotation(yticklabels_rot)

    if ofname is None:
        return ax
    dot.fig.savefig(ofname, bbox_inches="tight", dpi=300)


def ringplot(
    df: pd.DataFrame,
    column: str = "Adjusted P-value",
    group: Optional[str] = None,
    title: str = "",
    cutoff: float = 0.05,
    top_term: int = 10,
    size: float = 5,
    figsize: Tuple[float] = (4, 6),
    cmap: str = "viridis_r",
    ofname: Optional[str] = None,
    xticklabels_rot: Optional[float] = None,
    yticklabels_rot: Optional[float] = None,
    marker="o",
    show_ring: bool = True,
    **kwargs,
):
    """ringplot is deprecated, use dotplot instead

    :param df: GSEApy DataFrame results.
    :param group: the old `x`. Group by the variable in `df` that will produce categorical scatterplot.
    :param column: column name in `df` to map the dot colors. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: terms with `column` value < cut-off are shown. Work only for
                   ("Adjusted P-value", "P-value", "NOM p-val", "FDR q-val")
    :param top_term: number of enriched terms to show.
    :param size: float, scale the dot size to get proper visualization.
    :param figsize: tuple, matplotlib figure size.
    :param cmap: matplotlib colormap for mapping the `column` semantic.
    :param ofname: output file name. If None, don't save figure
    :param marker: the matplotlib.markers. See https://matplotlib.org/stable/api/markers_api.html
    :param show_ring bool: whether to show outer ring.

    :return: matplotlib.Axes. return None if given ofname.
             Only terms with `column` <= `cut-off` are plotted.
    """
    warnings.warn("ringplot is deprecated; use dotplot instead", DeprecationWarning, 2)
    if "x" in kwargs:
        warnings.warn("x is deprecated; use group", DeprecationWarning, 2)
        kwargs["group"] = kwargs["x"]
        del kwargs["x"]
    ax = dotplot(df, **kwargs)
    return ax


def barplot(
    df: pd.DataFrame,
    column: str = "Adjusted P-value",
    group: Optional[str] = None,
    title: str = "",
    cutoff: float = 0.05,
    top_term: int = 10,
    figsize: Tuple[float, float] = (4, 6),
    color: Union[str, List[str]] = "salmon",
    ofname: Optional[str] = None,
    **kwargs,
):
    """Visualize GSEApy Results.
    When multiple datasets exist in the input dataframe, the `group` argument is your friend.

    :param df: GSEApy DataFrame results.
    :param column: column name in `df` to map the x-axis data. Default: Adjusted P-value
    :param group: group by the variable in `df` that will produce bars with different colors.
    :param title: figure title.
    :param cutoff: terms with `column` value < cut-off are shown. Work only for
                   ("Adjusted P-value", "P-value", "NOM p-val", "FDR q-val")
    :param top_term: number of top enriched terms grouped by `hue` are shown.
    :param figsize: tuple, matplotlib figsize.
    :param color: color or list of matplotlib.colors. Must be reconigzed by matplotlib.
    :param ofname: output file name. If None, don't save figure

    :return: matplotlib.Axes. return None if given ofname.
             Only terms with `column` <= `cut-off` are plotted.
    """
    dot = DotPlot(
        df=df,
        x=group if group else None,  # x turns into hue in bar
        y="Term",
        hue=column,  # hue turns into x in bar
        title=title,
        thresh=cutoff,
        n_terms=int(top_term),
        figsize=figsize,
        cmap="viridis",  # placeholder only
        ofname=ofname,
    )
    if isinstance(color, str):
        color = [color]
    ax = dot.barh(color=color, group=group)

    if ofname is None:
        return ax
    dot.fig.savefig(ofname, bbox_inches="tight", dpi=300)


def traceplot(
    obj,
    terms: Optional[Union[str, List[str]]] = None,
    pheno_pos: str = "",
    pheno_neg: str = "",
    figsize: Tuple[float] = (6, 4),
    cmap: str = "seismic",
    ofname: Optional[str] = None,
    **kwargs,
):
    """Trace plot for terms

    :param obj: GSEA or Prerank Object.
    :param terms: terms to show in trace plot

    """
    # create bar plot
    if hasattr(sys, "ps1") and (ofname is None):
        # working inside python console, show (True) figure
        fig = plt.figure(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)

    if isinstance(terms, str):
        _terms = [terms]
    elif isinstance(terms, list):
        _terms = terms
    else:
        _terms = list(obj.keys())

    for t in _terms:
        if obj.res2d["Name"].nunique() > 1:
            for name, results in obj.results.items():
                if t in results:
                    RES = results[t]["RES"]
                    ax.plot(range(len(RES)), RES, label=name)
        else:
            results = obj.results
            if t in results:
                RES = results[t]["RES"]
                ax.plot(range(len(RES)), RES)
    ax.axhline(0, linewidth=1, linestyle="dashed", color="gray")
    ax.legend()
    ax.set_xlabel("Gene list ranking", fontsize=14, fontweight="bold")
    ax.set_ylabel("Enrichment Score", fontsize=14, fontweight="bold")
    if ofname is not None:
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches="tight", dpi=300)
        return
    return ax


def enrichment_map(
    df: pd.DataFrame,
    column: str = "Adjusted P-value",
    group: Optional[str] = None,
    cutoff: float = 0.05,
    top_term: int = 10,
    **kwargs,
):
    """Visualize GSEApy Results.
    Node size corresponds to the percentage of gene overlap in a certain term of interest.
    Colour of the node corresponds to the significance of the enriched terms.
    Edge size corresponds to the number of genes that overlap between the two connected nodes.
    Gray edges correspond to both nodes when it is the only colour edge.
    When there are two different edge colours, red corresponds to positve nodes and blue corresponds to negative nodes.

    :param df: GSEApy DataFrame results.
    :param column: column name in `df` to map the x-axis data. Default: Adjusted P-value
    :param group: group by the variable in `df` that will produce bars with different colors.
    :param title: figure title.
    :param cutoff: terms with `column` value < cut-off are shown. Work only for
                   ("Adjusted P-value", "P-value", "NOM p-val", "FDR q-val")
    :param top_term: number of top enriched terms grouped by `hue` are shown.
    :param figsize: tuple, matplotlib figsize.
    :param color: color or list of matplotlib.colors. Must be reconigzed by matplotlib.
    :param ofname: output file name. If None, don't save figure

    :return: matplotlib.Axes. return None if given ofname.
             Only terms with `column` <= `cut-off` are plotted.
    """
    # c = column
    if column not in df.columns:
        for c in ["Adjusted P-value", "P-value", "FDR q-val", "NOM p-val"]:
            if c in df:
                column = c
                break

    dot = DotPlot(
        df=df,
        x=group,  # x turns into hue in colors of nodes
        y="Term",  # node
        hue=column,  # node size
        thresh=cutoff,
        n_terms=int(top_term),
    )
    return dot.to_edgelist()
