# -*- coding: utf-8 -*-

import numpy as np
import logging, sys

from matplotlib.colors import Normalize
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from gseapy.parser import unique


class _MidpointNormalize(Normalize):
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
        # nnormalized to mean and std using entire matrix
        # z_scored = (data2d - data2d.values.mean()) / data2d.values.std(ddof=1)
        return data2d

    # if axis == 1:
    #     z_scored = data2d
    # else:
    #     z_scored = data2d.T

    # z_scored = (z_scored - z_scored.mean()) / z_scored.std(ddof=1)
    
    # if axis == 1:
    #     return z_scored
    # else:
    #     return z_scored.T
    return data2d.apply(lambda x: (x-x.mean())/x.std(ddof=1), axis=axis)

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def heatmap(df, z_score=None, title='', figsize=(5,5), cmap='RdBu_r', 
            xticklabels=True, yticklabels=True, ofname=None, **kwargs):
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
    # Get the positions and used label for the ticks
    nx, ny = df.T.shape
    xticks = np.arange(0, nx, 1) + .5
    yticks = np.arange(0, ny, 1) + .5

    # If working on commandline, don't show figure
    if hasattr(sys, 'ps1') and (ofname is None): 
        fig = plt.figure(figsize=figsize)
    else:
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    vmin = np.percentile(df.min(), 2)
    vmax =  np.percentile(df.max(), 98)
    matrix = ax.pcolormesh(df.values, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_ylim([0,len(df)])
    ax.set(xticks=xticks, yticks=yticks)
    ax.set_xticklabels(df.columns.values if xticklabels else '', fontsize=14, rotation=90)
    ax.set_yticklabels(df.index.values if yticklabels else '',  fontsize=14)
    ax.set_title("%s\nHeatmap of the Analyzed Geneset"%title, fontsize=20)
    ax.tick_params(axis='both', which='both', bottom=False, top=False,
                   right=False, left=False)
    # cax=fig.add_axes([0.93,0.25,0.05,0.20])
    # cbar = fig.colorbar(matrix, cax=cax)
    cbar = colorbar(matrix)
    cbar.ax.tick_params(axis='both', which='both', bottom=False, top=False,
                        right=False, left=False)
    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(False)
        cbar.ax.spines[side].set_visible(False)
    # cbar.ax.set_title('',loc='left')

    if ofname is not None: 
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
    return

def gseaplot(rank_metric, term, hits_indices, nes, pval, fdr, RES,
              pheno_pos='', pheno_neg='', figsize=(6,5.5), 
              cmap='seismic', ofname=None, **kwargs):
    """This is the main function for reproducing the gsea plot.

    :param rank_metric: pd.Series for rankings, rank_metric.values.
    :param term: gene_set name
    :param hits_indices: hits indices of rank_metric.index presented in gene set S.
    :param nes: Normalized enrichment scores.
    :param pval: nominal p-value.
    :param fdr: false discovery rate.
    :param RES: running enrichment scores.
    :param pheno_pos: phenotype label, positive correlated.
    :param pheno_neg: phenotype label, negative correlated.
    :param figsize: matplotlib figsize.
    :param ofname: output file name. If None, don't save figure 

    """
    # plt.style.use('classic')
    # center color map at midpoint = 0
    norm = _MidpointNormalize(midpoint=0)

    #dataFrame of ranked matrix scores
    x = np.arange(len(rank_metric))
    rankings = rank_metric.values
    # figsize = (6,6)
    phenoP_label = pheno_pos + ' (Positively Correlated)'
    phenoN_label = pheno_neg + ' (Negatively Correlated)'
    zero_score_ind = np.abs(rankings).argmin()
    z_score_label = 'Zero score at ' + str(zero_score_ind)
    nes_label = 'NES: '+ "{:.3f}".format(float(nes))
    pval_label = 'Pval: '+ "{:.3f}".format(float(pval))
    fdr_label = 'FDR: '+ "{:.3f}".format(float(fdr))
    im_matrix = np.tile(rankings, (2,1))

    # output truetype
    plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42})
    # in most case, we will have many plots, so do not display plots
    # It's also usefull to run this script on command line.

    # GSEA Plots
    gs = plt.GridSpec(16,1)
    if hasattr(sys, 'ps1') and (ofname is None):
        # working inside python console, show figure
        fig = plt.figure(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    # Ranked Metric Scores Plot
    ax1 =  fig.add_subplot(gs[11:])
    module = 'tmp' if ofname is  None else ofname.split(".")[-2]
    if module == 'ssgsea':
        nes_label = 'ES: '+ "{:.3f}".format(float(nes))
        pval_label='Pval: '
        fdr_label='FDR: '
        ax1.fill_between(x, y1=np.log(rankings), y2=0, color='#C9D3DB')
        ax1.set_ylabel("log ranked metric", fontsize=14)
    else:
        ax1.fill_between(x, y1=rankings, y2=0, color='#C9D3DB')
        ax1.set_ylabel("Ranked list metric", fontsize=14)
    ax1.text(.05, .9, phenoP_label, color='red',
             horizontalalignment='left', verticalalignment='top',
             transform=ax1.transAxes)
    ax1.text(.95, .05, phenoN_label, color='Blue',
             horizontalalignment='right', verticalalignment='bottom',
             transform=ax1.transAxes)
    # the x coords of this transformation are data, and the y coord are axes
    trans1 = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    if module != 'ssgsea':
        ax1.vlines(zero_score_ind, 0, 1, linewidth=.5, transform=trans1, linestyles='--', color='grey')
        ax1.text(zero_score_ind, 0.5, z_score_label,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=trans1)
    ax1.set_xlabel("Rank in Ordered Dataset", fontsize=14)
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(axis='both', which='both', top=False, right=False, left=False)
    ax1.locator_params(axis='y', nbins=5)
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  '{:.1f}'.format(tick_loc) ))

    # use round method to control float number
    # ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  round(tick_loc, 1) ))

    # gene hits
    ax2 = fig.add_subplot(gs[8:10], sharex=ax1)

    # the x coords of this transformation are data, and the y coord are axes
    trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    ax2.vlines(hits_indices, 0, 1,linewidth=.5,transform=trans2)
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(axis='both', which='both', bottom=False, top=False,
                    labelbottom=False, right=False, left=False, labelleft=False)
    # colormap
    ax3 =  fig.add_subplot(gs[10], sharex=ax1)
    ax3.imshow(im_matrix, aspect='auto', norm=norm, cmap=cmap, interpolation='none') # cm.coolwarm
    ax3.spines['bottom'].set_visible(False)
    ax3.tick_params(axis='both', which='both', bottom=False, top=False,
                    labelbottom=False, right=False, left=False,labelleft=False)

    # Enrichment score plot
    ax4 = fig.add_subplot(gs[:8], sharex=ax1)
    ax4.plot(x, RES, linewidth=4, color ='#88C544')
    ax4.text(.1, .1, fdr_label, transform=ax4.transAxes)
    ax4.text(.1, .2, pval_label, transform=ax4.transAxes)
    ax4.text(.1, .3, nes_label, transform=ax4.transAxes)

    # the y coords of this transformation are data, and the x coord are axes
    trans4 = transforms.blended_transform_factory(ax4.transAxes, ax4.transData)
    ax4.hlines(0, 0, 1, linewidth=.5, transform=trans4, color='grey')
    ax4.set_ylabel("Enrichment score (ES)", fontsize=14)
    ax4.set_xlim(min(x), max(x))
    ax4.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False)
    ax4.locator_params(axis='y', nbins=5)
    # FuncFormatter need two argument, I don't know why. this lambda function used to format yaxis tick labels.
    ax4.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  '{:.1f}'.format(tick_loc)) )

    # fig adjustment
    fig.suptitle(term, fontsize=16, fontweight='bold')
    fig.subplots_adjust(hspace=0)
    # fig.tight_layout()
    if ofname is not None: 
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
    return


def dotplot(df, column='Adjusted P-value', title='', cutoff=0.05, top_term=10, 
            sizes=None, norm=None, legend=True, figsize=(6, 5.5), 
            cmap='RdBu_r', ofname=None, **kwargs):
    """Visualize enrichr results.

    :param df: GSEApy DataFrame results.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: p-adjust cut-off.
    :param top_term: number of enriched terms to show.
    :param ascending: bool, the order of y axis.
    :param sizes: tuple, (min, max) scatter size. Not functional for now
    :param norm: maplotlib.colors.Normalize object.
    :param legend: bool, whether to show legend.
    :param figsize: tuple, figure size. 
    :param cmap: matplotlib colormap
    :param ofname: output file name. If None, don't save figure 

    """


    colname = column    
    # sorting the dataframe for better visualization
    if colname in ['Adjusted P-value', 'P-value']: 
        df = df[df[colname] <= cutoff]
        if len(df) < 1: 
            msg = "Warning: No enrich terms when cutoff = %s"%cutoff
            return msg
        df = df.assign(logAP=lambda x: - x[colname].apply(np.log10))
        colname='logAP'
    df = df.sort_values(by=colname).iloc[-top_term:,:]
    # 
    temp = df['Overlap'].str.split("/", expand=True).astype(int)
    df = df.assign(Hits=temp.iloc[:,0], Background=temp.iloc[:,1])
    df = df.assign(Hits_ratio=lambda x:x.Hits / x.Background)
    # x axis values
    x = df.loc[:, colname].values
    combined_score = df['Combined Score'].round().astype('int')
    # y axis index and values
    y = [i for i in range(0,len(df))]
    ylabels = df['Term'].values
    # Normalise to [0,1]
    # b = (df['Count']  - df['Count'].min())/ np.ptp(df['Count'])
    # area = 100 * b
    
    # control the size of scatter and legend marker
    levels = numbers = np.sort(df.Hits.unique())
    if norm is None:
        norm = Normalize()
    elif isinstance(norm, tuple):
        norm = Normalize(*norm)
    elif not isinstance(norm, Normalize):
        err = ("``size_norm`` must be None, tuple, "
                "or Normalize object.")
        raise ValueError(err)
    min_width, max_width = np.r_[20, 100] * plt.rcParams["lines.linewidth"]
    norm.clip = True
    if not norm.scaled():
        norm(np.asarray(numbers))
    size_limits = norm.vmin, norm.vmax
    scl = norm(numbers)
    widths = np.asarray(min_width + scl * (max_width - min_width))
    if scl.mask.any():
        widths[scl.mask] = 0
    sizes = dict(zip(levels, widths))
    df['sizes'] = df.Hits.map(sizes)
    area = df['sizes'].values

    # creat scatter plot
    if hasattr(sys, 'ps1') and (ofname is None):
        # working inside python console, show figure
        fig, ax = plt.subplots(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
    vmin = np.percentile(combined_score.min(), 2)
    vmax =  np.percentile(combined_score.max(), 98)
    sc = ax.scatter(x=x, y=y, s=area, edgecolors='face', c=combined_score,
                    cmap=cmap, vmin=vmin, vmax=vmax)

    if column in ['Adjusted P-value', 'P-value']:
        xlabel = "-log$_{10}$(%s)"%column
    else:
        xlabel = column 
    ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
    ax.yaxis.set_major_locator(plt.FixedLocator(y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(ylabels))
    ax.set_yticklabels(ylabels, fontsize=16)
    
    # ax.set_ylim([-1, len(df)])
    ax.grid()
    # colorbar
    cax=fig.add_axes([0.95,0.20,0.03,0.22])
    cbar = fig.colorbar(sc, cax=cax,)
    cbar.ax.tick_params(right=True)
    cbar.ax.set_title('Combined\nScore',loc='left', fontsize=12)

    # for terms less than 3
    if len(df) >= 3:
        # find the index of the closest value to the median
        idx = [area.argmax(), np.abs(area - area.mean()).argmin(), area.argmin()]
        idx = unique(idx)
    else:
        idx = df.index.values
    label = df.iloc[idx, df.columns.get_loc('Hits')]
    
    if legend:
        handles, _ = ax.get_legend_handles_labels()
        legend_markers = []
        for ix in idx: 
            legend_markers.append(ax.scatter([],[], s=area[ix], c='b'))
        # artist = ax.scatter([], [], s=size_levels,) 
        ax.legend(legend_markers, label, title='Hits')
    ax.set_title(title, fontsize=20, fontweight='bold')
    
    if ofname is not None: 
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
    return ax

def barplot(df, column='Adjusted P-value', title="", cutoff=0.05, top_term=10,
            figsize=(6.5,6), color='salmon', ofname=None, **kwargs):
    """Visualize enrichr results.

    :param df: GSEApy DataFrame results.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title.
    :param cutoff: cut-off of the cloumn you've chose.
    :param top_term: number of top enriched terms to show.
    :param figsize: tuple, matplotlib figsize.
    :param color: color for bars.
    :param ofname: output file name. If None, don't save figure    
    
    """

    colname = column   
    if colname in ['Adjusted P-value', 'P-value']: 
        df = df[df[colname] <= cutoff]
        if len(df) < 1: 
            msg = "Warning: No enrich terms using library %s when cutoff = %s"%(title, cutoff)
            return msg
        df = df.assign(logAP = lambda x: - x[colname].apply(np.log10))
        colname = 'logAP' 
    dd = df.sort_values(by=colname).iloc[-top_term:,:]
    # dd = d.head(top_term).sort_values('logAP')
    # create bar plot
    if hasattr(sys, 'ps1') and (ofname is None):
        # working inside python console, show (True) figure
        fig = plt.figure(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    bar = dd.plot.barh(x='Term', y=colname, color=color, 
                       alpha=0.75, fontsize=20, ax=ax)
    
    if column in ['Adjusted P-value', 'P-value']:
        xlabel = "-log$_{10}$(%s)"%column
    else:
        xlabel = column 
    bar.set_xlabel(xlabel, fontsize=20, fontweight='bold')
    bar.set_ylabel("")
    bar.set_title(title, fontsize=28, fontweight='bold')
    bar.xaxis.set_major_locator(MaxNLocator(integer=True))
    bar.legend_.remove()
    adjust_spines(ax, spines=['left','bottom'])

    if ofname is not None: 
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
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
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
