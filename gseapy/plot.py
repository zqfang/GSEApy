# -*- coding: utf-8 -*-

import numpy as np
import logging, sys
from matplotlib.colors import Normalize
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
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


def z_score(data2d, axis=0):
    """Standarize the mean and variance of the data axis Parameters.

    :param data2d: DataFrame to normalize.
    :param axis: int, Which axis to normalize across. If 0, normalize across rows, 
                  if 1, normalize across columns.

    
    :Returns: Normalized DataFrame. Noramlized data with a mean of 0 and variance of 1 
              across the specified axis.

    """
    if axis == 1:
        z_scored = data2d
    else:
        z_scored = data2d.T

    z_scored = (z_scored - z_scored.mean()) / z_scored.std()

    if axis == 1:
        return z_scored
    else:
        return z_scored.T


def heatmap(df, term, outdir, axis=0, figsize=(5,5), format='png'):
    """Visualize the dataframe. 
    
    :param df: DataFrame from expression table.
    :param term: gene set name.
    :param outdir: path to save heatmap.
    :param axis: z_score axis.
    :param figsize: heatmap figsize.
    :param format: Matplotlib supported figure formats.
     
    """
    df = z_score(df, axis=axis)
    df = df.iloc[::-1]
    # Get the positions and used label for the ticks
    nx, ny = df.T.shape
    xticks = np.arange(0, nx, 1) + .5
    yticks = np.arange(0, ny, 1) + .5

    #If working on commandline, don't show figure
    fig = Figure(figsize=figsize)
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)        
    vmin = np.percentile(df.min(), 2)
    vmax =  np.percentile(df.max(), 98)
    matrix = ax.pcolormesh(df.values, cmap=plt.cm.RdBu_r, vmin=vmin, vmax=vmax)
    ax.set_ylim([0,len(df)])
    ax.set(xticks=xticks, yticks=yticks)
    ax.set_xticklabels(df.columns.values, fontsize=18, rotation=90)
    ax.set_yticklabels(df.index.values,  fontsize=18)
    ax.set_title("%s\nHeatmap of the Analyzed GeneSet"%term, fontsize=24)
    ax.tick_params(axis='both', which='both', bottom='off', top='off', 
                   right='off', left='off')
    

    
    #fig.colorbar(matrix, ax=ax)
    cax=fig.add_axes([0.93,0.25,0.05,0.20])
    cbar = fig.colorbar(matrix, cax=cax)
    cbar.ax.tick_params(axis='both', which='both', bottom='off', top='off', 
                        right='off', left='off')
    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(False)
        cbar.ax.spines[side].set_visible(False)
    #cbar.ax.set_title('',loc='left')
    term = term.replace('/','_').replace(":","_")
    canvas.print_figure("{a}/{b}.heatmap.{c}".format(a=outdir, b=term, c=format),
                        bbox_inches='tight')

def gsea_plot(rank_metric, enrich_term, hit_ind, nes, pval, fdr, RES,
              phenoPos, phenoNeg, figsize, format, outdir, module):
    """This is the main function for reproducing the gsea plot.
    
    :param rank_metric: rankings, rank_metric['rank'].values.
    :param enrich_term: gene_set name
    :param hit_ind: hit indexs of rank_metric['gene_name'] presented in gene set S.
    :param nes: Normalized enrichment scores.
    :param pval: nominal p-value.
    :param fdr: false discoveray rate.
    :param RES: ranking enrichment scores of all genes in rank_metric['gene_name'].
    :param phenoPos: phenotype lable, positive correlated.
    :param phenoNeg: phenotype lable, negative correlated.
    :param figsize: matplotlib figsize.
    :return: 
    """    
    # plt.style.use('classic')
    # center color map at midpoint = 0
    norm = _MidpointNormalize(midpoint=0)
    
    #dataFrame of ranked matrix scores     
    x = rank_metric.index.values   
    #figsize = (6,6)
    phenoP_label = phenoPos + ' (Positively Correlated)'
    phenoN_label = phenoNeg + ' (Negatively Correlated)'
    zero_score_ind = np.abs(rank_metric['rank']).argmin()
    z_score_label = 'Zero score at ' + str(zero_score_ind)
    nes_label = 'NES: '+ "{:.3f}".format(float(nes))
    pval_label = 'Pval: '+ "{:.3f}".format(float(pval))
    fdr_label = 'FDR: '+ "{:.3f}".format(float(fdr))   
    im_matrix = rank_metric.ix[:,1:].T

    #in most case, we will have mangy plots, so do not display plots
    #It's also convinient to run this script on command line.         
    #plt.ioff()    
    #GSEA Plots
    gs = plt.GridSpec(16,1)
    #fig = plt.figure(figsize=figsize)
    fig = Figure(figsize=figsize)
    canvas = FigureCanvas(fig)
    #Ranked Metric Scores Plot
    ax1 =  fig.add_subplot(gs[11:])
    ax1.fill_between(x, y1= rank_metric['rank'], y2=0, color='#C9D3DB')
    ax1.set_ylabel("Ranked list metric", fontsize=14)    
    ax1.text(.05, .9, phenoP_label, color='red', horizontalalignment='left', verticalalignment='top',
         transform=ax1.transAxes)
    ax1.text(.95, .05, phenoN_label, color='Blue', horizontalalignment='right', verticalalignment='bottom',
         transform=ax1.transAxes)

    # the x coords of this transformation are data, and the y coord are axes
    trans1 = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax1.vlines(zero_score_ind, 0, 1, linewidth=.5, transform=trans1, linestyles='--', color='grey')
    ax1.text(zero_score_ind, 0.5, z_score_label, horizontalalignment='center', verticalalignment='center',
             transform=trans1)    
    ax1.set_xlabel("Rank in Ordered Dataset", fontsize=14)
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(axis='both', which='both', top='off', right='off', left='off')
    ax1.locator_params(axis='y', nbins=5)    
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  '{:.1f}'.format(tick_loc) ))
    
    # use round method to control float number
    #ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  round(tick_loc, 1) ))
    
    #gene hits
    ax2 = fig.add_subplot(gs[8:10], sharex=ax1)

    # the x coords of this transformation are data, and the y coord are axes
    trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    ax2.vlines(hit_ind, 0, 1,linewidth=.5,transform=trans2)
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(axis='both', which='both', bottom='off', top='off', 
                    labelbottom='off', right='off', left='off', labelleft='off')
    #colormap
    ax3 =  fig.add_subplot(gs[10], sharex=ax1)
    ax3.imshow(im_matrix, aspect='auto', norm=norm, cmap=plt.cm.seismic, interpolation='none') # cm.coolwarm
    ax3.spines['bottom'].set_visible(False)
    ax3.tick_params(axis='both', which='both', bottom='off', top='off', 
                    labelbottom='off', right='off', left='off',labelleft='off')

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
    ax4.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off')
    ax4.locator_params(axis='y', nbins=5)
    # FuncFormatter need two argment, I don't know why. this lambda function used to format yaxis tick labels.
    ax4.yaxis.set_major_formatter(plt.FuncFormatter(lambda tick_loc,tick_num :  '{:.1f}'.format(tick_loc)) )
    
    #fig adjustment
    fig.suptitle(enrich_term, fontsize=16)
    fig.subplots_adjust(hspace=0)
    #fig.tight_layout()
    #plt.close(fig)
    enrich_term = enrich_term.replace('/','_').replace(":","_")
    canvas.print_figure('{0}/{1}.{2}.{3}'.format(outdir, enrich_term, module, format), 
                         bbox_inches='tight', dpi=300,)
    return 

def dotplot(df, cutoff=0.05, figsize=(3.5,6), top_term=10, scale=1):
    """Visualize enrichr results.
    
    :param df: GSEApy DataFrame results. 
    :param cutoff: p-adjust cut-off. 
    :param top_term: number of enriched terms to show.
    :param scale: dotplot point size scale.
    :return:  a dotplot for enrichr terms. 

    """
    
    if 'fdr' in df.columns:
        #gsea results
        df.rename(columns={'fdr':'Adjusted P-value',}, inplace=True)
        df['hits_ratio'] =  df['matched_size'] / df['gene_set_size']
    else:
        #enrichr results
        df['Count'] = df['Overlap'].str.split("/").str[0].astype(int)
        df['Background'] = df['Overlap'].str.split("/").str[1].astype(int)
        df['hits_ratio'] =  df['Count'] / df['Background'] 

    # pvalue cut off
    df = df[df['Adjusted P-value'] <= cutoff]
    
    if len(df) < 1:
        logging.warning("Warning: No enrich terms when cuttoff = %s"%cutoff )
        return None
    #sorting the dataframe for better visualization
    df = df.sort_values(by='Adjusted P-value', ascending=False)
    df = df.head(top_term)
    # x axis values
    padj = df['Adjusted P-value']
    combined_score = df['Combined Score'].round().astype('int')
    x = - padj.apply(np.log10)
    # y axis index and values
    y=  [i for i in range(0,len(df))]
    labels = df.Term.values
    
    area = np.pi * (df['Count'] *scale) **2 
    
    #creat scatter plot
    if hasattr(sys, 'ps1'):
        #working inside python console, show figure
        fig, ax = plt.subplots(figsize=figsize)
    else:
        #If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
    vmin = np.percentile(combined_score.min(), 2)
    vmax =  np.percentile(combined_score.max(), 98)        
    sc = ax.scatter(x=x, y=y, s=area, edgecolors='face', c=combined_score,  
                    cmap = plt.cm.RdBu, vmin=vmin, vmax=vmax)
    ax.set_xlabel("-log$_{10}$(Adjust P-value)", fontsize=16)
    ax.yaxis.set_major_locator(plt.FixedLocator(y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(labels))
    ax.set_yticklabels(labels, fontsize=16)
    #ax.set_ylim([-1, len(df)])
    ax.grid()
       
    #colorbar
    cax=fig.add_axes([0.93,0.20,0.07,0.22])
    cbar = fig.colorbar(sc, cax=cax,)
    cbar.ax.tick_params(right='off')
    cbar.ax.set_title('Com-\nscore',loc='left', fontsize=12)


    
    #for terms less than 3
    if len(df) >= 3:
        
        # find the index of the closest value to the median 
        idx = [area.argmax(), np.abs(area - area.mean()).argmin(), area.argmin()]
        idx = unique(idx)
        x2 = [0]*len(idx)
    else:        
        x2 =  [0]*len(df)
        idx = df.index
    #scale of dots
    ax2 =fig.add_axes([0.93,0.55,0.09,0.06*len(idx)])  
    #s=area[idx]  
    l1 = ax2.scatter([],[], s=10, edgecolors='none')
    l2 = ax2.scatter([],[], s=50, edgecolors='none')
    l3 = ax2.scatter([],[], s=100, edgecolors='none')
    labels = df['Count'][idx]
    leg = ax.legend([l1, l2, l3], labels, nrow=3, frameon=True, fontsize=12,
                     handlelength=2, loc = 8, borderpad = 1.8,
                     handletextpad=1, title='Gene\nRatio', scatterpoints = 1)
    
    #canvas.print_figure('test', bbox_inches='tight')
    return fig 

def barplot(df, cutoff=0.05, figsize=(6.5,6), top_term=10):
    """ barplot for enrichr results"""

    # pvalue cut off
    d = df[df['Adjusted P-value'] <= cutoff]
    
    if len(d) < 1:
        return None
    d = d.assign(logAP = - np.log10(d.loc[:,'Adjusted P-value']).values )
    d = d.sort_values('logAP', ascending=False)
    dd = d.head(top_term).sort_values('logAP')
    fig = Figure(figsize=figsize)
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    bar = dd.plot.barh(x='Term', y='logAP', color="salmon", alpha=0.75, edgecolor='none',fontsize=32, ax=ax)
    bar.set_xlabel("-log$_{10}$ Adjust P-value", fontsize=32)
    bar.set_ylabel("")
    #bar.set_title("Enrichr",fontsize=32)
    bar.legend(loc=4)
    #fig.savefig(png, bbox_inches='tight')
    #fig.savefig(pdf, bbox_inches='tight')
    return fig 
    
def adjust_spines(ax, spines):
    """function for removing spines and ticks.
    
    :param ax: axes object
    :param spines: a list of spines names to keep. e.g [left, right, top, bottom]
                    if spines = []. remove all spines and ticks.
    
    """
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
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
