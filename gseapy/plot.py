# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.colors import Normalize



class _MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))






def gsea_plot(rank_metric, enrich_term, hit_ind, nes, pval, fdr, RES,
              phenoPos=None, phenoNeg=None, figsize =(6.5,6), **kwarg):
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
    :return: fig object of gsea plot.
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
    plt.ioff()    
    #GSEA Plots
    gs = plt.GridSpec(16,1)
    fig = plt.figure(figsize=figsize)
    #Ranked Metric Scores Plot
    ax1 =  fig.add_subplot(gs[11:])
    ax1.fill_between(x, y1= rank_metric['rank'], y2=0, color='#C9D3DB')
    ax1.set_ylabel("Ranked list metric",fontsize=14)    
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
                    labelbottom='off', right='off', left='off',labelleft='off')
    #colormap
    ax3 =  fig.add_subplot(gs[10],sharex=ax1)
    ax3.imshow(im_matrix, aspect='auto', norm=norm, cmap=plt.cm.seismic, interpolation='none') # cm.coolwarm
    ax3.spines['bottom'].set_visible(False)
    ax3.tick_params(axis='both', which='both', bottom='off', top='off', 
                    labelbottom='off', right='off', left='off',labelleft='off')

    # Enrichment score plot
    ax4 = fig.add_subplot(gs[:8],sharex=ax1)
    ax4.plot(x,RES,linewidth=4,color ='#88C544')
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
    plt.close(fig)
    
    return fig

def dotplot(df, cutoff=0.05, figsize=(3,6)):
    """Visualize enrichr results
    :param df: GSEApy DataFrame results 
    :return: dotplot
    """
    # pvalue cut off
    df = df[df['Adjusted P-value'] <= cutoff]
    
    #sorting the dataframe for better visualization
    df = df.sort_values(by='Adjusted P-value', ascending=False)
    
    # x axis values
    padj = df['Adjusted P-value']
    x = - padj.apply(np.log10)
    # y axis index and values
    y=  [i for i in range(0,len(df))]
    labels = df.Term.values
    
    #gene ratio      
    df['Count'] = df['Overlap'].str.split("/").str[0].astype(int)
    df['Background'] = df['Overlap'].str.split("/").str[1].astype(int)
      
    hits_ratio =  df['Count'] / df['Background'] 
    area = np.pi * (hits_ratio *50) **2 
    
    #creat scatter plot
    fig, ax = plt.subplots(figsize=figsize)
    sc = ax.scatter(x=x, y=y, s=area, edgecolors='face', c = padj,  
                    cmap = plt.cm.bwr_r,vmin=padj.min(), vmax=padj.max())
    ax.set_xlabel("-log$_{10}$(Adjust P-value)")
    ax.yaxis.set_major_locator(plt.FixedLocator(y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(labels))
    ax.set_ylim([-1, len(df)])
    ax.grid()
    

    
    #colorbar
    cax=fig.add_axes([0.93,0.25,0.05,0.20])
    cbar = fig.colorbar(sc, cax=cax,)
    cbar.ax.tick_params(right='off')
    cbar.ax.set_title('Padj',loc='left')

    #scale of dots
    ax2 =fig.add_axes([0.93,0.55,0.05,0.12])
  
     
    # find the index of the closest value to the median 
    idx = [area.argmax(), np.abs(area - area.median()).argmin(), area.argmin()]

    ax2.scatter(x=[0,0,0], y=y[:3],s=area[idx], c='black', edgecolors='face')
    for i, index in enumerate(idx):
        ax2.text(x=0.8, y=y[i], s=hits_ratio[index].round(2), verticalalignment='center', horizontalalignment='left')
    ax2.set_title("Ratio",loc='left')
    
    #turn off all spines and ticks
    ax2.axis('off')

    #plt.tight_layout()
    
    plt.show()
    
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
