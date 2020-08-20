import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import linregress
import sys, os, argparse, time,getpass,warnings
import itertools as it


Intro='''
    From meta-results remove one contributing study.

    Applicable to: 
        1. Inverse variance weighted meta-analysis
        2. Sample size weighted meta-analysis
    Figures:
        1. Manhattan plots for before and after removing the sub-study
        2. QQ-plots for before and after removing the sub-study
    Notes:
        1. Only common SNPs in the two file will be analyzed. The SNPs that
        exist only in the original meta-analysis results should be added back.
        

    Author: Yunpeng Wang, yunpeng.wng@gmail.com;
            Jiangming Sun, jiangming.sun@med.lu.se
    Data: 1st Aug, 2020

'''

# Bellow used for flip-strand, borrowed from ldsc
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
BASES = COMPLEMENT.keys()
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
        for x in it.product(BASES, BASES) if x[0] != x[1]}
VALID_SNPS = {x for x in map(lambda y: ''.join(y),
    it.product(BASES, BASES)) if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
MATCH_ALLELES = {x for x in map(lambda y: ''.join(y),
    it.product(VALID_SNPS, VALID_SNPS))
    # strand and ref match 
    if ((x[0] == x[2]) and (x[1] == x[3])) or
    # ref match, strand flip
    ((x[0] == COMPLEMENT[x[2]]) and (x[1] ==
    COMPLEMENT[x[3]])) or
    # ref flip, strand match
    ((x[0] == x[3]) and (x[1] == x[2])) or
    # strand and ref flip
    ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}
FLIP_ALLELES = {''.join(x):
        ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
        # strand flip
        ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
        for x in MATCH_ALLELES}


def flip_snps(alleleData):
    '''
    Flip to match with reference

    Input:
    ------
    alleleData, DataFrame contains at least columns
                    mA1, effect alelle of the meta-analysis result
                    mA2, the other alelle of the meta-analysis result
                    sA1, effect alelle of the sub-study result
                    sA2, the other alelle of the sub-study result

    Return:
    ------
    index,      SNPs need to be flipped 
    A DataFrame for alleles that are bad 

    Note:
    -----
    **   
    '''
    alleles = alleleData.mA1+alleleData.mA2+alleleData.sA1+alleleData.sA2
    idx = alleles.map(FLIP_ALLELES)
    badIdx = idx.isnull()
    alleleData.loc[:,'alleles'] = alleles
    tmp_alleles = FLIP_ALLELES.copy()
    bad_allele_flip = {'ATAT': False,
            'ATTA': True,
            'CGCG': False,
            'CGGC': True,
            'TAAT':True,
            'TATA':False,
            'TAAT':True,
            'GCCG': True,
            'GCGC': False}
    tmp_alleles.update(bad_allele_flip)
    return(alleles.map(tmp_alleles),badIdx)

def qqplot(p1, p2, outf, alpha=0.95, n_quantiles=100, log10conv=True, 
        color=['#440154ff', '#fde725ff', '#09e9f9'], fill_dens=[0.1, 0.1, 0.1]):

    '''
    Double QQ plots with confidence interval, adapted from AssocPlot

    Input:
    ------
    p1:             Numpy 1D array for p values of original meta-analysis
    p2:             Numpy 1D array for p values of after-subtraction analysis
    n_quantiles:    number of quantiles to plot
    alpha:          confidence interval
    log10conv:      conversion to -log10(p) for the figure

    Output:
    -------
    save figure
    '''
    pvals = [p1, p2]
    labels = ['Orignial-Meta', 'Subtracted']
    xmax = 0; ymax = 0
    for j in range(len(pvals)):
        # define quantiles positions:
        q_pos = np.concatenate([np.arange(99.)/len(pvals[j]), 
            np.logspace(-np.log10(len(pvals[j]))+2, 0, n_quantiles)])
        # define quantiles in data
        q_data = mquantiles(pvals[j], prob=q_pos, alphap=0, betap=1, 
                limit=(0, 1)) # linear interpolation
        # define theoretical predictions
        q_th = q_pos.copy()
        # evaluate errors
        q_err = np.zeros([len(q_pos),2])
        if np.sum(alpha) > 0:
            for i in range(0, len(q_pos)):
                q_err[i, :] = beta.interval(alpha, len(pvals[j])*q_pos[i], 
                        len(pvals[j]) - len(pvals[j])*q_pos[i])
            q_err[i, q_err[i, :] < 0] = 1e-15
        slope, intercept, r_value, p_value, std_err = linregress(q_th, q_data)
        plt.plot(-np.log10(q_th[n_quantiles-1:]), 
                -np.log10(q_data[n_quantiles-1:]), '-', color=color[j])
        plt.plot(-np.log10(q_th[:n_quantiles]), 
                -np.log10(q_data[:n_quantiles]), '.', color=color[j], 
                label=labels[j])
        xmax = np.max([xmax, -np.log10(q_th[1])])
        ymax = np.max([ymax, -np.log10(q_data[0])])
    plt.fill_between(-np.log10(q_th), -np.log10(q_err[:,0]), 
        -np.log10(q_err[:,1]), color=color[2], 
        alpha=0.5, label='%1.3f CI'%alpha)
    plt.legend(loc=2)
    plt.xlabel('Theoretical -log10 (P)')
    plt.ylabel('Observed -log10 (P)')
    plt.plot([0, 100], [0, 100],'--k')
    plt.xlim([0, np.ceil(xmax)])
    plt.ylim([0, np.ceil(ymax*1.05)])
    plt.tight_layout()
    plt.savefig(outf)
    plt.close()

def sorted_nicely( l ):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def manhattan(p1, p2, chr, pos, outf, chrs_plot=None, chrs_names=None,
               cut = 2, colors = ['#440154ff', '#fde725ff'], top1 = 0,
               top2 = 0, lines = [10, 15], lines_colors = ['g', 'r'],
               zoom = None):
    '''
    Double manhattan plot, adopted from AssocPlot.

    Input:
    ------
    p1:             p-values for the original meta-analysis results 
    p2:             p-values for the subtracted results
    chr:            chromosomes numbers
    pos:            genomic positions
    chrs_plot:      list of chromosomes that should be plotted. 
                        If empty [] all chromosomes will be plotted
    cut:            lower cut (default 2)
    colors:         sequence of colors (default: black/gray)
    top1:           largest -log10(p)or p1 to display 
    top2:           largest -log10(p)or p2 to display 
    lines:          Horizontal lines to plot.
    lines_colors:   Colors for the horizontal lines.
    zoom:           [chromosome, position, range] Zooms into a region.

    Output:
    Generate a figure
    '''

    # Setting things up
    shift=np.array([0.0])
    plt.clf()
    xlabel='Chromosome'
    ylabel='Subtracted   -log10(p-value)   Original-Meta'

    # If chrs_plot is empty, we need to generate a list of chromosomes
    if chrs_plot is None:
        chrs_list = np.unique(chr)
        if isinstance(chrs_list[0], str):
            chrs_list = sorted_nicely(chrs_list)
        else:
            chrs_list.sort()
    else:
        chrs_list = chrs_plot

    # If chrs_names is empty, we need to generate a list of names for chromosomes
    if chrs_names is None:
        chrs_names = [str(chrs_list[i]) for i in range(len(chrs_list))]

    plot_positions = False
    if len(chrs_list) == 1:
        plot_positions = True

    for ii, i in enumerate(chrs_list):
        ax1 = plt.subplot(2,1,1)
        filt = np.where(chr==i)[0]
        x = shift[-1]+pos[filt]
        y = -np.log10(p1[filt])
        plt.plot(x[y>cut], y[y>cut], '.', color=colors[ii % len(colors)])
        shift_f = np.max(x)

        if zoom is not None:
            if zoom[0] == i:
                zoom_shift = zoom[1] + shift[-1]

        plt.subplot(2,1,2)#, sharex=ax1)
        filt = np.where(chr==i)[0]
        x = shift[-1]+pos[filt]
        y = -np.log10(p2[filt])
        plt.plot(x[y>cut], y[y>cut], '.', color=colors[ii % len(colors)])
        shift_m = np.max(x)

        shift = np.append(shift, np.max([shift_f, shift_m]))

        plt.subplot(2,1,1)
        plt.xlim([0, shift[-1]])

        plt.subplot(2,1,2)
        plt.xlim([0, shift[-1]])

    # Defining top boundary of a plot
    if top1 == 0:
        if sys.version_info[0] == 2:
            top1 = np.nanmax(-np.log10(p1))[0]
        elif sys.version_info[0] == 3:
            top1 = np.nanmax(-np.log10(p1))
        else:
            raise 
    if top2 == 0: 
        if sys.version_info[0] == 2:
            top2 = np.nanmax(-np.log10(p1))[0]
        elif sys.version_info[0] == 3:
            top2 = np.nanmax(-np.log10(p1))
        else:
            raise 

    # Setting up the position of labels:
    shift_label = shift[-1]
    shift = (shift[1:]+shift[:-1])/2.
    labels = chrs_names

    # Plotting horizontal lines
    for i, y in enumerate(lines):
        plt.subplot(2,1,1)
        plt.plot([0, shift_label], [y, y], color = lines_colors[i])
        plt.subplot(2,1,2)
        plt.plot([0, shift_label], [y, y], color = lines_colors[i])

    plt.subplot(2,1,1)
    if not plot_positions:
        plt.xticks(shift, labels,fontsize=6,fontweight='bold')
    plt.ylim([cut+0.05, top1])
    plt.setp(plt.gca().get_xticklabels(), visible=False)
    plt.gca().invert_yaxis()
    if not plot_positions:
        plt.xticks(shift)

    plt.gca().spines['right'].set_color('none')
    plt.subplot(2,1,2)
    plt.ylim([cut, top2])
    if not plot_positions:
        plt.xticks(shift, labels,fontsize=6,fontweight='bold',rotation=45)
    plt.ylabel(ylabel)
    plt.gca().yaxis.set_label_coords(-0.065,1.)
    plt.xlabel(xlabel)
    plt.subplots_adjust(hspace=0.00)

    if zoom is not None:
        plt.subplot(2,1,1)
        plt.xlim([zoom_shift-zoom[2], zoom_shift+zoom[2]])
        plt.subplot(2,1,2)
        plt.xlim([zoom_shift-zoom[2], zoom_shift+zoom[2]])
    plt.gca().spines['right'].set_color('none')
    plt.savefig(outf, dpi=300)
    plt.close()

def sub_IVW(mBeta, mVar, sBeta,sVar):
    """
    Performing subtract meta for inverse variance weighted meta-results.
    
    Input:
    ------
        mBeta,  Beta vector of the meta-analysis result
        mVar,   Varance of Beta vector of the meta-analysis result
        sBeta,  Beta vector of the sub-study result
        sVar,   Varance of Beta vector of the sub-study result
    
    Output:
    -------
        Beta,   Beta vector of the subtracted result
        Se,     Standard error for Beta of the subtracted result
        z,      Z-score for Beta of the subtracted result
        p,      p for Beta of the subtracted result
    Note:
    -----
    There are some cases when the size of sub-study is very close to the
    original meta-analysis, ie. very similar variance of effect, we approximate
    such situation.
    In addition, in case the variance from the substudy is smaller than the
    meta-analysis, we ignore those SNPs.

    """

    m_eps = np.finfo(np.float64).resolution
    mVar = np.where(mVar <m_eps, m_eps, mVar)
    sVar = np.where(sVar <m_eps, m_eps, sVar)

    excl_weighted_betas = sBeta/sVar
    tmp = np.abs(1/mVar - 1/sVar)
    varvec = 1/np.where(tmp < m_eps, m_eps, tmp)
    beta = (mBeta/mVar - excl_weighted_betas)*varvec

    se = np.sqrt(varvec)
    z = beta / se
    p = 2 * norm.cdf(-np.abs(z))
    return (beta, se, z, p)

def sub_SSW(mZ, mN, sZ,sN):
    """
    Performing subtract meta for sample size weighted meta-results.
    
    Input:
    ------
        mZ,     Z-score vector of the meta-analysis result
        mN,     Sample size  vector of the meta-analysis result
        sZ,     Z-score vector of the sub-study result
        sN,     Sample size vector of the sub-study result
    
    Output:
    -------
        z,      Z-score vector of the subtracted result
        p,      p value vector of the subtracted result

    """
    denmvec = np.sqrt(mN - sN)
    bad_idx = mN < sN
    denmvec[bad_idx] = np.nan
    z = (mZ * np.sqrt(mN) -  sZ * np.sqrt(sN)) / denmvec
    p = 2 * norm.cdf(-np.abs(z))
    return (z, p)
       
def _inspect_data(args, logf):
    try:
        if not os.path.exists(args.masf):
            logf.write("!!! "+ args.masf + "Not found !\n")
            raise (ValueError, args.masf + "Not found !")
        if not os.path.exists(args.ssf):
            logf.write("!!! "+ args.ssf + "Not found !\n")
            raise (ValueError, args.ssf + "Not found !")
        mas = pd.read_table(args.masf, delim_whitespace=True, nrows=1, 
                na_values=[' ', '#N/A','\\N','N/A','NA','NULL','NaN', 'nan'])
        ss = pd.read_table(args.ssf, delim_whitespace=True, nrows=1,
                na_values=[' ', '#N/A','\\N','N/A','NA','NULL','NaN', 'nan'])
        if args.masSNP not in mas.columns: 
            logf.write("!!! No SNP ID column: "+args.masSNP+" in masf!\n")
            raise (ValueError, "Need SNP ID column for masf !")
        if args.masP not in mas.columns: 
            logf.write("!!! No p value column: "+args.masSNP+" in masf!\n")
            raise (ValueError, "Need p value column for masf !")
        if args.masCHR not in mas.columns: 
            logf.write("!!! No CHR ID column: "+args.masCHR+" in masf!\n")
            raise (ValueError, "Need CHR number column for masf !")
        if args.masPOS not in mas.columns: 
            logf.write("!!! No genomic position column: "+args.masPOS+" in masf!\n")
            raise (ValueError, "Need genomic position column for masf !")
        if args.masEff not in mas.columns: 
            if args.masOR not in mas.columns: 
                logf.write("!!! No effect column: "+args.masSNP+" in masf!\n")
                raise (ValueError, "Need effect column for masf !")
            else:
                args.isORm = True
        if args.ssSNP not in ss.columns: 
            logf.write("!!! No SNP ID column: "+args.ssSNP+" in ssf!\n")
            raise (ValueError, "Need SNP ID column for ssf !")
        if args.ssP not in ss.columns: 
            logf.write("!!! No p value column: "+args.ssSNP+" in ssf!\n")
            raise (ValueError, "Need p value column for ssf !")
        if args.ssEff not in ss.columns: 
            if args.ssOR not in ss.columns: 
                logf.write("!!! No effect column: "+args.ssSNP+" in ssf!\n")
                raise (ValueError, "Need effect column for ssf !")
            else:
                args.isORs = True
        if args.IVW:
            logf.write("Using inverse variance meta file: "+
                    args.masf+" and \n"+args.ssf+"\n")
            if (args.masSE not in mas.columns) and (args.masVar not in mas.columns):
                logf.write("!!! No standard error column: "+args.masSE+" in masf!\n")
                logf.write("!!! No variance column: "+args.masVar+" in masf !\n")
                raise (ValueError, "Need standard error or variance for masf !")
            if (args.ssSE not in ss.columns) and (args.ssV not in ss.columns):
                logf.write("!!! No standard error: "+args.ssSE+" in ssf Not found!\n")
                logf.write("!!! No variance : "+args.ssV+" in ssf Not found!\n")
                raise (ValueError, "Need standard error or variance for ssf !")
        else:
            logf.write("Using sample size weighted meta file: "+
                    args.masf+" and \n"+args.ssf+"\n")
            if args.masN not in mas.columns:
                logf.write("!!! No sample size column: "+args.masN+" in masf!\n")
                raise (ValueError, "Need sample size for masf !")
            if args.ssN not in ss.columns:
                logf.write("!!! No sample size column : "+args.ssSE+" in ssf!\n")
                raise (ValueError, "Need sample size column for ssf !")
        if args.flip:
            logf.write("Flip strand using: "+ args.masf+" as reference \n")
            if (args.masA1 not in mas.columns) or (args.masA2 not in mas.columns):
                logf.write("!!! No effect allele column: "+args.masA1+" in masf!\n")
                logf.write("!!! No the other allele column: "+args.masA2+" in masf!\n")
                raise (ValueError, "Need SNP alleles columns for masf !")
            if (args.ssA1 not in ss.columns) or (args.ssA2 not in ss.columns):
                logf.write("!!! No effect allele column: "+args.sA1+" in ssf!\n")
                logf.write("!!! No the other allele column: "+args.ssA2+" in ssf!\n")
                raise (ValueError, "Need SNP alleles columns for ssf !")
        return args
    except:
        raise
        logf.write('Input error')
        logf.close()

def OR_2_Beta(ORvec, varVec):
    """
    Get  Beta from OR and OR scale variance
    ref:
    https://www.andrewheiss.com/blog/2016/04/25/convert-logistic-regression-standard-errors-to-odds-ratios-with-r/
    """
    Beta = np.log(ORvec.astype('float64'))
    betaVarVec = np.true_divide(varVec, np.power(ORvec,2))
    return(Beta, betaVarVec)

def sub_meta(args):
    """
    Performing requested subtract meta analysis.

    """
    logf = open(args.out+".log",'w')
    logf.write(time.ctime() + '\n' + getpass.getuser()+'\n'+"*"*40+'\n\n')
    logf.write(Intro+'\n')
    logf.write("*"*40); logf.write('\nParameters will be used:\n')
    for arg in vars(args):
        logf.write('--'+arg+':   '+str(getattr(args, arg))+"\n")
    logf.write("*"*40); logf.write("\n")
    args.isORm = False; args.isORs = False
    args = _inspect_data(args, logf)
    mas = pd.read_table(args.masf, delim_whitespace=True,  
            na_values=[' ', '#N/A','\\N','N/A','NA','NULL','NaN', 'nan'])
    logf.write("** Read in "+str(mas.shape[0])+" SNPs from \n\t" + args.masf+"\n")
    ss = pd.read_table(args.ssf, delim_whitespace=True, 
            na_values=[' ', '#N/A','\\N','N/A','NA','NULL','NaN', 'nan'])
    logf.write("** Read in "+str(ss.shape[0])+" SNPs from \n\t" + args.ssf+"\n")
    mas.rename(columns={args.masSNP:'SNP', args.masP:'mP', args.masCHR:'mCHR',
        args.masPOS:'mPOS'},inplace=True)
    ss.rename(columns={args.ssSNP:'SNP', args.ssP:'sP'},inplace=True)
    if args.isORm:
        mas.loc[:,'mEff'] = np.log(mas.loc[:, args.masOR].astype('float64'))
    else:
        mas.rename(columns={args.masEff:'mEff'},inplace=True)
    if args.isORs:
        ss.loc[:,'sEff'] = np.log(ss.loc[:,args.ssOR].astype('float64'))
    else:
        ss.rename(columns={args.ssEff:'sEff'},inplace=True)
    if args.flip:
        mas.rename(columns={args.masA1:'mA1', args.masA2:'mA2'},inplace=True)
        mas.loc[:, 'mA1'] = mas.mA1.str.upper()
        mas.loc[:,'mA2'] = mas.mA2.str.upper()
        ss.rename(columns={args.ssA1:'sA1', args.ssA2:'sA2'},inplace=True)
        ss.loc[:, 'sA1'] = ss.sA1.str.upper()
        ss.loc[:,'sA2'] = ss.sA2.str.upper()
    if args.IVW:
        if args.masSE in mas.columns:
            if not args.masisORSE:
                mas.loc[:, 'mVar'] = np.power(mas.loc[:,args.masSE],2)
            else:
                tmp, mas.loc[:, 'mVar'] = OR_2_Beta(mas.loc[:,args.masOR],
                        np.power(mas.loc[:,args.masSE],2))
        else:
            mas.rename(columns={args.masVar:'mVar'},inplace=True)
        if args.ssSE in ss.columns:
            if not args.ssisORSE:
                ss.loc[:, 'sVar'] = np.power(ss.loc[:,args.ssSE],2)
            else:
                tmp,ss.loc[:, 'sVar'] = OR_2_Beta(ss.loc[:,args.ssOR],
                        np.power(ss.loc[:,args.ssSE],2))
        else:
            ss.rename(columns={args.ssV:'sVar'},inplace=True)
    else: 
        if args.masZ in mas.columns:
            mas.rename(columns={args.masZ:'mZ'},inplace=True)
        else:
            if args.masSE not in mas.columns:
                raise(ValueError, 'Need effect size and SE to compute zscore')
        if args.ssZ in ss.columns:
            ss.rename(columns={args.ssZ:'sZ'},inplace=True)
        else:
            if args.ssSE not in ss.columns:
                raise(ValueError, 'Need effect size and SE to compute zscore')
        mas.rename(columns={args.masN:'mN'},inplace=True)
        ss.rename(columns={args.ssN:'sN'},inplace=True)
    validx_mas = (np.isfinite(mas.loc[:, 'mEff'])) & (
            np.isfinite(mas.loc[:, 'mP'])) & (mas.loc[:,'mP'] >0) & (
                    mas.loc[:,'mP'] <=1) 
    mas.loc[~validx_mas,:].to_csv(args.out+".badP_or_badEff_masf", index=False)
    logf.write("** Write bad " + str(np.sum(~validx_mas)
        ) + " SNPs,ie. bad p values or bad effects, in "+ 
            args.masf+"to *.badP_or_badEff_masf\n" )
    mas = mas.loc[validx_mas,:]
    validx_ss = (np.isfinite(ss.loc[:, 'sEff'])) & (np.isfinite(
        ss.loc[:, 'sP'])) & (ss.loc[:,'sP'] >=0) & (ss.loc[:,'sP'] <=1) 
    ss.loc[~validx_ss,:].to_csv(args.out+".badP_or_badEff_ssf", index=False)
    logf.write("** Write bad " + str(np.sum(~validx_ss)
        ) + " SNPs,ie. bad p values or bad effects, in "+ 
            args.ssf+"to *.badP_or_badEff_ssf\n" )
    ss = ss.loc[validx_ss,:]

    dat = pd.merge(mas, ss, on=['SNP'], sort=False)
    if args.flip:
        flip_idx,ambi_snp_idx = flip_snps(dat)
        dat.loc[:,'origsEff'] = dat.sEff
        dat.loc[flip_idx, 'sEff'] = -dat.loc[flip_idx,'origsEff']
        dat.loc[:,'origsA1'] = dat.sA1
        dat.loc[:,'origsA2'] = dat.sA2
        dat.loc[flip_idx, 'sA1'] = dat.loc[flip_idx,'mA2']
        dat.loc[flip_idx, 'sA2'] = dat.loc[flip_idx,'mA1']
    if args.IVW:
        logf.write("** Performing IVW subtract meta for  " + str(dat.shape[0]) + " SNPs\n")
        ivw_eff, ivw_se, ivw_z, ivw_p = sub_IVW(dat.mEff, dat.mVar, dat.sEff, 
                dat.sVar)
        dat.loc[:,'subBeta' ]= ivw_eff
        dat.loc[:,'subP' ]= ivw_p
        dat.loc[:,'subSE' ]= ivw_se
    else:
        logf.write("** Performing sample size weighted subtract meta for  " + str(
            dat.shape[0]) + " SNPs, assuming provide effect is Z-score\n")
        sam_z, sam_p = sub_SSW(dat.mEff/dat.loc[:,args.masSE], dat.mN, 
                dat.sEff/dat.loc[:,args.ssSE], dat.sN)
        dat.loc[:,'subZ' ]= sam_z
        dat.loc[:,'subP' ]= sam_p
    logf.write(
        """\n** Write merged data along with subtracted meta-results to\n\t %s;
        subBeta, Subtracted effect for the allele A1 in orignial meta-analysis results (inverse variance weighted)
        subP,    Subtracted P for the subBeta or subZ
        subSE,   Subtracted SE for the subBeta
        subZ,    Subtracted z-score for the allele A1 in the original meta-analysis results (sample size weighted)
            """ % (args.out+'.csv',))
    dat.to_csv(args.out+".tsv", index=False, sep="\t")
    idx_bad = dat.subP.isnull()
    dat.loc[idx_bad,:].to_csv(args.out+".bad_subP_SNPS.tsv", index=False, sep="\t")
    logf.write("!!! There are "+str(np.sum(idx_bad))+""" SNPs, having subtracted meta P (subP) NAN. This is most likely that the original meta results having a larger SE than that of the sub-study. The same case apply to sample size weighted methods, large N in meta than in substudy. Check *.bad_subP_SNPs.tsv\n""")
    if args.flip:
        dat.loc[ambi_snp_idx,:].to_csv(args.out+".ambi_SNPs.tsv", index=False, sep="\t")
        logf.write("!!! There are "+str(np.sum(ambi_snp_idx))+""" SNPs, having ambiguous strand, ie. AT or  CG as A1/A2 or A2/A1. check *.ambi_SNPs.tsv.\n""")

    dat = dat.loc[~idx_bad,:]
    logf.close()
    qqplot(np.array(dat.mP), np.array(dat.subP), args.out+"_qq.pdf")
    manhattan(dat.mP, dat.subP, dat.mCHR, dat.mPOS,args.out+"_manhattan.png",
            top1=args.top1, top2=args.top2)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser(
            prog="Subtract Meta",
            formatter_class=argparse.RawTextHelpFormatter,
            description=Intro)
    parser.add_argument('--masf', type=str,required=True, 
            help='Meta-analysis result file name')  
    parser.add_argument('--masA1', type=str, default='A1',
            help='Meta-analysis result effect allele column name')  
    parser.add_argument('--masA2', type=str, default='A2',
            help='Meta-analysis result the other allele column name')  
    parser.add_argument('--masEff', type=str, default='Beta',
            help='Meta-analysis result effect (of A1) column name')  
    parser.add_argument('--masOR', type=str, default='OR',
            help='Meta-analysis result Odds ratio (of A1) column name')  
    parser.add_argument('--masisORSE', action='store_true',
            help='Is the Meta-analysis result SE on OR scale, default on ln(OR) scale')  
    parser.add_argument('--masP', type=str, default='P',
            help='Meta-analysis result p value (of A1) column name')  
    parser.add_argument('--masSE', type=str, default='SE',
            help='Meta-analysis result standard error (of Beta) column name')  
    parser.add_argument('--masV', type=str, default='VAR',
            help='Meta-analysis result variance (of Beta) column name')  
    parser.add_argument('--masCHR', type=str, default='CHR',
            help='Chromosome number column name in original Meta-analysis results')  
    parser.add_argument('--masPOS', type=str, default='POS',
            help='Genomic position column name in original Meta-analysis results')  
    parser.add_argument('--masN', type=str, default='N',
            help='Meta-analysis result Sample size column name')  
    parser.add_argument('--masZ', type=str, default='Zscore',
            help='Meta-analysis result Z score column name')  
    parser.add_argument('--masSNP', type=str, default='SNP',
            help='Meta-analysis result SNP column name')  
    parser.add_argument('--ssf', type=str, required=True,
            help='Sub-study result file name')  
    parser.add_argument('--ssA1', type=str, default='A1',
            help='sub-study result effect allele column name')  
    parser.add_argument('--ssA2', type=str, default='A2',
            help='sub-study result the other allele column name')  
    parser.add_argument('--ssEff', type=str, default='Beta',
            help='sub-study result effect (of A1) column name')  
    parser.add_argument('--ssOR', type=str, default='OR',
            help='sub-study result Odds ratio (of A1) column name')  
    parser.add_argument('--ssP', type=str, default='P',
            help='sub-study result p value (of A1) column name')  
    parser.add_argument('--ssSE', type=str, default='SE',
            help='sub-study result standard error (of Beta) column name')  
    parser.add_argument('--ssisORSE', action='store_true',
            help='Is the sub-study result SE on OR scale, default on ln(OR) scale')  
    parser.add_argument('--ssV', type=str, default='VAR',
            help='sub-study result variance (of Beta) column name')  
    parser.add_argument('--ssN', type=str, default='N',
            help='sub-study result Sample size column name')  
    parser.add_argument('--ssZ', type=str, default='Zscore',
            help='sub-study result Z score column name')  
    parser.add_argument('--ssSNP', type=str, default='SNP',
            help='sub-study result SNP column name')  
    parser.add_argument('--top1', type=float, default=0.0,
            help='max -log10(P) for original meta-analysis to plot in the manhattan (default from the data)')  
    parser.add_argument('--top2', type=float, default=0.0,
            help='max -log10(P) for subtracted results to plot in the manhattan (default from the data)')  
    parser.add_argument('--flip', action='store_true', 
            help='whether flip strand, using meta-analysis result as reference')
    parser.add_argument('--noIVW', action='store_true',  
            help='whether meta-analysis result is inverse variance weighted? Otherwise using sample size weighted')  
    parser.add_argument('--out', type=str, required=True,
            help='Result file prefix')  
    args = parser.parse_args()
    if args.noIVW:
        args.IVW = False
    else:
        args.IVW = True
    print(Intro)
    print('Parameters will be used:')
    for arg in vars(args):
        print ('--'+arg+':   ',getattr(args, arg))
    sub_meta(args)
