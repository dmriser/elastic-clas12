#!/usr/bin/env python 

import argparse 
import numpy as np

# Trick for docker install to run headless
# and never use plt.show() 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

from scipy.optimize import minimize 

from array import array 
from ROOT import (TH1F, TH2F, TH1D, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors, TLine)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def numpify(histo):
    """ TH1F to np.arrays. """

    if type(histo) not in [TH1F, TH1D]:
        raise NotImplementedError('Can not numpify type {}'.format(type(histo)))

    nbins = histo.GetNbinsX()

    # Setup output 
    x_lows = np.zeros(nbins)
    x_highs = np.zeros(nbins)
    values = np.zeros(nbins)
    errors = np.zeros(nbins)
    
    for i in range(1, nbins + 1):
        x_lows[i-1] = histo.GetBinLowEdge(i)
        x_highs[i-1] = histo.GetBinLowEdge(i) + histo.GetBinWidth(i)
        values[i-1] = histo.GetBinContent(i)
        errors[i-1] = histo.GetBinError(i)

    return x_lows, x_highs, values, errors

def chi2(data, theory, err):
    return np.sum((data-theory)**2 / (0.00001 + err**2)) 

def model(x,p):
    return p[0] * np.exp( -0.5 * (x - p[1])**2 / p[2]**2 )

def scipy_fit_slice(x, y, err, bounds):        
    
    p0 = np.random.uniform(-1, 1, 3)    
    p0[0] = np.max(y)
    p0[1] = np.mean(x*y)
    p0[2] = np.std(x*y)
    
    res = minimize(lambda p: chi2(y, model(x,p), err), x0=p0, bounds=bounds, method='Nelder-Mead')
    return res.x

def fit_slices(histo, x_range, x_bin_step, y_fit_range):

    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)

        fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus')
        fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
        fit.SetName(histo.GetTitle() + '_fit{}'.format(i))

        if y_fit_range:
            fit.SetParameter(1, 0.5 * (y_fit_range[0] + y_fit_range[1]))

        projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])

        slices.append(projec)
        fits.append(fit)

        x_low = histo.GetXaxis().GetBinCenter(x_bin)
        x_high = histo.GetXaxis().GetBinCenter(x_bin + x_bin_step)
        x_values.append(0.5 * (x_high + x_low))

    means = array('d')
    means_err = array('d')
    stds = array('d')
    stds_err = array('d')
    zeros = array('d')

    for f in fits:
        means.append(f.GetParameter(1))
        means_err.append(f.GetParError(1))
        stds.append(f.GetParameter(2))
        stds_err.append(f.GetParError(2))
        zeros.append(0.0)

    graph = TGraphErrors(len(x_values), x_values, means, zeros, stds)
    graph.SetName('g_' + histo.GetName())

    # return graph, slices, fits
    return np.array(x_values), np.array(means), np.array(stds), slices, fits
    
def plot_fits_mpl(histos1, histos2, config1, config2,
                  x_range, x_bin_step, title_formatter, y_fit_range,
                  save_name, y_range=None, title=None,
                  xtitle=None, ytitle=None, max_errorbar=0.5,
                  hline=None, x_shift=None):
    
    fig = plt.figure(figsize=(16,12))

    # Plot options for each type. 
    opts1 = {'marker':'o', 'linestyle':'', 'color':'k'}
    opts2 = {'marker':'o', 'linestyle':'', 'color':'red'}
    
    for i in range(1,7):

        ax = fig.add_subplot(2, 3, i)

        # Get histogram slices for plotting. 
        tit = title_formatter.format(i)
        print('Fitting: ', tit)

        # Somebody wants that we should get the histograms ourselves
        if isinstance(histos1, dict) and isinstance(histos2, dict):
            x1, mu1, sig1, slices1, fits1 = fit_slices(histos1.get(tit, default_histo2d), x_range, x_bin_step, y_fit_range)
            x2, mu2, sig2, slices2, fits2 = fit_slices(histos2.get(tit, default_histo2d), x_range, x_bin_step, y_fit_range)
            
        # Somebody passed actual histograms, should be six of them 
        elif isinstance(histos1, list) and isinstance(histos2, list):
            x1, mu1, sig1, slices1, fits1 = fit_slices(histos1[i-1], x_range, x_bin_step, y_fit_range)
            x2, mu2, sig2, slices2, fits2 = fit_slices(histos2[i-1], x_range, x_bin_step, y_fit_range)

        x1, mu1, sig1 = remove_bad_points(x1, mu1, sig1, max_errorbar)
        x2, mu2, sig2 = remove_bad_points(x2, mu2, sig2, max_errorbar)

        if x_shift:
            x2 += 0.5 * (x2[1]-x2[0])
        
        label1 = 'Sector {} ({})'.format(i, config1)
        label2 = 'Sector {} ({})'.format(i, config2)

        # Draw things 
        ax.errorbar(x1, mu1, sig1, label=label1, **opts1)
        ax.errorbar(x2, mu2, sig2, label=label2, **opts2)

        if y_range:
            ax.set_ylim(y_range)

        # Add center line
        if hline:
            ax.axhline(hline, linestyle='--', linewidth=1, color='k', alpha=0.85)

        # Add a grid
        ax.grid(alpha=0.2)

        # Legend
        ax.legend(frameon=False)
        
        if title:
            ax.set_title(title)

        if xtitle:
            ax.set_xlabel(xtitle)

        if ytitle:
            ax.set_ylabel(ytitle)
            
    fig.tight_layout()
    fig.savefig(save_name, bbox_inches='tight')

    output_slice_pdf(histos1, title_formatter, 'data', x_range, x_bin_step, y_fit_range)
    output_slice_pdf(histos2, title_formatter, 'sim', x_range, x_bin_step, y_fit_range)    

def output_slice_pdf(histos, title_formatter, sub_title, x_range, x_bin_step, y_fit_range):
    """ Debugging of the fits for slices. """

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextSize(0.03)
    
    slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
    slice_pdfname = title_formatter.split('_{}')[0] + '_' + sub_title + '_slices.pdf'
    slice_can.Print(slice_pdfname + '[')
    for i in range(1,7):
        title = title_formatter.format(i)

        if isinstance(histos, dict):
            x, mu, sig, slices, fits = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range)
        elif isinstance(histos, list):
            x, mu, sig, slices, fits = fit_slices(histos[i-1], x_range, x_bin_step, y_fit_range)
            
        nrows = 5
        ncols = int(np.ceil(len(slices) / nrows) + 1)
        slice_can.Clear()
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices, fits)):
            slice_can.cd(j+1)
            s.Draw()
            lab.DrawLatex(0.15, 0.88, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(f.GetParameter(1), f.GetParameter(2)))

        slice_can.Print(slice_pdfname)
    slice_can.Print(slice_pdfname + ']')
    
    
def remove_bad_points(x, mu, sig, max_errorbar):
    condition = np.logical_and(mu != 0, sig < max_errorbar)
    idx = np.where(condition)[0]
    return x[idx], mu[idx], sig[idx]

def plot_sector_page_single(canvas, histos, title_formatter, save_name, xtitle=None,
                            ytitle=None, title=None, x_range=None, log=False, hline=None):
    """ Plot one histogram for each sector. """
    
    label = TLatex()
    label.SetNDC()
    label.SetTextSize(0.045)
    
    canvas.Clear()
    canvas.Divide(3,2)

    root_garbage_can = [] 
    for i in range(1,7):
        canvas.cd(i)

        if x_range:
            histos.get(title_formatter.format(i), default_histo2d).GetXaxis().SetRangeUser(x_range[0], x_range[1])

        histos.get(title_formatter.format(i), default_histo2d).Draw('colz')
        if log:
            gPad.SetLogz()
        
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.45, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.0325, 0.65, ytitle)
            label.SetTextAngle(0)

        if hline:
            xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax()
            line = TLine(xmin, hline, xmax, hline)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
            
    canvas.Print(save_name)
    

def plot_sector_page(canvas, histos1, histos2, config1, config2, title_formatter,
                     save_name, xtitle=None, ytitle=None, title=None, x_range=None):
    """ Compare two distributions in one plot, using root. """
    
    label = TLatex()
    label.SetNDC()
    label.SetTextSize(0.045)
    
    canvas.Clear()
    canvas.Divide(3,2)

    for i in range(1,7):
        canvas.cd(i)

        #_, __, values1, errors1 = numpify(histos1.get(title_formatter.format(i), default_histo))
        #_, __, values2, errors2 = numpify(histos2.get(title_formatter.format(i), default_histo))

        #scale = np.sum(values1 * values2) / np.sum(values2**2)
        
        histos1.get(title_formatter.format(i), default_histo).SetLineColor(1)
        histos1.get(title_formatter.format(i), default_histo).Scale(1 / histos1.get(title_formatter.format(i), default_histo).GetMaximum())
        
        histos2.get(title_formatter.format(i), default_histo).SetLineColor(99)
        histos2.get(title_formatter.format(i), default_histo).Scale(1 / histos2.get(title_formatter.format(i), default_histo).GetMaximum())
        #histos2.get(title_formatter.format(i), default_histo).Scale(scale)

        hmax = max(histos1.get(title_formatter.format(i), default_histo).GetMaximum(),
                   histos2.get(title_formatter.format(i), default_histo).GetMaximum())
        
        histos1.get(title_formatter.format(i), default_histo).SetMaximum(hmax * 1.1)
        histos2.get(title_formatter.format(i), default_histo).SetMaximum(hmax * 1.1)

        if x_range:
            histos1.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
            histos2.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
        
        histos1.get(title_formatter.format(i), default_histo).Draw('hist')
        histos2.get(title_formatter.format(i), default_histo).Draw('histsame')

        # Annotate the figure with the configuration names.
        caption1 = '#color[1]({0})'.format(config1).replace('(', '{').replace(')', '}')
        caption2 = '#color[99]({0})'.format(config2).replace('(', '{').replace(')', '}')

        label.DrawLatex(0.75, 0.86, caption1)
        label.DrawLatex(0.75, 0.825, caption2)
        
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.45, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.02, 0.65, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)
    
if __name__ == '__main__':

    # Parse command line arguments. 
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--data_file', required=True)
    ap.add_argument('-s', '--sim_file', required=True)
    ap.add_argument('-o', '--output_prefix', required=True)
    args = ap.parse_args()

    # Setup files
    files = {}
    files['data'] = TFile(args.data_file)
    files['sim'] = TFile(args.sim_file)

    output_pdfname = args.output_prefix + '.pdf'

    # Load histograms
    histos = {}
    for config_type, file in files.items():
        histos[config_type] = load_histos(file)
        print(histos[config_type].keys())

    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
                  x_range=[0.40,2.00], y_range=[-1.0,1.0], x_bin_step=2, title_formatter='histos_p_pro_dp_pro_from_angles_CTOF_{}',
                  save_name='p_pro_dp_pro_fit_rad_{}.png'.format(args.output_prefix),
                  title='Proton Momentum Resolution', xtitle='$P_p$', ytitle='$\Delta P_{p}$',
                  max_errorbar = 0.8, y_fit_range=[-0.8, 0.8], hline=0.00001, x_shift=True
    ) 
 
    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
                  x_range=[1.5,4.00], y_range=[-1.0,1.0], x_bin_step=2, title_formatter='histos_p_ele_dp_ele_from_angles_CTOF_{}',
                  save_name='p_ele_dp_ele_fit_rad_{}.png'.format(args.output_prefix),
                  title='Electron Momentum Resolution', xtitle='$P_e$', ytitle='$\Delta P_{e}$',
                  max_errorbar = 0.8, y_fit_range=[-1.4, 1.4], hline=0.00001, x_shift=True
    ) 

    """

    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
                  x_range=[7,30], y_range=[-3.0,3.0], x_bin_step=4, title_formatter='histos_theta_ele_dtheta_ele_CTOF_{}',
                  save_name='theta_ele_dtheta_ele_fit_rad_{}.png'.format(args.output_prefix),
                  title='Electron Angular Resolution', xtitle='$\theta_e$', ytitle='$\Delta \theta_{e}$',
                  max_errorbar = 0.8, y_fit_range=[-3, 3], hline=0.00001, x_shift=False
    ) 

    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
                  x_range=[50,70], y_range=[-3.0,3.0], x_bin_step=4, title_formatter='histos_theta_pro_dtheta_pro_CTOF_{}',
                  save_name='theta_pro_dtheta_pro_fit_rad_{}.png'.format(args.output_prefix),
                  title='Proton Angular Resolution', xtitle='$\theta_p$', ytitle='$\Delta \theta_{p}$',
                  max_errorbar = 0.8, y_fit_range=[-3, 3], hline=0.00001, x_shift=False
    ) 
    """
    
    # ------------------------------------------------------
    # Plot some distribution compare 
    # ------------------------------------------------------
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)

    can = TCanvas('can', 'can', 1100, 800)

    plot_sector_page(can, histos['data'], histos['sim'],
                     config1='Data', config2='Sim', title_formatter='histos_p_pro_CTOF_{}',
                     save_name='histos_p_proton_ctof_compare_rad.pdf',
                     xtitle='P_{p} (GeV/c)', title='P_{p}')

    plot_sector_page(can, histos['data'], histos['sim'],
                     config1='Data', config2='Sim', title_formatter='histos_p_pro_FTOF_{}',
                     save_name='histos_p_proton_ftof_compare_rad.pdf',
                     xtitle='P_{p} (GeV/c)', title='P_{p}')

    plot_sector_page(can, histos['data'], histos['sim'],
                     config1='Data', config2='Sim', title_formatter='histos_p_ele_CTOF_{}',
                     save_name='histos_p_electron_ctof_compare_rad.pdf',
                     xtitle='P_{e} (GeV/c)', title='P_{e}')

    plot_sector_page(can, histos['data'], histos['sim'],
                     config1='Data', config2='Sim', title_formatter='histos_p_ele_FTOF_{}',
                     save_name='histos_p_electron_ftof_compare_rad.pdf',
                     xtitle='P_{e} (GeV/c)', title='P_{e}')
