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
                  gPad, gStyle, TLatex, TGraphErrors, TLine,
                  kGray, kRed, kBlue, kBlack)

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

def get_min_max(histos):
#    return min([h.GetMinimum() for h in histos]), max([h.GetMaximum() for h in histos])
    return histos[0].GetMinimum(), max([h.GetMaximum() for h in histos])

def color_draw(histo, color=kGray, opts=""):
    histo.SetLineColor(kBlack)
    histo.SetFillColorAlpha(color,1.0)
    histo.Draw(opts)
    
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

def plot_phase_space(can, histos, keyword, output_prefix):
    """ Docstring. """

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_p_ele_theta_ele_'+keyword+'_CTOF'].Draw('colz')
    latex.DrawLatex(0.3, 0.92, 'Data (proton in CTOF)')
    
    can.cd(2)
    histos['data']['histos_p_ele_theta_ele_'+keyword+'_FTOF'].Draw('colz')
    latex.DrawLatex(0.3, 0.92, 'Data (proton in FTOF)')
    
    can.cd(3)
    histos['sim']['histos_p_ele_theta_ele_'+keyword+'_CTOF'].Draw('colz')
    latex.DrawLatex(0.3, 0.92, 'Sim (proton in CTOF)')
    
    can.cd(4)
    histos['sim']['histos_p_ele_theta_ele_'+keyword+'_FTOF'].Draw('colz')
    latex.DrawLatex(0.3, 0.92, 'Sim (proton in FTOF)')

    for ican in range(1, 5):
        can.cd(ican)
        latex.DrawLatex(0.45, 0.02, 'p_{e} (GeV)')
        latex.SetTextAngle(90.0)
        latex.DrawLatex(0.03, 0.45, '#theta_{e} (Deg.)')
        latex.SetTextAngle(0.0)
            
    can.Print('phase_space_' + keyword + '_' + output_prefix + '.pdf')
    can.Clear()
    
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

    # Setup Cuts
    cuts = {}
    cuts['w'] = [1.25, 6.0]
    cuts['angle_ep'] = [178.0, 180.0]
    cuts['theta_gamma'] = [0.0, 3.0]
    cuts['missing_mass'] = [-0.4, 0.4]
    
    # Global opts
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    root_garbage_can = [] 
    
    can = TCanvas('can', 'can', 1100, 800)

    plot_phase_space(can, histos,
                     keyword='isr',
                     output_prefix=args.output_prefix)
    plot_phase_space(can, histos,
                     keyword='elastic',
                     output_prefix=args.output_prefix)

    # -----------------------------------------------------------
    # Plot W
    # -----------------------------------------------------------
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    
    can.Divide(2,2)

    can.cd(1)
    color_draw(histos['data']['histos_w_inclusive_'], kGray, "")
    color_draw(histos['data']['histos_w_CTOF'], kBlue, "same")
    color_draw(histos['data']['histos_w_pass_angle_CTOF'], kRed, "same")
    latex.DrawLatex(0.45, 0.02, 'W (GeV/c^{2})')
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    latex.SetTextColor(kBlue)
    latex.DrawLatex(0.15, 0.85, 'w/ proton')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.80, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    
    can.cd(2)
    color_draw(histos['data']['histos_w_inclusive_'], kGray, "")
    color_draw(histos['data']['histos_w_FTOF'], kBlue, "same")
    color_draw(histos['data']['histos_w_pass_angle_FTOF'], kRed, "same")
    latex.DrawLatex(0.45, 0.02, 'W (GeV/c^{2})')
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF')
    latex.SetTextColor(kBlue)
    latex.DrawLatex(0.15, 0.85, 'w/ proton')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.80, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    
    can.cd(3)
    color_draw(histos['sim']['histos_w_inclusive_'], kGray, "")
    color_draw(histos['sim']['histos_w_CTOF'], kBlue, "same")
    color_draw(histos['sim']['histos_w_pass_angle_CTOF'], kRed, "same")
    latex.DrawLatex(0.45, 0.02, 'W (GeV/c^{2})')
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    latex.SetTextColor(kBlue)
    latex.DrawLatex(0.75, 0.85, 'w/ proton')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.75, 0.80, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    
    can.cd(4)
    color_draw(histos['sim']['histos_w_inclusive_'], kGray, "")
    color_draw(histos['sim']['histos_w_FTOF'], kBlue, "same")
    color_draw(histos['sim']['histos_w_pass_angle_FTOF'], kRed, "same")
    latex.DrawLatex(0.45, 0.02, 'W (GeV/c^{2})')
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF')
    latex.SetTextColor(kBlue)
    latex.DrawLatex(0.75, 0.85, 'w/ proton')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.75, 0.80, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)

    can.Print('w_' + args.output_prefix + '.pdf')

    # -----------------------------------------------------------
    # Plot delta phi
    # -----------------------------------------------------------

    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_theta_gamma_CTOF'].SetLineColor(kBlack)
    histos['data']['histos_theta_gamma_CTOF'].SetFillColorAlpha(kGray,1.0)
    histos['data']['histos_theta_gamma_CTOF'].Draw()
    histos['data']['histos_theta_gamma_pass_angle_CTOF'].SetLineColor(kBlack)
    histos['data']['histos_theta_gamma_pass_angle_CTOF'].SetFillColorAlpha(kRed,1.0)
    histos['data']['histos_theta_gamma_pass_angle_CTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#theta_{#gamma} (deg)')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')

    ymin, ymax = get_min_max([
        histos['data']['histos_theta_gamma_CTOF'],
        histos['data']['histos_theta_gamma_pass_angle_CTOF']]
    )
    line = TLine(3.0, ymin, 3.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(2)
    histos['data']['histos_theta_gamma_FTOF'].SetLineColor(kBlack)
    histos['data']['histos_theta_gamma_FTOF'].SetFillColorAlpha(kGray,1.0)
    histos['data']['histos_theta_gamma_FTOF'].Draw()
    histos['data']['histos_theta_gamma_pass_angle_FTOF'].SetLineColor(kBlack)
    histos['data']['histos_theta_gamma_pass_angle_FTOF'].SetFillColorAlpha(kRed,1.0)
    histos['data']['histos_theta_gamma_pass_angle_FTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#theta_{#gamma} (deg)')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.72, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF')

    ymin, ymax = get_min_max([
        histos['data']['histos_theta_gamma_FTOF'],
        histos['data']['histos_theta_gamma_pass_angle_FTOF']]
    )
    line = TLine(3.0, ymin, 3.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
        
    can.cd(3)
    histos['sim']['histos_theta_gamma_CTOF'].SetLineColor(kBlack)
    histos['sim']['histos_theta_gamma_CTOF'].SetFillColorAlpha(kGray,1.0)
    histos['sim']['histos_theta_gamma_CTOF'].Draw()
    histos['sim']['histos_theta_gamma_pass_angle_CTOF'].SetLineColor(kBlack)
    histos['sim']['histos_theta_gamma_pass_angle_CTOF'].SetFillColorAlpha(kRed,1.0)
    histos['sim']['histos_theta_gamma_pass_angle_CTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#theta_{#gamma} (deg)')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    ymin, ymax = get_min_max([
        histos['sim']['histos_theta_gamma_CTOF'],
        histos['sim']['histos_theta_gamma_pass_angle_CTOF']]
    )
    line = TLine(3.0, ymin, 3.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
    
    can.cd(4)
    histos['sim']['histos_theta_gamma_FTOF'].SetLineColor(kBlack)
    histos['sim']['histos_theta_gamma_FTOF'].SetFillColorAlpha(kGray,1.0)
    histos['sim']['histos_theta_gamma_FTOF'].Draw()
    histos['sim']['histos_theta_gamma_pass_angle_FTOF'].SetLineColor(kBlack)
    histos['sim']['histos_theta_gamma_pass_angle_FTOF'].SetFillColorAlpha(kRed,1.0)
    histos['sim']['histos_theta_gamma_pass_angle_FTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#theta_{#gamma} (deg)')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.72, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF')
    ymin, ymax = get_min_max([
        histos['sim']['histos_theta_gamma_FTOF'],
        histos['sim']['histos_theta_gamma_pass_angle_FTOF']]
    )
    line = TLine(3.0, ymin, 3.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
    
    can.Print('theta_gamma_' + args.output_prefix + '.pdf')
     
    # -----------------------------------------------------------
    # Plot angle EP 
    # -----------------------------------------------------------

    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_angle_ep_CTOF'].SetLineColor(kBlack)
    histos['data']['histos_angle_ep_CTOF'].SetFillColorAlpha(kGray,1.0)
    histos['data']['histos_angle_ep_CTOF'].Draw()
    histos['data']['histos_angle_ep_pass_w_elastic_CTOF'].SetLineColor(kBlack)
    histos['data']['histos_angle_ep_pass_w_elastic_CTOF'].SetFillColorAlpha(kRed,1.0)
    histos['data']['histos_angle_ep_pass_w_elastic_CTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#Delta#phi (deg)')
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.68, 0.85, 'Elastic Peak')
    latex.SetTextColor(kBlack)

    ymin, ymax = get_min_max([
        histos['data']['histos_angle_ep_CTOF'],
        histos['data']['histos_angle_ep_pass_w_elastic_CTOF']]
    )
    line = TLine(178.0, ymin, 178.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
    
    can.cd(2)
    histos['data']['histos_angle_ep_FTOF'].SetLineColor(kBlack)
    histos['data']['histos_angle_ep_FTOF'].SetFillColorAlpha(kGray,1.0)
    histos['data']['histos_angle_ep_FTOF'].Draw()
    histos['data']['histos_angle_ep_pass_w_elastic_FTOF'].SetLineColor(kBlack)
    histos['data']['histos_angle_ep_pass_w_elastic_FTOF'].SetFillColorAlpha(kRed,1.0)
    histos['data']['histos_angle_ep_pass_w_elastic_FTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#Delta#phi (deg)')
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.68, 0.85, 'Elastic Peak')
    latex.SetTextColor(kBlack)

    ymin, ymax = get_min_max([
        histos['data']['histos_angle_ep_FTOF'],
        histos['data']['histos_angle_ep_pass_w_elastic_FTOF']]
    )
    line = TLine(178.0, ymin, 178.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
    
    can.cd(3)
    histos['sim']['histos_angle_ep_CTOF'].SetLineColor(kBlack)
    histos['sim']['histos_angle_ep_CTOF'].SetFillColorAlpha(kGray,1.0)
    histos['sim']['histos_angle_ep_CTOF'].Draw()
    histos['sim']['histos_angle_ep_pass_w_elastic_CTOF'].SetLineColor(kBlack)
    histos['sim']['histos_angle_ep_pass_w_elastic_CTOF'].SetFillColorAlpha(kRed,1.0)
    histos['sim']['histos_angle_ep_pass_w_elastic_CTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#Delta#phi (deg)')
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.68, 0.85, 'Elastic Peak')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['sim']['histos_angle_ep_CTOF'],
        histos['sim']['histos_angle_ep_pass_w_elastic_CTOF']]
    )
    line = TLine(178.0, ymin, 178.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(4)
    histos['sim']['histos_angle_ep_FTOF'].SetLineColor(kBlack)
    histos['sim']['histos_angle_ep_FTOF'].SetFillColorAlpha(kGray,1.0)
    histos['sim']['histos_angle_ep_FTOF'].Draw()
    histos['sim']['histos_angle_ep_pass_w_elastic_FTOF'].SetLineColor(kBlack)
    histos['sim']['histos_angle_ep_pass_w_elastic_FTOF'].SetFillColorAlpha(kRed,1.0)
    histos['sim']['histos_angle_ep_pass_w_elastic_FTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, '#Delta#phi (deg)')
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.68, 0.85, 'Elastic Peak')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['sim']['histos_angle_ep_FTOF'],
        histos['sim']['histos_angle_ep_pass_w_elastic_FTOF']]
    )
    line = TLine(178.0, ymin, 178.0, ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.Print('angle_ep_' + args.output_prefix + '.pdf')

        
    # -----------------------------------------------------------
    # Plot all the plots 
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_theta_e_theta_gamma_pass_angle_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, '#theta_{e} (deg)')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{#gamma} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    
    can.cd(2)
    histos['data']['histos_theta_e_theta_gamma_pass_angle_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, '#theta_{e} (deg)')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{#gamma} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF')
    
    can.cd(3)
    histos['sim']['histos_theta_e_theta_gamma_pass_angle_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, '#theta_{e} (deg)')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{#gamma} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    
    can.cd(4)
    histos['sim']['histos_theta_e_theta_gamma_pass_angle_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, '#theta_{e} (deg)')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{#gamma} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF')
        
    can.Print('theta_theta_' + args.output_prefix + '.pdf')


    # -----------------------------------------------------------
    # Plot Angular Sum
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_w_theta_sum_pass_angle_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')

    can.cd(2)
    histos['data']['histos_w_theta_sum_pass_angle_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF')

    can.cd(3)
    histos['sim']['histos_w_theta_sum_pass_angle_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')

    can.cd(4)
    histos['sim']['histos_w_theta_sum_pass_angle_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF')

    can.Print('w_theta_sum_' + args.output_prefix + '.pdf')

    # -----------------------------------------------------------
    # Plot Missing Mass
    # -----------------------------------------------------------

    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_missing_mass_CTOF'].SetLineColor(kBlack)
    histos['data']['histos_missing_mass_CTOF'].SetFillColorAlpha(kGray,1.0)
    histos['data']['histos_missing_mass_CTOF'].Draw()
    histos['data']['histos_missing_mass_pass_angle_CTOF'].SetLineColor(kBlack)
    histos['data']['histos_missing_mass_pass_angle_CTOF'].SetFillColorAlpha(kRed,1.0)
    histos['data']['histos_missing_mass_pass_angle_CTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, 'M_{X}^{2}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    ymin, ymax = get_min_max([
        histos['data']['histos_missing_mass_CTOF'],
        histos['data']['histos_missing_mass_pass_angle_CTOF']]
    )
    lline = TLine(-0.08, ymin, -0.08, ymax)
    lline.SetLineColor(1)
    lline.SetLineStyle(1)
    lline.Draw('same')
    root_garbage_can.append(lline)
    rline = TLine(0.08, ymin, 0.08, ymax)
    rline.SetLineColor(1)
    rline.SetLineStyle(1)
    rline.Draw('same')
    root_garbage_can.append(rline)

    can.cd(2)
    histos['data']['histos_missing_mass_FTOF'].SetLineColor(kBlack)
    histos['data']['histos_missing_mass_FTOF'].SetFillColorAlpha(kGray,1.0)
    histos['data']['histos_missing_mass_FTOF'].Draw()
    histos['data']['histos_missing_mass_pass_angle_FTOF'].SetLineColor(kBlack)
    histos['data']['histos_missing_mass_pass_angle_FTOF'].SetFillColorAlpha(kRed,1.0)
    histos['data']['histos_missing_mass_pass_angle_FTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, 'M_{X}^{2}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF')
    ymin, ymax = get_min_max([
        histos['data']['histos_missing_mass_FTOF'],
        histos['data']['histos_missing_mass_pass_angle_FTOF']]
    )
    lline = TLine(-0.08, ymin, -0.08, ymax)
    lline.SetLineColor(1)
    lline.SetLineStyle(1)
    lline.Draw('same')
    root_garbage_can.append(lline)
    rline = TLine(0.08, ymin, 0.08, ymax)
    rline.SetLineColor(1)
    rline.SetLineStyle(1)
    rline.Draw('same')
    root_garbage_can.append(rline)

    can.cd(3)
    histos['sim']['histos_missing_mass_CTOF'].SetLineColor(kBlack)
    histos['sim']['histos_missing_mass_CTOF'].SetFillColorAlpha(kGray,1.0)
    histos['sim']['histos_missing_mass_CTOF'].Draw()
    histos['sim']['histos_missing_mass_pass_angle_CTOF'].SetLineColor(kBlack)
    histos['sim']['histos_missing_mass_pass_angle_CTOF'].SetFillColorAlpha(kRed,1.0)
    histos['sim']['histos_missing_mass_pass_angle_CTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, 'M_{X}^{2}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    ymin, ymax = get_min_max([
        histos['sim']['histos_missing_mass_CTOF'],
        histos['sim']['histos_missing_mass_pass_angle_CTOF']]
    )
    lline = TLine(-0.08, ymin, -0.08, ymax)
    lline.SetLineColor(1)
    lline.SetLineStyle(1)
    lline.Draw('same')
    root_garbage_can.append(lline)
    rline = TLine(0.08, ymin, 0.08, ymax)
    rline.SetLineColor(1)
    rline.SetLineStyle(1)
    rline.Draw('same')
    root_garbage_can.append(rline)

    
    can.cd(4)
    histos['sim']['histos_missing_mass_FTOF'].SetLineColor(kBlack)
    histos['sim']['histos_missing_mass_FTOF'].SetFillColorAlpha(kGray,1.0)
    histos['sim']['histos_missing_mass_FTOF'].Draw()
    histos['sim']['histos_missing_mass_pass_angle_FTOF'].SetLineColor(kBlack)
    histos['sim']['histos_missing_mass_pass_angle_FTOF'].SetFillColorAlpha(kRed,1.0)
    histos['sim']['histos_missing_mass_pass_angle_FTOF'].Draw('same')
    latex.DrawLatex(0.45, 0.02, 'M_{X}^{2}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'w/ #Delta #phi cut')
    latex.SetTextColor(kBlack)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF')
    ymin, ymax = get_min_max([
        histos['sim']['histos_missing_mass_FTOF'],
        histos['sim']['histos_missing_mass_pass_angle_FTOF']]
    )
    lline = TLine(-0.4, ymin, -0.4, ymax)
    lline.SetLineColor(1)
    lline.SetLineStyle(1)
    lline.Draw('same')
    root_garbage_can.append(lline)
    rline = TLine(0.4, ymin, 0.4, ymax)
    rline.SetLineColor(1)
    rline.SetLineStyle(1)
    rline.Draw('same')
    root_garbage_can.append(rline)

    
    can.Print('missing_mass_' + args.output_prefix + '.pdf')
    
    # -----------------------------------------------------------
    # Plot Angular Sum
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_w_theta_sum_pass_missing_mass_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF (Pass M_{X})')

    can.cd(2)
    histos['data']['histos_w_theta_sum_pass_missing_mass_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF (Pass M_{X})')

    can.cd(3)
    histos['sim']['histos_w_theta_sum_pass_missing_mass_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF (Pass M_{X})')

    can.cd(4)
    histos['sim']['histos_w_theta_sum_pass_missing_mass_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF (Pass M_{X})')

    can.Print('w_theta_sum_pass_missing_mass_' + args.output_prefix + '.pdf')

    # -----------------------------------------------------------
    # Plot Angular Sum
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_w_theta_sum_pass_missing_mass_angles_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF (Pass M_{X}, #Delta #phi)')

    can.cd(2)
    histos['data']['histos_w_theta_sum_pass_missing_mass_angles_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF (Pass M_{X}, #Delta #phi)')

    can.cd(3)
    histos['sim']['histos_w_theta_sum_pass_missing_mass_angles_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF (Pass M_{X}, #Delta #phi)')

    can.cd(4)
    histos['sim']['histos_w_theta_sum_pass_missing_mass_angles_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF (Pass M_{X}, #Delta #phi)')

    can.Print('w_theta_sum_pass_missing_mass_angles_' + args.output_prefix + '.pdf')

    # -----------------------------------------------------------
    # Plot Angular Sum
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    histos['data']['histos_w_theta_sum_isr_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF (ISR)')

    can.cd(2)
    histos['data']['histos_w_theta_sum_isr_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in FTOF (ISR)')

    can.cd(3)
    histos['sim']['histos_w_theta_sum_isr_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF (ISR)')

    can.cd(4)
    histos['sim']['histos_w_theta_sum_isr_FTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, '#theta_{e} + #theta_{p} (deg)')
    latex.SetTextAngle(0.0)
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in FTOF (ISR)')

    can.Print('w_theta_sum_isr_' + args.output_prefix + '.pdf')

    # -----------------------------------------------------------
    # Plot Drift Chambers
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(3,2)

    can.cd(1)
    histos['data']['histos_dc1_elastic_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'X')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, 'Y')
    latex.SetTextAngle(0.0) 
    latex.DrawLatex(0.15, 0.95, 'Data w/ Proton in CTOF (Elastic) DC1')

    can.cd(2)
    histos['data']['histos_dc2_elastic_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'X')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, 'Y')
    latex.SetTextAngle(0.0) 
    latex.DrawLatex(0.15, 0.95, 'Data w/ Proton in CTOF (Elastic) DC2')

    can.cd(3)
    histos['data']['histos_dc3_elastic_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'X')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, 'Y')
    latex.SetTextAngle(0.0) 
    latex.DrawLatex(0.15, 0.95, 'Data w/ Proton in CTOF (Elastic) DC3')

    can.cd(4)
    histos['data']['histos_dc1_isr_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'X')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, 'Y')
    latex.SetTextAngle(0.0) 
    latex.DrawLatex(0.15, 0.95, 'Data w/ Proton in CTOF (ISR) DC1')

    can.cd(5)
    histos['data']['histos_dc2_isr_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'X')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, 'Y')
    latex.SetTextAngle(0.0) 
    latex.DrawLatex(0.15, 0.95, 'Data w/ Proton in CTOF (ISR) DC2')

    can.cd(6)
    histos['data']['histos_dc3_isr_CTOF'].Draw('colz')
    gPad.SetLogz()
    latex.DrawLatex(0.45, 0.02, 'X')
    latex.SetTextAngle(90.0)
    latex.DrawLatex(0.02, 0.4, 'Y')
    latex.SetTextAngle(0.0) 
    latex.DrawLatex(0.15, 0.95, 'Data w/ Proton in CTOF (ISR) DC3')

    can.Print('dc_summary_' + args.output_prefix + '.pdf')


    # -----------------------------------------------------------
    # Plot Summary of DATA/CTOF event selection
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    color_draw(histos['data']['histos_w_CTOF'], kGray, "")
    color_draw(histos['data']['histos_w_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['data']['histos_w_CTOF'],
        histos['data']['histos_w_eep_CTOF']]
    )
    line = TLine(cuts['w'][0], ymin, cuts['w'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(2)
    color_draw(histos['data']['histos_theta_gamma_CTOF'], kGray, "")
    color_draw(histos['data']['histos_theta_gamma_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, '#theta_{#gamma}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['data']['histos_theta_gamma_CTOF'],
        histos['data']['histos_theta_gamma_eep_CTOF']]
    )
    line = TLine(cuts['theta_gamma'][1], ymin, cuts['theta_gamma'][1], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(3)
    color_draw(histos['data']['histos_angle_ep_CTOF'], kGray, "")
    color_draw(histos['data']['histos_angle_ep_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, '#Delta#phi_{ep}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['data']['histos_angle_ep_CTOF'],
        histos['data']['histos_angle_ep_eep_CTOF']]
    )
    line = TLine(cuts['angle_ep'][0], ymin, cuts['angle_ep'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(4)
    color_draw(histos['data']['histos_missing_mass_CTOF'], kGray, "")
    color_draw(histos['data']['histos_missing_mass_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Data w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, 'M_{X}^{2}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['data']['histos_missing_mass_CTOF'],
        histos['data']['histos_missing_mass_eep_CTOF']]
    )
    line = TLine(cuts['missing_mass'][0], ymin, cuts['missing_mass'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
    line = TLine(cuts['missing_mass'][1], ymin, cuts['missing_mass'][1], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.Print('es_summary_data_ctof_' + args.output_prefix + '.pdf')

    # -----------------------------------------------------------
    # Plot Summary of SIM/CTOF event selection
    # -----------------------------------------------------------
    
    can.Clear()
    can.Divide(2,2)

    can.cd(1)
    color_draw(histos['sim']['histos_w_CTOF'], kGray, "")
    color_draw(histos['sim']['histos_w_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, 'W')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['sim']['histos_w_CTOF'],
        histos['sim']['histos_w_eep_CTOF']]
    )
    line = TLine(cuts['w'][0], ymin, cuts['w'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(2)
    color_draw(histos['sim']['histos_theta_gamma_CTOF'], kGray, "")
    color_draw(histos['sim']['histos_theta_gamma_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, '#theta_{#gamma}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['sim']['histos_theta_gamma_CTOF'],
        histos['sim']['histos_theta_gamma_eep_CTOF']]
    )
    line = TLine(cuts['theta_gamma'][0], ymin, cuts['theta_gamma'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(3)
    color_draw(histos['sim']['histos_angle_ep_CTOF'], kGray, "")
    color_draw(histos['sim']['histos_angle_ep_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, '#Delta#phi_{ep}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['sim']['histos_angle_ep_CTOF'],
        histos['sim']['histos_angle_ep_eep_CTOF']]
    )
    line = TLine(cuts['angle_ep'][0], ymin, cuts['angle_ep'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.cd(4)
    color_draw(histos['sim']['histos_missing_mass_CTOF'], kGray, "")
    color_draw(histos['sim']['histos_missing_mass_eep_CTOF'], kRed, "same")
    latex.DrawLatex(0.3, 0.95, 'Sim w/ Proton in CTOF')
    latex.DrawLatex(0.45, 0.02, 'M_{X}^{2}')
    latex.SetTextColor(kRed)
    latex.DrawLatex(0.15, 0.85, 'Pass Others')
    latex.SetTextColor(kBlack)
    ymin, ymax = get_min_max([
        histos['sim']['histos_missing_mass_CTOF'],
        histos['sim']['histos_missing_mass_eep_CTOF']]
    )
    line = TLine(cuts['missing_mass'][0], ymin, cuts['missing_mass'][0], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)
    line = TLine(cuts['missing_mass'][1], ymin, cuts['missing_mass'][1], ymax)
    line.SetLineColor(1)
    line.SetLineStyle(1)
    line.Draw('same')
    root_garbage_can.append(line)

    can.Print('es_summary_sim_ctof_' + args.output_prefix + '.pdf')
