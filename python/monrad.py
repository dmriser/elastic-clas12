#!/usr/bin/env python 

import argparse 
import numpy as np

# Trick for docker install to run headless
# and never use plt.show() 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

from array import array 
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors,
                  TLine)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def setup_global_options():
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)
    
def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False,
                     y_fit_range=None, landscape=False, x_range=None, vline=None,
                     hline=None):

    root_garbage_can = []
    
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    for i in range(1,7):
        canvas.cd(i)
        
        #if isinstance(histos[title_formatter.format(i)], TH1F):
        if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
 
            if y_fit_range:
                fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])

            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                
            histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
            histos.get(title_formatter.format(i), default_histo).Draw()

            if y_fit_range:
                label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
        #elif isinstance(histos[title_formatter.format(i)], TH2F):
        elif isinstance(histos.get(title_formatter.format(i), default_histo), TH2F):
            # histos[title_formatter.format(i)].Draw('colz')
            histos.get(title_formatter.format(i), default_histo).Draw('colz')
            if log:
                gPad.SetLogz() 
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass

        if vline:
            line = TLine(vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmin(),
                         vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmax())
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)

        if hline:
            xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 
            line = TLine(xmin, hline, xmax, hline)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

def plot_page(canvas, histos, histo_title, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False):
    
    canvas.Clear() 
        
    if isinstance(histos[histo_title], TH1F):
        histos[histo_title].Draw()
    elif isinstance(histos[histo_title], TH2F):
        histos[histo_title].Draw('colz')
        if log:
            gPad.SetLogz() 
    else:
        #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
        pass
    
    if title:
        label.DrawLatex(0.1, 0.925, title)

    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    canvas.Print(save_name)


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

    return graph, slices, fits
    # return np.array(x_values), np.array(means), np.array(stds), slices, fits 

 
def plot_fits(canvas, histos, x_range, x_bin_step, title_formatter,
              save_name, label, y_fit_range, y_range=None,
              title=None, xtitle=None, ytitle=None, hline=None):

    canvas.Clear()
    canvas.Divide(2,3)

    root_is_dumb = []
    for i in range(1,7):
        canvas.cd(i)

        title = title_formatter.format(i)
        graph, slices, fits = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range=y_fit_range)
        graph.SetMarkerStyle(8)
        graph.SetMarkerSize(1)
        
        if y_range:
            graph.GetHistogram().SetMinimum(y_range[0])
            graph.GetHistogram().SetMaximum(y_range[1])
            graph.Draw('AP')
            root_is_dumb.append(graph)
        else:
            graph.Draw('AP')
            root_is_dumb.append(graph)
            
        if hline:
            line = TLine(x_range[0], hline, x_range[1], hline)
            line.SetLineStyle(8)
            line.SetLineWidth(1)
            line.Draw()
            root_is_dumb.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

    # For the slices 
    slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
    slice_pdfname = title_formatter.split('_{}')[0] + '_slices.pdf'
    slice_can.Print(slice_pdfname + '[')
    for i in range(1,7):
        title = title_formatter.format(i)
        graph, slices, fits = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range)

        # Size of slices page
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

    
def plot_fits_mpl(histos, x_range, x_bin_step, title_formatter,
                  save_name, y_range=None, title=None,
                  xtitle=None, ytitle=None, max_errorbar=0.5):
    
    fig = plt.figure(figsize=(12,16))

    opts = {'marker':'o', 'color':'k', 'linestyle':''}
    
    for i in range(1,7):
        ax = fig.add_subplot(3, 2, i)

        tit = title_formatter.format(i)
        x, mu, sig, slices, fits = fit_slices(histos.get(tit, default_histo2d), x_range, x_bin_step)

        # Here we remove the points with huge resolution, or
        # no events to fit.  Just leave the slices alone. 
        condition = np.logical_and(mu != 0, sig < max_errorbar)
        indices = np.where(condition)[0]
        x = x[indices]
        mu = mu[indices]
        sig = sig[indices]
        
        label = 'Sector {}'.format(i)
        if y_range:
            ax.errorbar(x, mu, sig, label=label, **opts)
            ax.set_ylim(y_range)
        else:
            ax.errorbar(x, mu, sig, label=label, **opts)

        # Add center line
        ax.axhline(0.0, linestyle='--', linewidth=1, color='k', alpha=0.85)

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

def add_text_page(can, label, text, save_name):
    """ Write some text onto a page. """
    can.Clear()
    can.cd(1)
    label.DrawLatex(0.1, 0.5, text)
    can.Print(save_name)
    
    
if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-i',
        '--input_file',
        required=True
    )
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    input_rootfile = args.input_file 
    output_pdfname = args.output_prefix + '.pdf'
    rootfile = TFile(input_rootfile)
    histos = load_histos(rootfile)

    for k,v in histos.items():
        print(k, v)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 800, 1100)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    # add_text_page(can, lab, text='W Monitoring Plots', save_name=output_pdfname)

    plot_sector_page(can, histos, 'histos_w_CTOF_{}', lab, save_name=output_pdfname,
                     title='(CTOF)', xtitle='W')
    
    plot_sector_page(can, histos, 'histos_w_FTOF_{}', lab, save_name=output_pdfname,
                     title='(FTOF)', xtitle='W')

    plot_sector_page(can, histos, 'histos_w_pass_angle_CTOF_{}', lab, save_name=output_pdfname,
                     title='Pass #Delta #phi Cut (CTOF)', xtitle='W')
    
    plot_sector_page(can, histos, 'histos_w_pass_angle_FTOF_{}', lab, save_name=output_pdfname,
                     title='Pass #Delta #phi Cut (FTOF)', xtitle='W')

    plot_sector_page(can, histos, 'histos_w_pass_all_CTOF_{}', lab, save_name=output_pdfname,
                     title='Pass All (CTOF)', xtitle='W')
    
    plot_sector_page(can, histos, 'histos_w_pass_all_FTOF_{}', lab, save_name=output_pdfname,
                     title='Pass All (FTOF)', xtitle='W')

    #    add_text_page(can, lab, text='Event Selection', save_name=output_pdfname)
    
    plot_sector_page(can, histos, 'histos_angle_ep_CTOF_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward), Proton (CTOF)', xtitle='#phi_{ep}')
    
    plot_sector_page(can, histos, 'histos_angle_ep_FTOF_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward), Proton (FTOF)', xtitle='#phi_{ep}')

    plot_sector_page(can, histos, 'histos_angle_ep_pass_theta_gamma_CTOF_{}', lab, save_name=output_pdfname,
                     title='(CTOF)  Pass #theta_{gamma}', xtitle='#phi_{ep}')
    
    plot_sector_page(can, histos, 'histos_angle_ep_pass_theta_gamma_FTOF_{}', lab, save_name=output_pdfname,
                     title='(FTOF) Pass #theta_{gamma}', xtitle='#phi_{ep}')

    plot_sector_page(can, histos, 'histos_theta_gamma_CTOF_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward), Proton (CTOF)', xtitle='#theta_{#gamma}')
    
    plot_sector_page(can, histos, 'histos_theta_gamma_FTOF_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward), Proton (FTOF)', xtitle='#theta_{#gamma}')
    
    plot_sector_page(can, histos, 'histos_theta_gamma_pass_angle_CTOF_{}', lab, save_name=output_pdfname,
                     title='Pass #Delta #phi (CTOF)', xtitle='#theta_{#gamma}')
     
    plot_sector_page(can, histos, 'histos_theta_gamma_pass_angle_FTOF_{}', lab, save_name=output_pdfname,
                     title='Pass #Delta #phi (FTOF)', xtitle='#theta_{#gamma}')

    plot_sector_page(can, histos, 'histos_p_ele_theta_ele_CTOF_{}', lab, save_name=output_pdfname,
                     title='(CTOF)', xtitle='p_{e}', ytitle='#theta_{e}', log=True)

    plot_sector_page(can, histos, 'histos_p_ele_theta_ele_FTOF_{}', lab, save_name=output_pdfname,
                     title='(FTOF)', xtitle='p_{e}', ytitle='#theta_{e}', log=True)

    plot_sector_page(can, histos, 'histos_p_pro_theta_pro_CTOF_{}', lab, save_name=output_pdfname,
                     title='(CTOF)', xtitle='p_{p}', ytitle='#theta_{p}', log=True)
 
    plot_sector_page(can, histos, 'histos_p_pro_theta_pro_FTOF_{}', lab, save_name=output_pdfname,
                     title='(FTOF)', xtitle='p_{p}', ytitle='#theta_{p}', log=True)
 
    plot_sector_page(can, histos, 'histos_theta_e_theta_gamma_pass_angle_CTOF_{}', lab, save_name=output_pdfname,
                     title='(CTOF)', xtitle='#theta_{e}', ytitle='#theta_{#gamma}', log=True)

    plot_sector_page(can, histos, 'histos_theta_e_theta_gamma_pass_angle_FTOF_{}', lab, save_name=output_pdfname,
                     title='(FTOF)', xtitle='#theta_{e}', ytitle='#theta_{#gamma}', log=True)

    add_text_page(can, lab, text='Resolutions', save_name=output_pdfname)

    plot_sector_page(can, histos, 'histos_p_ele_dp_ele_CTOF_{}', lab, save_name=output_pdfname,
                     title='Momentum Resolution (CTOF)', xtitle='p_{e}', ytitle='#Delta p_{e}',
                     log=True)

    plot_fits(can, histos, x_range=[1.5, 9.5], x_bin_step=6, title_formatter='histos_p_ele_dp_ele_CTOF_{}',
              save_name=output_pdfname, label=lab, y_range=[-4, 4],
              title='#Delta p_{e} vs. p_{e} (CTOF)', xtitle='p_{e}', ytitle='#Delta p_{e}', hline=0.0001, y_fit_range=[-4, 4])

    plot_sector_page(can, histos, 'histos_p_ele_dp_ele_FTOF_{}', lab, save_name=output_pdfname,
                     title='Momentum Resolution (FTOF)', xtitle='p_{e}', ytitle='#Delta p_{e}',
                     log=True)

    plot_fits(can, histos, x_range=[1.5, 9.5], x_bin_step=6, title_formatter='histos_p_ele_dp_ele_FTOF_{}',
              save_name=output_pdfname, label=lab, y_range=[-4, 4],
              title='#Delta p_{e} vs. p_{e} (FTOF)', xtitle='p_{e}', ytitle='#Delta p_{e}', hline=0.0001, y_fit_range=[-4, 4])

    plot_sector_page(can, histos, 'histos_p_pro_dp_pro_CTOF_{}', lab, save_name=output_pdfname,
                     title='Momentum Resolution (CTOF)', xtitle='p_{p}', ytitle='#Delta p_{p}',
                     log=True)

    plot_fits(can, histos, x_range=[0.0, 5.5], x_bin_step=6, title_formatter='histos_p_pro_dp_pro_CTOF_{}',
              save_name=output_pdfname, label=lab, y_range=[-4, 4],
              title='#Delta p_{p} vs. p_{p} (CTOF)', xtitle='p_{p}', ytitle='#Delta p_{p}', hline=0.0001, y_fit_range=[-4, 4])
    plot_sector_page(can, histos, 'histos_p_pro_dp_pro_FTOF_{}', lab, save_name=output_pdfname,
                     title='Momentum Resolution (FTOF)', xtitle='p_{p}', ytitle='#Delta p_{p}',
                     log=True)

    plot_fits(can, histos, x_range=[0.0, 5.5], x_bin_step=6, title_formatter='histos_p_pro_dp_pro_FTOF_{}',
              save_name=output_pdfname, label=lab, y_range=[-4, 4],
              title='#Delta p_{p} vs. p_{p} (FTOF)', xtitle='p_{p}', ytitle='#Delta p_{p}', hline=0.0001, y_fit_range=[-4, 4])
    
    plot_sector_page(can, histos, 'histos_theta_ele_dtheta_ele_CTOF_{}', lab, save_name=output_pdfname,
                     title='#Theta Resolution (CTOF)', xtitle='#theta_{e}', ytitle='#Delta #theta_{e}',
                     log=True)

    plot_sector_page(can, histos, 'histos_theta_ele_dtheta_ele_FTOF_{}', lab, save_name=output_pdfname,
                     title='#Theta Resolution (FTOF)', xtitle='#theta_{e}', ytitle='#Delta #theta_{e}',
                     log=True)

    plot_sector_page(can, histos, 'histos_theta_pro_dtheta_pro_CTOF_{}', lab, save_name=output_pdfname,
                     title='#Theta Resolution (CTOF)', xtitle='#theta_{p}', ytitle='#Delta #theta_{p}',
                     log=True)

    plot_sector_page(can, histos, 'histos_theta_pro_dtheta_pro_FTOF_{}', lab, save_name=output_pdfname,
                     title='#Theta Resolution (FTOF)', xtitle='#theta_{p}', ytitle='#Delta #theta_{p}',
                     log=True)
    
    #plot_fits(can, histos, x_range=[6.0,12.0], x_bin_step=6, title_formatter='histos_theta_w_ele_{}',
    #          save_name=output_pdfname, label=lab, y_range=[0.6, 1.5],
    #          title='W vs. #theta_{e}', xtitle='#theta_{e}', ytitle='W', hline=0.938, y_fit_range=[0.85, 1.1])

    can.Print('{}]'.format(output_pdfname))


    # single plots on landscape canvas                                                                                                                                    
    lcan = TCanvas('can', 'can', 1100, 800)
    """
    plot_sector_page(lcan, histos, 'histos_theta_ele_dtheta_ele_CTOF_{}', lab,
                     save_name='theta_ele_dtheta_ele_CTOF_{}.pdf'.format(args.output_prefix), 
                     title='#Delta #theta_{e} vs #theta_{e}', ytitle='#Delta #theta_{e}',
                     xtitle='#theta_{e}', log=True, landscape=True)

    plot_sector_page(lcan, histos, 'histos_theta_pro_dtheta_pro_CTOF_{}', lab,
                     save_name='theta_pro_dtheta_pro_CTOF_{}.pdf'.format(args.output_prefix), 
                     title='#Delta #theta_{p} vs #theta_{p}', ytitle='#Delta #theta_{p}',
                     xtitle='#theta_{p}', log=True, landscape=True)
    """
    plot_sector_page(lcan, histos, 'histos_p_ele_dp_ele_from_angles_CTOF_{}', lab,
                     save_name='p_ele_dp_ele_CTOF_{}.pdf'.format(args.output_prefix), 
                     title='#Delta p_{e} vs p_{e}', ytitle='#Delta p_{e}',
                     xtitle='p_{e}', log=True, landscape=True)

    plot_sector_page(lcan, histos, 'histos_p_pro_dp_pro_from_angles_CTOF_{}', lab,
                     save_name='p_pro_dp_pro_CTOF_{}.pdf'.format(args.output_prefix), 
                     title='#Delta p_{p} vs p_{p}', ytitle='#Delta p_{p}',
                     xtitle='p_{p}', log=True, landscape=True)

    plot_sector_page(lcan, histos, 'histos_p_ele_dp_ele_from_angles_CTOF_{}', lab,
                     save_name='p_ele_dp_ele_from_angles_CTOF_{}.pdf'.format(args.output_prefix), 
                     title='#Delta p_{e} vs p_{e}', ytitle='#Delta p_{e}',
                     xtitle='p_{e}', log=True, landscape=True)

    plot_sector_page(lcan, histos, 'histos_p_pro_dp_pro_from_angles_CTOF_{}', lab,
                     save_name='p_pro_dp_pro_from_angles_CTOF_{}.pdf'.format(args.output_prefix), 
                     title='#Delta p_{p} vs p_{p}', ytitle='#Delta p_{p}',
                     xtitle='p_{p}', log=True, landscape=True)
