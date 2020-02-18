import matplotlib.pyplot as plt
import numpy as np

from ROOT import TFile, TCanvas, TH1F

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

if __name__ == "__main__":

    input_filename = '../groovy/out.hipo.root'
    root_file = TFile(input_filename)
    histos = load_histos(root_file)

    data = {}
    data['CTOF'] = {}
    data['FTOF'] = {}
    for dtype in data:
        data[dtype]['n'] = []
        data[dtype]['run'] = []
        data[dtype]['mean'] = []
        data[dtype]['std'] = []
        data[dtype]['max'] = []


    can = TCanvas("can", "can", 800, 600)
    can.Print('elastic.pdf[')
    for name, hist in histos.items():
        if '_w_' in name:
            detector, run = name.split('_')[-2:]        
            data[detector]['run'].append(run)
            data[detector]['n'].append(sum([hist.GetBinContent(i) for i in range(1,hist.GetNbinsX())]))
            data[detector]['mean'].append(hist.GetMean())
            data[detector]['std'].append(hist.GetStdDev())
            data[detector]['max'].append(hist.GetBinCenter(hist.GetMaximumBin()))

            can.Clear()
            hist.Draw()
            can.Print('elastic.pdf')
    
    can.Print('elastic.pdf]')
    

    # Matplot 
    ylims = {}
    #ylims['n'] = [0,6000]
    ylims['mean'] = [0.85,1.1]
    ylims['max'] = [0.85,1.1]
    ylims['std'] = [0.0, 0.3]

    for dtype in data:
        for var in data[dtype]:
            if var is not 'run':
                plt.figure(figsize=(12,5))
                plt.plot(data[dtype]['run'], data[dtype][var], 
                         marker='o', linestyle='')
                plt.xlabel('Run Number')
                plt.ylabel(var)
                plt.title(dtype + ":" + var)
                plt.grid(alpha=0.2)
                
                if var in ylims:
                    plt.ylim(ylims[var])

                plt.savefig('timeline_' + dtype + '_' + var + '.pdf', 
                            bbox_inches='tight')
                

    
    
