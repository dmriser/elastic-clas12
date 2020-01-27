#!/usr/bin/groovy
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.math.F1D
import com.google.gson.Gson
import java.util.zip.ZipOutputStream
import java.util.zip.ZipEntry
import java.io.File
import java.io.FileOutputStream
import org.jlab.groot.data.Directory
import org.jlab.groot.data.TDirectory


def convertH1F = { h1 ->
  def data = []
  h1.getDataSize(0).times{
    if(h1.getDataY(it)!=0)
      data.add([h1.getDataX(it), h1.getDataY(it)])
  }
  return [ type: 'H1F',
    name: h1.getName(),
    title: h1.getTitle(),
    xtitle: h1.getTitleX(),
    ytitle: h1.getTitleY(),
    xbins: h1.getXaxis().getNBins(),
    xmin: h1.getXaxis().min(),
    xmax: h1.getXaxis().max(),
    data: data ]
}

def convertH2F = { h2 ->
  def data = []
  h2.getDataSize(0).times{ ix ->
    h2.getDataSize(1).times{ iy ->
      if(h2.getData(ix,iy)!=0)
        data.add([h2.getDataX(ix), h2.getDataY(iy), h2.getData(ix,iy)])
    }
  }
  return [ type: 'H2F',
    name: h2.getName(),
    title: h2.getTitle(),
    xtitle: h2.getTitleX(),
    ytitle: h2.getTitleY(),
    xbins: h2.getXAxis().getNBins(),
    xmin: h2.getXAxis().min(),
    xmax: h2.getXAxis().max(),
    ybins: h2.getYAxis().getNBins(),
    ymin: h2.getYAxis().min(),
    ymax: h2.getYAxis().max(),
    data: data ]
}


def convertH3F = { h3 ->
  def data = []
  def xa = h3.getXAxis()
  def ya = h3.getYAxis()
  def za = h3.getZAxis()
  xa.getNBins().times{ ix ->
    ya.getNBins().times{ iy ->
      za.getNBins().times{ iz ->
        if(h3.getBinContent(ix,iy,iz)!=0)
          data.add([xa.getBinCenter(ix), ya.getBinCenter(iy), za.getBinCenter(iz), h3.getBinContent(ix,iy,iz)])
      }
    }
  }
  return [ type: 'H3F',
    name: h3.getName(),
    title: h3.getTitle(),
    xtitle: h3.getTitleX(),
    ytitle: h3.getTitleY(),
    xbins: h3.getXAxis().getNBins(),
    xmin: h3.getXAxis().min(),
    xmax: h3.getXAxis().max(),
    ybins: h3.getYAxis().getNBins(),
    ymin: h3.getYAxis().min(),
    ymax: h3.getYAxis().max(),
    zbins: h3.getZAxis().getNBins(),
    zmin: h3.getZAxis().min(),
    zmax: h3.getZAxis().max(),
    data: data ]
}


def convertGraph = { gr ->
  def data = []
  gr.getDataSize(0).times{
    data.add([gr.getDataX(it), gr.getDataY(it), gr.getDataEX(it), gr.getDataEY(it)])
  }
  return [ type: 'Graph',
    name: gr.getName(),
    title: gr.getTitle(),
    xtitle: gr.getTitleX(),
    ytitle: gr.getTitleY(),
    data: data ]
}


def convertF1D = { f1 ->
  return [ type: 'F1',
    name: f1.getName(),
    expression: f1.getExpression(),
    xmin: f1.getRange().getMin(),
    xmax: f1.getRange().getMax() ]
}



def gson = new Gson()

def fout = File.createTempFile('j2root','.zip')
def zipFile = new ZipOutputStream(new FileOutputStream(fout))

def out = new TDirectory()
out.readFile(args[0])

def addDS(dir, path='', dsets=[:]) {
  dir.getDirectoryList().each{
    addDS(dir.getDir(it), [path,it].join('/'), dsets)
  }
  dsets << dir.getObjectMap().collectEntries{ [([path,it.key].join('/')): it.value] }
  return dsets
}

def dsout = [:]
def dsets = addDS(out)

println(dsets.values()*.name)
dsets.each{
  def className = it.value.getClass().getName()
  if(className.contains('H1F')) {
    dsout[it.key] = convertH1F(it.value)
  } else if(className.contains('H2F')) {
    dsout[it.key] = convertH2F(it.value)
  } else if(className.contains('F1D')) {
    dsout[it.key] = convertF1D(it.value)
  } else if(className.contains('Graph')) {
    dsout[it.key] = convertGraph(it.value)
  }
}

zipFile.putNextEntry(new ZipEntry('datasets.json'))
zipFile << gson.toJson(dsout)
zipFile.closeEntry()
zipFile.close()

def pythonCode = '''
import json, sys
#sys.argv.append('-b')
import ROOT
from zipfile import ZipFile

def convertGraph(name, gr):
  grout = ROOT.TGraphErrors()
  grout.SetName(name)
  grout.SetTitle(';'.join([gr[kk] for kk in ['title','xtitle','ytitle']]))
  for x,y,dx,dy in gr['data']:
    ip = grout.GetN()
    grout.SetPoint(ip,x,y)
    grout.SetPointError(ip,dx,dy)
  return grout

def convertH1F(name, h1):
  hh = ROOT.TH1F(name, *[h1[kk] for kk in ['title','xbins','xmin','xmax']])
  hh.SetTitle(';'.join([h1[kk] for kk in ['title','xtitle','ytitle']]))
  for x,y in h1['data']:
    ib = hh.Fill(x)
    hh.SetBinContent(ib,y)
  return hh

def convertH2F(name, h2):
  hh = ROOT.TH2F(name, *[h2[kk] for kk in ['title','xbins','xmin','xmax','ybins','ymin','ymax']])
  hh.SetTitle(';'.join([h2[kk] for kk in ['title','xtitle','ytitle']]))
  for x,y,z in h2['data']:
    ib = hh.Fill(x,y)
    hh.SetBinContent(ib,z)
  return hh

def convertH3F(name, h3):
  hh = ROOT.TH3F(name, *[h3[kk] for kk in ['title','xbins','xmin','xmax','ybins','ymin','ymax','zbins','zmin','zmax']])
  hh.SetTitle(';'.join([h3[kk] for kk in ['title','xtitle','ytitle']]))
  for x,y,z,u in h3['data']:
    ib = hh.Fill(x,y,z)
    hh.SetBinContent(ib,u)
  return hh


with ZipFile(sys.argv[1], 'r') as zip:
  with zip.open('datasets.json') as ff:
    data = json.load(ff)


ff = ROOT.TFile(sys.argv[2]+'.root', "recreate")
dsets = []
for name in data:
  rootname = name[1:].replace('/','_').replace(' ','_')
  if 'Graph' in data[name]['type']:
    gr = convertGraph(rootname,data[name])
    dsets.append(gr)
  elif 'H1F' in data[name]['type']:
    h1 = convertH1F(rootname,data[name])
    dsets.append(h1)
  elif 'H2F' in data[name]['type']:
    h2 = convertH2F(rootname,data[name])
    dsets.append(h2)
  elif 'H3F' in data[name]['type']:
    h3 = convertH3F(rootname,data[name])
    dsets.append(h3)

for ds in dsets:
  print(ds)
  ds.Write()

ff.Close()
'''

def pythonFile = File.createTempFile('j2root','.py')
pythonFile<<pythonCode

def pycmd = 'python '+pythonFile.getAbsolutePath()+' '+fout.getAbsolutePath()+' '+args[0]
println pycmd
def proc = pycmd.execute()
println proc.err.text
println proc.text

