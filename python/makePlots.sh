#!/bin/bash

# takes data/mc both
#python data_sim.py -d=rga.root -s=esepp.root -o=compare-in
#python data_sim_rad.py -d=rga-rad.root -s=esepp-rad.root -o=compare-in
python es.py -d=es-rga.root -s=es-esepp.root -o=es

# take only data or mc individually
#python mon.py -i=rga.root -o=rga
#python mon.py -i=esepp.root -o=esepp

#python monrad.py -i=rga-rad.root -o=rga-rad
#python monrad.py -i=esepp-rad.root -o=esepp-rad

