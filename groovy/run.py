"""
Start the analysis off right.
"""

import argparse
import subprocess
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--event_selection', required=False, action='store_true')
    parser.add_argument('--monitor', required=False, action='store_true')
    parser.add_argument('--monitor_rad', required=False, action='store_true')
    args = parser.parse_args()
    
    run_cmd = 'run-groovy'
    es_out = 'event-selection.hipo'
    mon_out = 'monitor.hipo'
    monrad_out = 'mon-rad.hipo'
    es_exe = 'event-selection.groovy'
    mon_exe = 'monitor.groovy'
    monrad_exe = 'monitor-rad.groovy'
    convert_exe = 'j2root.groovy'

    mc_files = '/volatile/clas12/dmriser/esepp/cooked/inbend_pass3/*'

    with open('smallFiles.txt', 'r') as inp:
        lines = inp.readlines()
        data_files = "".join([d.replace('\n',' ') for d in lines])

    if args.event_selection:
        subprocess.call([run_cmd, es_exe, data_files])
        subprocess.call([run_cmd, convert_exe, es_out])
        subprocess.call(['mv', es_out + '.root', 'es-rga.root'])

        subprocess.call([run_cmd, es_exe, mc_files])
        subprocess.call([run_cmd, convert_exe, es_out])
        subprocess.call(['mv', es_out + '.root', 'es-esepp.root'])

    if args.monitor:
        #subprocess.call([run_cmd, mon_exe, data_files])
        #subprocess.call([run_cmd, convert_exe, mon_out])
        #subprocess.call(['mv', mon_out + '.root', 'rga.root'])

        subprocess.call([run_cmd, mon_exe, mc_files])
        subprocess.call([run_cmd, convert_exe, mon_out])
        subprocess.call(['mv', mon_out + '.root', 'esepp.root'])

    if args.monitor_rad:
        subprocess.call([run_cmd, monrad_exe, data_files])
        subprocess.call([run_cmd, convert_exe, monrad_out])
        subprocess.call(['mv', monrad_out + '.root', 'rga-rad.root'])

        subprocess.call([run_cmd, monrad_exe, mc_files])
        subprocess.call([run_cmd, convert_exe, monrad_out])
        subprocess.call(['mv', monrad_out + '.root', 'esepp-rad.root'])
