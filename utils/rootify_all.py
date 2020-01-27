#!/usr/bin/env python 

import glob 
import subprocess 

if __name__ == "__main__":

    files = glob.glob('*.hipo')

    for f in files:
        subprocess.call(['run-groovy', 'j2root.groovy', f])
