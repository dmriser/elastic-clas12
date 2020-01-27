"""
A utility to convert the ROOT files output 
from the event generator into LUND format. 

"""

from collections import namedtuple

import argparse 
import numpy as np 

from ROOT import TFile, TNtuple, TTree

PROTON_MASS = 0.938

LundHeader = namedtuple(
    'LundHeader', 
    "npart mass_target z_target pol_target pol_beam beam_type beam_energy nucleon_id process_id event_weight"
)

LundParticle = namedtuple(
    'LundParticle',
    "index lifetime type_id particle_id parent first_daughter px py pz energy mass vx vy vz"
)

def get_px_py_pz(p, theta, phi):
    """ Getting particle momenta. """
    pz = p * np.cos(theta)
    px = p * np.sin(theta) * np.cos(phi)
    py = p * np.sin(theta) * np.sin(phi)
    return px, py, pz

def get_particles(event, vertex):
    """ Return LUND particles for e, p, g. """
    ele_p = event.E_l / 1000.0
    pro_p = np.sqrt(event.E_p**2 - (1000 * PROTON_MASS)**2) / 1000.0
    gam_p = event.E_g / 1000.0

    # Assume point vertex 
    px, py, pz = get_px_py_pz(ele_p, event.theta_l, event.phi_l)
    ele = LundParticle(
        index=1, lifetime=0, type_id=1, particle_id=11, 
        px=px, py=py, pz=pz, energy=event.E_l/1000.0,
        mass=0, vx=vertex[0], vy=vertex[1], vz=vertex[2], 
        parent=0, first_daughter=0
    )

    px, py, pz = get_px_py_pz(pro_p, event.theta_p, event.phi_p)
    pro = LundParticle(
        index=2, lifetime=0, type_id=1, particle_id=2212, 
        px=px, py=py, pz=pz, energy=event.E_p/1000.0,
        mass=PROTON_MASS, vx=vertex[0], vy=vertex[1], vz=vertex[2], 
        parent=0, first_daughter=0
    )

    px, py, pz = get_px_py_pz(gam_p, event.theta_g, event.phi_g)
    gam = LundParticle(
        index=3, lifetime=0, type_id=1, particle_id=22, 
        px=px, py=py, pz=pz, energy=event.E_g/1000.0,
        mass=0, vx=vertex[0], vy=vertex[1], vz=vertex[2],  
        parent=0, first_daughter=0
    )
    return ele, pro, gam

def get_dummy_header(npart=3):
    return LundHeader(
        npart=npart,
        mass_target=1,
        z_target=1,
        pol_target=0,
        pol_beam=0,
        beam_type=11,
        beam_energy=10.646,
        nucleon_id=0,
        process_id=0,
        event_weight=1
    )

def particle_to_string(part):
    return ' {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(
        part.index, part.lifetime, part.type_id, part.particle_id,
        part.parent, part.first_daughter, part.px, part.py, part.pz,
        part.energy, part.mass, part.vx, part.vy, part.vz
    )
    
def header_to_string(header):
    return ' {} {} {} {} {} {} {} {} {} {}\n'.format(
        header.npart, header.mass_target, header.z_target, 
        header.pol_target, header.pol_beam, header.beam_type,
        header.beam_energy, header.nucleon_id, 
        header.process_id, header.event_weight
    )

if __name__ == '__main__':

    ap = argparse.ArgumentParser() 
    ap.add_argument('-i', '--input_file', required=True)
    ap.add_argument('-o', '--output_file_prefix', required=True)
    ap.add_argument('-b', '--batch_size', type=int, default=1000, 
                    help='Events per LUND file.')
    args = ap.parse_args()

    rootfile = TFile(args.input_file)
    ntuple = rootfile.Get('ntp')

    file_id = 0 
    current_file = args.output_file_prefix + str(file_id) + '.lund'
    output = open(current_file, 'w')
    for iev, ev in enumerate(ntuple):
        
        # Generate RG-A vertex 
        vertex = np.array([0.0, 0.0, np.random.uniform(-5, 0)])

        # Get header and particles 
        header = get_dummy_header()
        ele, pro, gam = get_particles(ev, vertex)

        output.write(header_to_string(header))
        output.write(particle_to_string(ele))
        output.write(particle_to_string(pro))
        output.write(particle_to_string(gam))

        if iev != 0 and iev % args.batch_size == 0:
            print('Finished {} events ({} files).'.format(iev, file_id))
                
            file_id += 1
            current_file = args.output_file_prefix + str(file_id) + '.txt'
            output.close() 
            output = open(current_file, 'w')

    output.close() 
