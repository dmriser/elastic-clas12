import event.Event
import event.EventConverter
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.io.hipo.HipoDataSource

// Additional class members to make the object
// more useful for CLAS12 analysis.
Particle.metaClass.pindex = null
Particle.metaClass.sphi = null
Particle.metaClass.sector = null

def beam = new Particle(11, 0.0, 0.0, 10.646)
def target = new Particle(2212, 0.0, 0.0, 0.0)

cuts = [
        w: [0.8, 1.15],
        w_loose: [0.8, 1.30],
        angle: [178, 180]
]

new_kin_bounds = [
        theta_ele : [4, 20],
        theta_pro : [20, 70],
        p_ele     : [7, 10.5],
        p_pro     : [0.5, 4.5],
        w         : [0.6, 2.2],
        q2        : [1.2, 4.5],
        phi       : [-30, 330],
        vz        : [-20, 15],
        dp_ele    : [-1.2, 1.2],
        dp_pro    : [-1.2, 1.2],
        dtheta_pro: [-20, 20],
        de_beam   : [-2, 2],
        angle_ep  : [125, 180]
]

outbend_bounds = [
        theta_ele : [4, 20],
        theta_pro : [20, 70],
        p_ele     : [7, 10.5],
        p_pro     : [0.5, 4.5],
        w         : [0.6, 1.7],
        q2        : [0.5, 4.5],
        phi       : [-30, 330],
        vz        : [-20, 15],
        dp_ele    : [-1.2, 1.2],
        dp_pro    : [-1.2, 1.2],
        dtheta_pro: [-20, 20],
        de_beam   : [-2, 2],
        angle_ep  : [125, 180]
]

rgk_kin_bounds = [
        theta_ele : [5, 30],
        theta_pro : [20, 70],
        p_ele     : [4, 8.5],
        p_pro     : [0.5, 4.5],
        w         : [0.6, 1.7],
        q2        : [1.2, 4.5],
        phi       : [-30, 330],
        vz        : [-20, 15],
        dp_ele    : [-1.2, 1.2],
        dp_pro    : [-1.2, 1.2],
        dtheta_pro: [-20, 20],
        de_beam   : [-2, 2],
        angle_ep  : [125, 180]
]

tighter_kin_bounds = [
        theta_ele : [5, 15],
        theta_pro : [30, 70],
        p_ele     : [7.8, 10.5],
        p_pro     : [0.5, 3.5],
        w         : [0.6, 1.7],
        phi       : [-30, 330],
        dp_ele    : [-0.6, 0.6],
        dp_pro    : [-1.6, 1.6],
        dtheta_ele: [-2, 2],
        dtheta_pro: [-4, 4],
        angle_ep  : [160, 180],
        fracp_ele : [-0.1, 0.1],
        fracp_pro : [-0.5, 0.5],
        q2        : [1.2, 4.5],
        vz        : [-20, 15],
        de_beam   : [-2, 2]
]

lim = tighter_kin_bounds

def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}

histos = new ConcurrentHashMap()
histoBuilders = [
        w        : { title -> limited_h1(title, 200, lim.w) },
        theta_res: { title -> limited_h1(title, 200, lim.dtheta_pro) },
        p_res    : { title -> limited_h1(title, 200, lim.dp_ele) },
        p_ele    : { title -> limited_h1(title, 200, lim.p_ele) },
        p_pro    : { title -> limited_h1(title, 200, lim.p_pro) },
        vz       : { title -> limited_h1(title, 200, lim.vz) },
        de_beam  : { title -> limited_h1(title, 200, lim.de_beam) },
        angle_ep : { title -> limited_h1(title, 200, lim.angle_ep) },
        theta_p  : { title -> limited_h1(title, 200, lim.theta_pro) },
        theta_ele: { title -> limited_h1(title, 200, lim.theta_ele) }
]

histoBuilders2 = [
        de_beam_de_beam  : { title -> limited_h2(title, 200, 200, lim.de_beam, lim.de_beam) },
        p_pro_dp         : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dp_pro) },
        p_ele_dp         : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dp_ele) },
        p_ele_theta      : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.theta_ele) },
        p_pro_dtheta     : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dtheta_pro) },
        p_ele_dtheta     : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dtheta_ele) },
        p_ele_fracp      : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.fracp_ele) },
        p_pro_fracp      : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.fracp_pro) },
        p_w_ele          : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.w) },
        phi_vz           : { title -> limited_h2(title, 200, 200, lim.phi, lim.vz) },
        phi_dp           : { title -> limited_h2(title, 100, 100, lim.phi, lim.dp_ele) },
        phi_theta        : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_ele) },
        phi_w            : { title -> limited_h2(title, 200, 200, lim.phi, lim.w) },
        phi_theta_proton : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_pro) },
        theta_ele_de_beam: { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.de_beam) },
        theta_pro_de_beam: { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.de_beam) },
        theta_ele_dp     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dp_ele) },
        theta_ele_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dtheta_ele) },
        theta_pro_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dtheta_pro) },
        theta_pro_dp     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dp_ele) },
        theta_ele_vz     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.vz) },
        theta_pro_vz     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.vz) },
        theta_w_ele      : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.w) },
        theta_ele_dp     : { title -> limited_h2(title, 100, 100, lim.theta_ele, lim.dp_ele) },
        theta_pro_dp     : { title -> limited_h2(title, 100, 100, lim.theta_pro, lim.dp_pro) },
        w_q2             : { title -> limited_h2(title, 200, 200, lim.w, lim.q2) },
]


def getGeneratedSector(phi){
    return Math.ceil(phi / 60) + 3
}

def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def getKin(beam, target, electron) {

    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    return [x: x, y: y, w: w, nu: nu, q2: q2]
}

def getPKin(beam, target, electron, proton) {

    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()
    missing.combine(proton, -1)
    def missing_mass = missing.mass2()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    def zaxis = new Vector3(0, 0, 1)
    def enorm = electron.vector().vect().cross(zaxis)
    def pnorm = proton.vector().vect().cross(zaxis)
    def phi = enorm.theta(pnorm)

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e()]
}

def predictElectronThetaFromMomentum(beam, mom) {
    def pred_theta = Math.acos(1 - PDGDatabase.getParticleMass(2212) * ((beam.e() - mom) / (beam.e() * mom)))
    return pred_theta
}

def predictElasticBasedOnElectronAngle(beam, alpha) {
    // alpha - electron angle in radians (known)
    // beta  - proton angle in radians (inferred)
    //
    // Don't confuse the beta in this function with the
    // physics beta, v/c.

    // This is stupid, but it makes the next lines shorter
    def proton_mass = PDGDatabase.getParticleMass(2212)
    def cos_alpha = Math.cos(alpha)
    def sin_alpha = Math.sin(alpha)

    // Kinematics prediction
    def pred_ele_p = beam.e() / (1 + (beam.e() / proton_mass) * (1 - cos_alpha))
    def pred_pro_p = Math.sqrt(beam.e()**2 - 2 * beam.e() * pred_ele_p * cos_alpha + pred_ele_p**2)
    def pred_beta = Math.asin((pred_ele_p * sin_alpha) / pred_pro_p)

    return [pred_ele_p, pred_beta, pred_pro_p]
}

def predictBeamEnergy(ele, pro){
    def pred_e_beam = ele.p() / (1 + (ele.p() / PDGDatabase.getParticleMass(2212)) * (Math.cos(ele.theta()) - 1))
    def a0 = 1 - 1 / (Math.cos(ele.theta()) - Math.sin(ele.theta()) / Math.tan(-1 * pro.theta()))
    def a1 = Math.sin(ele.theta()) / Math.sin(ele.theta() + pro.theta())
    def pred_e_beam_from_angles = 2 * PDGDatabase.getParticleMass(2212) * a0 / (a1**2 - a0**2)
    return [pred_e_beam, pred_e_beam_from_angles]
}

// Histogram filling methods
def fillInclusiveHistograms(histos, kin, sector){
    histos.computeIfAbsent('w_inclusive_' + sector, histoBuilders.w).fill(kin.w)
    histos.computeIfAbsent('w_q2_inclusive_' + sector, histoBuilders2.w_q2).fill(kin.w, kin.q2)
}

def fillProtonResolutions(histos, event, ele, pro, pred_pro_theta, pred_pro_p, delta_beam_energy){

    histos.computeIfAbsent("theta_proton_" + "_" + ele.sector, histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
    histos.computeIfAbsent("theta_electron_" + "_" + ele.sector, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
    histos.computeIfAbsent("p_ele_" + "_" + ele.sector, histoBuilders.p_ele).fill(ele.p())
    histos.computeIfAbsent("p_pro_" + "_" + ele.sector, histoBuilders.p_pro).fill(pro.p())
    histos.computeIfAbsent('delta_p_proton_' + ele.sector, histoBuilders.p_res).fill(pro.p() - pred_pro_p)

    histos.computeIfAbsent('delta_theta_proton_' + ele.sector, histoBuilders.theta_res).fill(
            Math.toDegrees(pro.theta() - pred_pro_theta))

    histos.computeIfAbsent('phi_proton_vz_proton', histoBuilders2.phi_vz).fill(
            pro.sphi, event.vz[pro.pindex])

    histos.computeIfAbsent('phi_proton_delta_vz', histoBuilders2.phi_vz).fill(
            pro.sphi, event.vz[ele.pindex] - event.vz[pro.pindex])

    histos.computeIfAbsent('vz_electron_' + ele.sector, histoBuilders.vz).fill(event.vz[ele.pindex])
    histos.computeIfAbsent('vz_proton_' + ele.sector, histoBuilders.vz).fill(event.vz[pro.pindex])
    histos.computeIfAbsent('delta_vz_' + ele.sector, histoBuilders.vz).fill(event.vz[ele.pindex] - event.vz[pro.pindex])

    histos.computeIfAbsent('theta_proton_delta_p_proton_' + ele.sector,
            histoBuilders2.theta_pro_dp).fill(Math.toDegrees(pro.theta()), pro.p() - pred_pro_p)
    histos.computeIfAbsent('theta_proton_de_beam_' + ele.sector,
            histoBuilders2.theta_pro_dp).fill(Math.toDegrees(pro.theta()), delta_beam_energy)

    histos.computeIfAbsent('theta_proton_delta_theta_proton_' + ele.sector,
            histoBuilders2.theta_pro_dtheta).fill(
            Math.toDegrees(pro.theta()), Math.toDegrees(pro.theta() - pred_pro_theta))
    histos.computeIfAbsent('theta_proton_vz_proton_' + ele.sector, histoBuilders2.theta_pro_vz).fill(
            Math.toDegrees(pro.theta()), event.vz[pro.pindex])
    histos.computeIfAbsent('p_proton_delta_p_proton_' + ele.sector, histoBuilders2.p_pro_dp).fill(
            event.p[pro.pindex], event.p[pro.pindex] - pred_pro_p)

}

def fillElectronResolutions(histos, event, beam, ele, pro, pred_e_beam,
                            pred_e_beam_from_angles, pred_ele_p, pred_ele_theta){

    // One dimensional
    histos.computeIfAbsent('delta_p_electron_' + ele.sector, histoBuilders.p_res).fill(ele.p() - pred_ele_p)
    histos.computeIfAbsent('de_beam_' + ele.sector, histoBuilders.de_beam).fill(beam.e() - pred_e_beam)
    histos.computeIfAbsent('de_beam_from_angles' + ele.sector, histoBuilders.de_beam).fill(
            beam.e() - pred_e_beam_from_angles)

    // Two dimensional
    histos.computeIfAbsent('p_electron_delta_p_electron_' + ele.sector, histoBuilders2.p_ele_dp).fill(
            ele.p(), ele.p() - pred_ele_p
    )

    histos.computeIfAbsent('phi_electron_theta_electron', histoBuilders2.phi_theta).fill(
            ele.sphi, Math.toDegrees(ele.theta()))

    histos.computeIfAbsent('phi_electron_delta_p_electron', histoBuilders2.phi_dp).fill(
            ele.sphi, ele.p() - pred_ele_p)

    histos.computeIfAbsent('phi_electron_delta_vz', histoBuilders2.phi_vz).fill(
            ele.sphi, event.vz[ele.pindex] - event.vz[pro.pindex])

    histos.computeIfAbsent('theta_electron_vz_electron_' + ele.sector, histoBuilders2.theta_ele_vz).fill(
            Math.toDegrees(ele.theta()), event.vz[ele.pindex])

    histos.computeIfAbsent('theta_electron_delta_p_electron_' + ele.sector,
            histoBuilders2.theta_ele_dp).fill(Math.toDegrees(ele.theta()), ele.p() - pred_ele_p)

    histos.computeIfAbsent('theta_electron_delta_theta_electron_' + ele.sector,
            histoBuilders2.theta_ele_dtheta).fill(
            Math.toDegrees(ele.theta()), Math.toDegrees(ele.theta() - pred_ele_theta))

    histos.computeIfAbsent('phi_electron_vz_electron', histoBuilders2.phi_vz).fill(
            ele.sphi, event.vz[ele.pindex])
    histos.computeIfAbsent('phi_electron_vz_proton', histoBuilders2.phi_vz).fill(
            ele.sphi, event.vz[pro.pindex])
    histos.computeIfAbsent('theta_ele_de_beam_' + ele.sector, histoBuilders2.theta_ele_de_beam).fill(
            Math.toDegrees(ele.theta()), beam.e() - pred_e_beam)
    histos.computeIfAbsent('theta_ele_de_beam_from_angles_' + ele.sector, histoBuilders2.theta_ele_de_beam).fill(
            Math.toDegrees(ele.theta()), beam.e() - pred_e_beam_from_angles)
    histos.computeIfAbsent('de_beam_de_beam_from_angles' + ele.sector, histoBuilders2.de_beam_de_beam).fill(
            beam.e() - pred_e_beam, beam.e() - pred_e_beam_from_angles)

}

def fillSimulationResolutions(histos, sector, ele, pro, gen_ele, gen_pro){
    
    histos.computeIfAbsent('p_electron_dp_electron_simulation_' + sector, histoBuilders2.p_ele_dp).fill(
	gen_ele.p(), ele.p() - gen_ele.p()
    )
    histos.computeIfAbsent('theta_electron_dtheta_electron_simulation_' + sector, histoBuilders2.theta_ele_dtheta).fill(
	Math.toDegrees(gen_ele.theta()), Math.toDegrees(ele.theta() - gen_ele.theta())
    )
    histos.computeIfAbsent('p_proton_dp_proton_simulation_' + sector, histoBuilders2.p_pro_dp).fill(
	gen_pro.p(), pro.p() - gen_pro.p()
    )
    histos.computeIfAbsent('theta_proton_dtheta_proton_simulation_' + sector, histoBuilders2.theta_pro_dtheta).fill(
	Math.toDegrees(gen_pro.theta()), Math.toDegrees(pro.theta() - gen_pro.theta())
    )

}

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent()) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

	    // Reconstructed (for data and simulation)
            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
                def sector = event.dc_sector[idx]
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                ele.sector = sector
                ele.pindex = idx
                ele.sphi = shiftPhi(Math.toDegrees(ele.phi()))
                def kin = getKin(beam, target, ele)

                fillInclusiveHistograms(histos, kin, sector)

                // Predict kinematics based on electron measured variables.
                def (pred_ele_p, pred_pro_theta, pred_pro_p) = predictElasticBasedOnElectronAngle(beam, ele.theta())
                def pred_ele_theta = predictElectronThetaFromMomentum(beam, ele.p())

                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    pro.pindex = it
                    pro.sphi = shiftPhi(Math.toDegrees(pro.phi()))
                    def pkin = getPKin(beam, target, ele, pro)

                    def (pred_e_beam, pred_e_beam_from_angles) = predictBeamEnergy(ele, pro)

                    // One dimensional
                    histos.computeIfAbsent('w_' + sector, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('w', histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('angle_ep', histoBuilders.angle_ep).fill(pkin.angle)
                    histos.computeIfAbsent('angle_ep_' + sector, histoBuilders.angle_ep).fill(pkin.angle)

                    // For illustration of selection criteria
                    if (event.ctof_status.contains(it)) {
                        histos.computeIfAbsent('w_in_ctof', histoBuilders.w).fill(pkin.w)
                    }

                    // Pass phi but no cut on W
                    if (pkin.angle > cuts.angle[0] && event.ctof_status.contains(it)) {
                        histos.computeIfAbsent('w_pass_angle_in_ctof', histoBuilders.w).fill(pkin.w)
                        histos.computeIfAbsent('w_pass_angle_in_ctof_' + sector, histoBuilders.w).fill(pkin.w)
                        histos.computeIfAbsent('w_q2_pass_angle_in_ctof_' + sector, histoBuilders2.w_q2).fill(pkin.w, pkin.q2)
                        histos.computeIfAbsent('p_w_ele_' + sector, histoBuilders2.p_w_ele).fill(ele.p(), pkin.w)
                        histos.computeIfAbsent('theta_w_ele_' + sector, histoBuilders2.theta_w_ele).fill(
                                Math.toDegrees(ele.theta()), pkin.w)
			histos.computeIfAbsent('p_theta_ele_pass_angle_in_ctof_' + sector, histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta())
			)
                    }

                    // Pass W but no cut on phi
                    if (pkin.w > cuts.w[0] && pkin.w < cuts.w[1] && event.ctof_status.contains(it)) {
                        histos.computeIfAbsent('angle_ep_pass_w_in_ctof', histoBuilders.angle_ep).fill(pkin.angle)
                        histos.computeIfAbsent('angle_ep_pass_w_in_ctof_' + sector, histoBuilders.angle_ep).fill(pkin.angle)
                    }

                    // Elastic protons in forward and central.
                    if (pkin.w > cuts.w[0] && pkin.w < cuts.w[1] && pkin.angle > cuts.angle[0]) {
                        histos.computeIfAbsent('theta_p_combined', histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))

                        if (event.ctof_status.contains(it)) {
                            histos.computeIfAbsent('theta_p_ctof', histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
                        }

                        // layers = {1: ftof1a, 2: ftof1b, 3: ftof2}
                        else if (event.tof_status.contains(it)) {
                            def layer = event.tof_layer[it]
                            histos.computeIfAbsent('theta_p_tof_' + layer, histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
                        }

                    }

                    histos.computeIfAbsent('w_q2_' + sector, histoBuilders2.w_q2).fill(pkin.w, pkin.q2)
                    histos.computeIfAbsent('phi_electron_w', histoBuilders2.phi_w).fill(ele.sphi, pkin.w)

                    if (pkin.angle > cuts.angle[0] && pkin.w < cuts.w_loose[1] && event.ctof_status.contains(it)) {
                        fillElectronResolutions(histos, event, beam, ele, pro, pred_e_beam,
                                pred_e_beam_from_angles, pred_ele_p, pred_ele_theta)

                        // We can go tight on protons
                        if (pkin.w > cuts.w[0] && pkin.w < cuts.w[1]) {
                            fillProtonResolutions(histos, event, ele, pro, pred_pro_theta, pred_pro_p, beam.e() - pred_e_beam_from_angles)


			    // This is a good elastic event, let's see if this is simulation.
			    // If so, find the generated partciles and make more histograms.
			    (0..<event.mc_npart).find { j ->
				event.mc_pid[j] == 11
			    }?.each { gidx ->
				def gen_ele = new Particle(11, event.mc_px[gidx], event.mc_py[gidx], event.mc_pz[gidx])
				def gen_sector = (int) getGeneratedSector(Math.toDegrees(gen_ele.phi()))
				
				(0..<event.mc_npart).find { k ->
				    event.mc_pid[k] == 2212
				}?.each { pidx ->
				    def gen_pro = new Particle(2212, event.mc_px[pidx], event.mc_py[pidx], event.mc_pz[pidx])
				    def gen_pkin = getPKin(beam, target, gen_ele, gen_pro)
				    fillSimulationResolutions(histos, sector, ele, pro, gen_ele, gen_pro)
				}
			    }
                        }
                    }
                }
            }

            eventIndex++
        }
    }
}

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("monitor.hipo")
