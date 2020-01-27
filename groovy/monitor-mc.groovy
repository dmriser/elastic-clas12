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


// Kinematics from RG-A
def beam = new Particle(11, 0.0, 0.0, 10.646)
//def beam = new Particle(11, 0.0, 0.0, 7.5)
def target = new Particle(2212, 0.0, 0.0, 0.0)

orig_kin_bounds = [
        theta_ele : [5, 30],
        theta_pro : [50, 90],
        p_ele     : [0, 11],
        p_pro     : [0, 3],
        w         : [0.8, 1.7],
        q2        : [0, 5],
        phi       : [-30, 330],
        vz        : [-20, 15],
        dp_ele    : [-0.2, 0.2],
        dp_pro    : [-0.2, 0.2],
        dtheta_pro: [-10, 15],
        de_beam   : [-2, 2],
        angle_ep  : [125, 180]
]

new_kin_bounds = [
        theta_ele : [5, 20],
        theta_pro : [20, 70],
        p_ele     : [7, 10.5],
        p_pro     : [0.5, 4.5],
        w         : [0.6, 1.7],
        q2        : [1.2, 4.5],
        phi       : [-30, 330],
        vz        : [-20, 15],
        dp_ele    : [-1.2, 1.2],
        dp_pro    : [-1.2, 1.2],
        dtheta_pro: [-20, 20],
        de_beam    : [-2, 2],
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
        de_beam    : [-2, 2],
        angle_ep  : [125, 180]
]

lim = new_kin_bounds

def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}

histos = new ConcurrentHashMap()

histoBuilders = [
        w         : { title -> limited_h1(title, 400, lim.w) },
        theta_res : { title -> limited_h1(title, 400, lim.dtheta_pro) },
        p_res     : { title -> limited_h1(title, 400, lim.dp_ele) },
        vz        : { title -> limited_h1(title, 400, lim.vz) },
        de_beam   : { title -> limited_h1(title, 400, lim.de_beam) },
        angle_ep  : { title -> limited_h1(title, 400, lim.angle_ep) },
        theta_p   : { title -> new H1F("$title", "$title", 400, 5, 85) },
        theta_ele : { title -> limited_h1(title, 400, lim.theta_ele )}
]

histoBuilders2 = [
        w_p_ele           : { title -> limited_h2(title, 200, 200, lim.w, lim.p_ele) },
        w_q2              : { title -> limited_h2(title, 200, 200, lim.w, lim.q2) },
        phi_w             : { title -> limited_h2(title, 200, 200, lim.phi, lim.w) },
        theta_ele_vz      : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.vz) },
        phi_vz            : { title -> limited_h2(title, 200, 200, lim.phi, lim.vz) },
        theta_ele_dp      : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dp_ele) },
        theta_ele_dtheta  : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dtheta_pro) },
        theta_pro_dtheta  : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dtheta_pro) },
        theta_pro_dp      : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dp_ele) },
        theta_pro_vz      : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.vz) },
        phi_dp            : { title -> limited_h2(title, 200, 200, lim.phi, lim.dp_ele) },
        phi_theta         : { title -> limited_h2(title, 200, 200, lim.phi, lim.theta_ele) },
        phi_theta_proton  : { title -> limited_h2(title, 200, 200, lim.phi, lim.theta_pro) },
        p_pro_dp          : { title -> limited_h2(title, 200, 200, lim.p_pro, lim.dp_ele) },
        theta_ele_de_beam : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.de_beam) },
        de_beam_de_beam   : { title -> limited_h2(title, 200, 200, lim.de_beam, lim.de_beam) },
]


def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def relativePhi(phi, sector) {
    if (sector == 4) {
        return phi
    } else if (sector == 5) {
        return phi - 60
    } else if (sector == 6) {
        return phi - 120
    } else if (sector == 1) {
        return phi - 180
    } else if (sector == 2) {
        return phi - 240
    } else if (sector == 3) {
        return phi - 300
    }
}


def angleBetween(v1, v2) {
    v1.unit()
    v2.unit()
    return Math.toDegrees(
            Math.acos(v1.dot(v2))
    )
}

def getDeltaVertex(event, i, j) {
    return Math.sqrt((event.vx[j] - event.vx[i])**2 + (event.vy[j] - event.vy[i])**2 + (event.vz[j] - event.vz[i])**2)
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

    def missing_ele = new Particle(beam)
    missing_ele.combine(target, 1)
    missing_ele.combine(proton, -1)

    def missing_pro = new Particle(beam)
    missing_pro.combine(target, 1)
    missing_pro.combine(electron, -1)

    def dtheta_ele = Math.toDegrees(electron.theta() - missing_ele.theta())
    def dp_ele = electron.p() - missing_ele.p()
    def dtheta_pro = Math.toDegrees(proton.theta() - missing_pro.theta())
    def dp_pro = proton.p() - missing_pro.p()

    return [x           : x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e(),
            dtheta_ele  : dtheta_ele, dtheta_pro: dtheta_pro,
            dp_ele      : dp_ele, dp_pro: dp_pro]
}

def getElectronDeltas(beam, ele) {
    def calc_energy = beam.e() / (1 + (beam.e() / PDGDatabase.getParticleMass(2212))) * (1 - Math.cos(ele.theta()))
    def delta_energy = calc_energy - ele.e()
    def calc_theta = Math.toDegrees(Math.acos(1 + (PDGDatabase.getParticleMass(2212) / beam.e()) * (1 - beam.e() / ele.e())))
    def delta_theta = calc_theta - Math.toDegrees(ele.theta())
    return [delta_theta, delta_energy]
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

	    // Generated Stuff 
	    (0 ..< event.mc_npart).find{
		event.mc_pid[it] == 11
	    }?.each{ idx -> 
		def ele = new Particle(11, event.mc_px[idx], event.mc_py[idx], event.mc_pz[idx])
		def kin = getKin(beam, target, ele)
		def sphi_ele = shiftPhi(Math.toDegrees(ele.phi()))

		(0 ..< event.mc_npart).find{
		    event.mc_pid[it] == 2212
		}?.each{ pidx -> 
		    def pro = new Particle(2212, event.mc_px[pidx], event.mc_py[pidx], event.mc_pz[pidx])		    
                    def pkin = getPKin(beam, target, ele, pro)
                    def phi_pro = Math.toDegrees(pro.phi())
                    def sphi_pro = shiftPhi(phi_pro)
                    def pred_e_beam = ele.p() / (1 + (ele.p() / PDGDatabase.getParticleMass(2212)) * (Math.cos(ele.theta()) - 1))
                    def a0 = 1 - 1 / (Math.cos(ele.theta()) - Math.sin(ele.theta()) / Math.tan(-1 * pro.theta()))
                    def a1 = Math.sin(ele.theta()) / Math.sin(ele.theta() + pro.theta())
                    def pred_e_beam_from_angles = 2 * PDGDatabase.getParticleMass(2212) * a0 / (a1**2 - a0**2)

		    histos.computeIfAbsent('w_q2_gen', histoBuilders2.w_q2).fill(pkin.w, pkin.q2)
		    histos.computeIfAbsent('phi_proton_theta_proton_gen', histoBuilders2.phi_theta_proton).fill(
			sphi_pro, Math.toDegrees(pro.theta())
		    )
		    histos.computeIfAbsent('phi_electron_theta_electron_gen', histoBuilders2.phi_theta).fill(
			sphi_ele, Math.toDegrees(ele.theta())
		    )
		    histos.computeIfAbsent('theta_proton_gen', histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
		    histos.computeIfAbsent('theta_electron_gen', histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
		    histos.computeIfAbsent('angle_ep_gen', histoBuilders.angle_ep).fill(pkin.angle)
		}

	    }
	    

	    // Reconstructed Stuff 
            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                def kin = getKin(beam, target, ele)
                def sector = event.dc_sector[idx]
                def phi = Math.toDegrees(ele.phi())
                def sphi = shiftPhi(phi)

                // Some histograms need to be filled with all electrons, even if there
                // are no postive tracks in the central.
                histos.computeIfAbsent('w_inclusive_' + sector, histoBuilders.w).fill(kin.w)
                histos.computeIfAbsent('w_q2_inclusive_' + sector, histoBuilders2.w_q2).fill(kin.w, kin.q2)

                // Believe the electron angular measurement and infer the rest of the
                // event kinematics based on the assumption that it's an elastic scattering.
                def (pred_ele_p, pred_pro_theta, pred_pro_p) = predictElasticBasedOnElectronAngle(beam, ele.theta())

                (0..<event.npart).findAll { event.charge[it] > 0 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    def pkin = getPKin(beam, target, ele, pro)
                    def phi_pro = Math.toDegrees(pro.phi())
                    def sphi_pro = shiftPhi(phi_pro)
                    def pred_e_beam = ele.p() / (1 + (ele.p() / PDGDatabase.getParticleMass(2212)) * (Math.cos(ele.theta()) - 1))

                    def a0 = 1 - 1 / (Math.cos(ele.theta()) - Math.sin(ele.theta()) / Math.tan(-1 * pro.theta()))
                    def a1 = Math.sin(ele.theta()) / Math.sin(ele.theta() + pro.theta())
                    def pred_e_beam_from_angles = 2 * PDGDatabase.getParticleMass(2212) * a0 / (a1**2 - a0**2)

                    // One dimensional
                    histos.computeIfAbsent('w_' + sector, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('w', histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('angle_ep', histoBuilders.angle_ep).fill(pkin.angle)


		    // For illustration of selection criteria
                    if (event.ctof_status.contains(it)) {
                        histos.computeIfAbsent('w_in_ctof', histoBuilders.w).fill(pkin.w)
                    }
                    if (pkin.angle > 0 && event.ctof_status.contains(it)) {
                        histos.computeIfAbsent('w_pass_angle_in_ctof', histoBuilders.w).fill(pkin.w)
                    }
                    if (pkin.w > 0.8 && pkin.w < 1.08 && event.ctof_status.contains(it)) {
                        histos.computeIfAbsent('angle_ep_pass_w_in_ctof', histoBuilders.angle_ep).fill(pkin.angle)
                    }

                    // Elastic protons in forward and central.
                    if (pkin.w > 0.8 && pkin.w < 1.08 && pkin.angle > 0){
                        histos.computeIfAbsent('theta_p_combined', histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))

                        if (event.ctof_status.contains(it)){
                            histos.computeIfAbsent('theta_p_ctof', histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
                        }

                        // layers = {1: ftof1a, 2: ftof1b, 3: ftof2}
                        else if (event.tof_status.contains(it)) {
                            def layer = event.tof_layer[it]
                            histos.computeIfAbsent('theta_p_tof_' + layer, histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
                        }

                    }

                    // Two dimensional
                    histos.computeIfAbsent('w_q2_' + sector, histoBuilders2.w_q2).fill(pkin.w, pkin.q2)
                    histos.computeIfAbsent('phi_electron_w', histoBuilders2.phi_w).fill(sphi, pkin.w)

                    // Require that the proton is in central detector.
                    if (pkin.angle > 0 && pkin.w < 1.3 && event.ctof_status.contains(it)) {

                        // One dimensional
                        histos.computeIfAbsent('delta_p_electron_' + sector, histoBuilders.p_res).fill(ele.p() - pred_ele_p)
                        histos.computeIfAbsent('delta_p_proton_' + sector, histoBuilders.p_res).fill(pro.p() - pred_pro_p)
                        histos.computeIfAbsent('delta_theta_proton_' + sector, histoBuilders.theta_res).fill(
                                Math.toDegrees(pro.theta() - pred_pro_theta))
                        histos.computeIfAbsent('vz_electron_' + sector, histoBuilders.vz).fill(event.vz[idx])
                        histos.computeIfAbsent('vz_proton_' + sector, histoBuilders.vz).fill(event.vz[it])
                        histos.computeIfAbsent('delta_vz_' + sector, histoBuilders.vz).fill(event.vz[idx] - event.vz[it])
                        histos.computeIfAbsent('de_beam_' + sector, histoBuilders.de_beam).fill(beam.e() - pred_e_beam)
                        histos.computeIfAbsent('de_beam_from_angles' + sector, histoBuilders.de_beam).fill(
                                beam.e() - pred_e_beam_from_angles)

                        // Two dimensional
			histos.computeIfAbsent('w_p_ele_' + sector, histoBuilders2.w_p_ele).fill(pkin.w, ele.p())
                        histos.computeIfAbsent('phi_electron_theta_electron', histoBuilders2.phi_theta).fill(
                                sphi, Math.toDegrees(ele.theta()))
                        histos.computeIfAbsent('phi_electron_delta_p_electron', histoBuilders2.phi_dp).fill(
                                sphi, ele.p() - pred_ele_p)
                        histos.computeIfAbsent('phi_electron_delta_vz', histoBuilders2.phi_vz).fill(sphi, event.vz[idx] - event.vz[it])
                        histos.computeIfAbsent('theta_electron_vz_electron_' + sector, histoBuilders2.theta_ele_vz).fill(
                                Math.toDegrees(ele.theta()), event.vz[idx])
                        histos.computeIfAbsent('theta_electron_delta_p_electron_' + sector,
                                histoBuilders2.theta_ele_dp).fill(Math.toDegrees(ele.theta()), ele.p() - pred_ele_p)
                        histos.computeIfAbsent('theta_electron_delta_theta_proton_' + sector,
                                histoBuilders2.theta_ele_dtheta).fill(
                                Math.toDegrees(ele.theta()), Math.toDegrees(pro.theta() - pred_pro_theta))
                        histos.computeIfAbsent('theta_proton_delta_p_proton_' + sector,
                                histoBuilders2.theta_pro_dp).fill(Math.toDegrees(pro.theta()), pro.p() - pred_pro_p)

                        histos.computeIfAbsent('theta_proton_delta_theta_proton_' + sector,
                                histoBuilders2.theta_pro_dtheta).fill(
                                Math.toDegrees(pro.theta()), Math.toDegrees(pro.theta() - pred_pro_theta))
                        histos.computeIfAbsent('theta_proton_vz_proton_' + sector, histoBuilders2.theta_pro_vz).fill(
                                Math.toDegrees(pro.theta()), event.vz[it])
                        histos.computeIfAbsent('p_proton_delta_p_proton_' + sector, histoBuilders2.p_pro_dp).fill(
                                event.p[it], event.p[it] - pred_pro_p)
                        histos.computeIfAbsent('phi_electron_vz_electron', histoBuilders2.phi_vz).fill(
                                sphi, event.vz[idx])
                        histos.computeIfAbsent('phi_electron_vz_proton', histoBuilders2.phi_vz).fill(
                                sphi, event.vz[it])
                        histos.computeIfAbsent('phi_proton_vz_proton', histoBuilders2.phi_vz).fill(
                                sphi_pro, event.vz[it])
                        histos.computeIfAbsent('phi_proton_delta_vz', histoBuilders2.phi_vz).fill(
                                sphi_pro, event.vz[idx] - event.vz[it])
                        histos.computeIfAbsent('theta_ele_de_beam_' + sector, histoBuilders2.theta_ele_de_beam).fill(
                                Math.toDegrees(ele.theta()), beam.e() - pred_e_beam)
                        histos.computeIfAbsent('theta_ele_de_beam_from_angles_' + sector, histoBuilders2.theta_ele_de_beam).fill(
                                Math.toDegrees(ele.theta()), beam.e() - pred_e_beam_from_angles)
                        histos.computeIfAbsent('de_beam_de_beam_from_angles' + sector, histoBuilders2.de_beam_de_beam).fill(
                                beam.e() - pred_e_beam, beam.e() - pred_e_beam_from_angles)


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
out.writeFile("histos.hipo")
