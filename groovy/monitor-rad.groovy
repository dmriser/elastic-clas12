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
    high_w:[1.15,999.9],
    w_loose: [0.8, 1.30],
    angle: [178, 180],
    missing_pt: [0.0, 0.2],
    theta_gamma: [0, 3],
    p_ele:[1.5, 10.646],
    missing_mass: [-0.1, 0.1]
]

tighter_kin_bounds = [
        theta_ele : [5, 45],
        theta_pro : [5, 90],
        p_ele     : [0.1, 10.5],
        p_pro     : [0.1, 5.5],
        w         : [0.6, 4.7],
        x         : [0.0, 1.0],
        phi       : [-30, 330],
        dp_ele    : [-1, 1],
        dp_pro    : [-1, 1],
        dtheta_ele: [-15, 15],
        dtheta_pro: [-6, 6],
        angle_ep  : [120, 180],
        q2        : [1.2, 4.5],
    missing_pt    : [0, 1],
    e_gamma: [0, 11],
    theta_gamma:[0, 35],
    theta_egamma:[0, 35],
    chi2:[0,10]
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
        p_ele    : { title -> limited_h1(title, 200, lim.p_ele) },
        p_pro    : { title -> limited_h1(title, 200, lim.p_pro) },
        de_beam  : { title -> limited_h1(title, 200, lim.de_beam) },
        angle_ep : { title -> limited_h1(title, 200, lim.angle_ep) },
        theta_p  : { title -> limited_h1(title, 200, lim.theta_pro) },
        theta_ele: { title -> limited_h1(title, 200, lim.theta_ele) },
        theta_pro: { title -> limited_h1(title, 200, lim.theta_pro) },
        e_gamma: { title -> limited_h1(title, 200, lim.e_gamma) },
    theta_gamma: { title -> limited_h1(title, 200, lim.theta_gamma) }
]

histoBuilders2 = [
        p_pro_dp         : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dp_pro) },
        p_ele_dp         : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dp_ele) },
        p_ele_theta      : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.theta_ele) },
        p_pro_theta      : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.theta_pro) },
        theta_ele_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dtheta_ele) },
        theta_pro_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dtheta_pro) },
        w_q2             : { title -> limited_h2(title, 200, 200, lim.w, lim.q2) },
        x_q2             : { title -> limited_h2(title, 200, 200, lim.x, lim.q2) },
    theta_theta: { title -> limited_h2(title, 200, 200, lim.theta_egamma, lim.theta_gamma) },
    chi2_dp_ele      : { title -> limited_h2(title, 100, 100, lim.chi2, lim.dp_ele) },
    chi2_dp_pro      : { title -> limited_h2(title, 100, 100, lim.chi2, lim.dp_pro) },
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
    //def phi = Math.toDegrees(electron.phi() - proton.phi())
    def missing_pt = Math.sqrt(missing.px()**2 + missing.py()**2)
    def theta_gamma = Math.toDegrees(missing.theta())
    def norm = missing.vector().vect().mag() * electron.vector().vect().mag()
    def theta_egamma = Math.toDegrees(Math.acos(missing.vector().vect().dot(electron.vector().vect()) / norm))

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e(),
	    missing_pt: missing_pt, theta_gamma:theta_gamma, 
	    theta_egamma:theta_egamma]
}

def predictElectron(pro){
    //def a = pro.p()**2 
    //def b = -2 * pro.p() * Math.cos(pro.theta())
    //def c = PDGDatabase.getParticleMass(2212) - pro.e()
    //def beamEnergy = (c**2 - a) / (b - 2 * c)

    def mp = PDGDatabase.getParticleMass(2212)
    def beam_num = mp * (pro.p() + (-mp + Math.sqrt(mp**2 + pro.p()**2)) * Math.cos(pro.theta()))
    def beam_den = 2 * mp * Math.cos(pro.theta()) - pro.p() * Math.sin(pro.theta())**2
    def beamEnergy  = beam_num / beam_den

    def den = pro.p() * Math.sin(-1 * pro.theta())
    def num = beamEnergy - pro.p() * Math.cos(pro.theta())
    def alpha = Math.atan2(den,num)

    def kprime = pro.p() * Math.sin(-1 * pro.theta()) / Math.sin(alpha)

    return [momentum:kprime, theta:alpha]
}

def predictProton(ele){
    def beamEnergy = ele.p() / (1 + ele.p() / PDGDatabase.getParticleMass(2212) * (Math.cos(ele.theta()) - 1))
    def beta = Math.atan(ele.p() * Math.sin(ele.theta()) / (beamEnergy - ele.p() * Math.cos(ele.theta())))
    def pprime = ele.p() * Math.sin(ele.theta()) / Math.sin(beta)
    return [momentum:pprime, theta:beta]
}

def predictElectronMathematica(pro){
    // This one actually works! 
    def M = PDGDatabase.getParticleMass(2212)
    def Pp = pro.p()
    def Ep = Math.sqrt(M**2 + Pp**2)
    def beta = pro.theta()
    def cbeta = Math.cos(beta)
    def c2beta = Math.cos(2 * beta)
    def c3beta = Math.cos(3 * beta)
    def sbeta = Math.sin(-beta)
    
    // Angle 
    def expr_num = 2 * M**2 - 2 * M * Ep - Pp * (M + Ep) * cbeta
    expr_num += 2 * M * (M + Ep) * c2beta + M * Pp * c3beta + Pp * Ep * c3beta
    def expr_den = 2 * (2 * M**2 + Pp**2 - Pp**2 * c2beta)
    def alpha = Math.acos(-1 * expr_num/expr_den)
    
    // Mom
    expr_num = -2 * M * (-M + Ep)  * cbeta + Pp * (M + Ep) + (M - Ep) * c2beta
    expr_den = 4 * M * cbeta - 2 * Pp * sbeta**2
    def kprime = expr_num / expr_den

    return [momentum:kprime, theta:alpha]
}

def predictElectronMomentum(ele, pro){
    def M = PDGDatabase.getParticleMass(2212)
    return M * Math.cos(ele.theta()/2 + pro.theta()) / Math.sin(ele.theta()/2) / Math.sin(ele.theta() + pro.theta())
}

def predictProtonMomentum(ele, pro){
    def M = PDGDatabase.getParticleMass(2212)
    return 2 * M * Math.cos(ele.theta()/2) * Math.cos(ele.theta()/2 + pro.theta()) / Math.sin(pro.theta()) / Math.sin(ele.theta() + pro.theta())    
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
                event.pid[it] == 11 && event.status[it] < 0 && event.p[it] > cuts.p_ele[0]
            }?.each { idx ->
                def sector = event.dc_sector[idx]
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                ele.sector = sector
                ele.pindex = idx
                ele.sphi = shiftPhi(Math.toDegrees(ele.phi()))
                def kin = getKin(beam, target, ele)
             
                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    pro.pindex = it
                    pro.sphi = shiftPhi(Math.toDegrees(pro.phi()))
                    def pkin = getPKin(beam, target, ele, pro)

		    def ctof = event.ctof_status.contains(it).findResult{ stat -> stat ? "CTOF" : "FTOF"}

                    histos.computeIfAbsent('w_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('angle_ep_' + ctof + '_' + sector, histoBuilders.angle_ep).fill(pkin.angle)
                    histos.computeIfAbsent('theta_gamma_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_gamma)


		    if (pkin.theta_gamma < cuts.theta_gamma[1]){
			histos.computeIfAbsent('w_pass_theta_gamma_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('angle_ep_pass_theta_gamma_' + ctof + '_' + sector, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if (pkin.angle > cuts.angle[0] && pkin.angle < cuts.angle[1]){
			histos.computeIfAbsent('w_pass_angle_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
 			histos.computeIfAbsent('w_pass_angle_'  + ctof, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('theta_gamma_pass_angle_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
 			histos.computeIfAbsent('theta_egamma_pass_angle_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_egamma)
			histos.computeIfAbsent('theta_e_theta_gamma_pass_angle_' + ctof + '_' + sector, histoBuilders2.theta_theta).fill(
			    Math.toDegrees(ele.theta()), pkin.theta_gamma
			)

			def pass_missing_mass = pkin.missing_mass > cuts.missing_mass[0] && pkin.missing_mass < cuts.missing_mass[1]
			def pass_high_w = pkin.w > cuts.high_w[0]
			def pass_angle_ep = pkin.angle > cuts.angle[0] && pkin.angle < cuts.angle[1]

			// These are ISR events. 
			if (pkin.theta_gamma < cuts.theta_gamma[1] && pass_missing_mass && pass_high_w && pass_angle_ep) {
 			    histos.computeIfAbsent('w_pass_all_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
			    histos.computeIfAbsent('p_ele_theta_ele_' + ctof + '_' + sector, histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('p_pro_theta_pro_' + ctof + '_' + sector, histoBuilders2.p_pro_theta).fill(
			    pro.p(), Math.toDegrees(pro.theta()))

			    histos.computeIfAbsent('p_ele_' + ctof + '_' + sector, histoBuilders.p_ele).fill(ele.p())
			    histos.computeIfAbsent('theta_ele_' + ctof + '_' + sector, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('p_pro_' + ctof + '_' + sector, histoBuilders.p_pro).fill(pro.p())
			    histos.computeIfAbsent('theta_pro_' + ctof + '_' + sector, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))

			    // Resolutions 
			    def pred_ele = predictElectronMathematica(pro)
			    def pred_pro = predictProton(ele)
			    
			    def pred_ele_p = predictElectronMomentum(ele, pro)
			    def pred_pro_p = predictProtonMomentum(ele, pro)

			    histos.computeIfAbsent('p_ele_dp_ele_' + ctof + '_' + sector, histoBuilders2.p_ele_dp).fill(ele.p(), pred_ele.momentum - ele.p())
			    histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector, histoBuilders2.p_ele_dp).fill(
				ele.p(), pred_ele_p - ele.p())
			    histos.computeIfAbsent('p_pro_dp_pro_' + ctof + '_' + sector, histoBuilders2.p_pro_dp).fill(pro.p(), pred_pro.momentum - pro.p())
			    histos.computeIfAbsent('p_pro_dp_pro_from_angles_' + ctof + '_' + sector, histoBuilders2.p_pro_dp).fill(
				pro.p(), pred_pro_p - pro.p())

			    histos.computeIfAbsent('theta_ele_dtheta_ele_' + ctof + '_' + sector, histoBuilders2.theta_ele_dtheta).fill(
				Math.toDegrees(ele.theta()), Math.toDegrees(pred_ele.theta - ele.theta()))
			    histos.computeIfAbsent('theta_pro_dtheta_pro_' + ctof + '_' + sector, histoBuilders2.theta_pro_dtheta).fill(
				Math.toDegrees(pro.theta()), Math.toDegrees(pred_pro.theta - pro.theta()))
			    
			    histos.computeIfAbsent('chi2_dp_ele_' + ctof, histoBuilders2.chi2_dp_ele).fill(event.chi2pid[it], pred_ele_p - ele.p())
			    histos.computeIfAbsent('chi2_dp_pro_' + ctof, histoBuilders2.chi2_dp_ele).fill(event.chi2pid[it], pred_pro_p - pro.p())
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
out.writeFile("mon-rad.hipo")
