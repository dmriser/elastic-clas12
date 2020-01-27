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
    high_w: [1.15, 999.9],
    w_loose: [0.8, 1.30],
    angle: [178, 180],
    missing_pt: [0.0, 0.2],
    theta_gamma: [0, 3],
    p_ele:[1.5, 10.646],
    missing_mass:[-0.4, 0.4],
    missing_mass_ftof:[-0.1,0.1]
]

tighter_kin_bounds = [
        theta_ele   : [5, 45],
        theta_pro   : [5, 90],
    theta_sum : [0, 120],
        p_ele       : [0.1, 10.5],
        p_pro       : [0.1, 5.5],
        w           : [0.6, 4.7],
        x           : [0.0, 1.0],
        phi         : [-30, 330],
        dp_ele      : [-3, 3],
        dp_pro      : [-3, 3],
        dtheta_ele  : [-180, 180],
        dtheta_pro  : [-6, 6],
        angle_ep    : [120, 180],
        q2          : [1.2, 4.5],
        missing_pt  : [0, 1],
        e_gamma     : [0, 11],
        theta_gamma :[0, 35],
        theta_egamma:[0, 35],
    missing_mass: [-1, 1],
    chi2:[0, 10],
    dc1:[-150,150],
    dc2:[-250,250],
    dc3:[-350,350]
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
    theta_gamma: { title -> limited_h1(title, 200, lim.theta_gamma) },
    missing_mass: { title -> limited_h1(title, 200, lim.missing_mass) },
    chi2 : { title -> limited_h1(title, 200, lim.chi2) }
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
    w_theta_sum : { title -> limited_h2(title, 200, 200, lim.w, lim.theta_sum) },
    dc1 : {title -> limited_h2(title, 200, 200, lim.dc1, lim.dc1)},
    dc2 : {title -> limited_h2(title, 200, 200, lim.dc2, lim.dc2)},
    dc3 : {title -> limited_h2(title, 200, 200, lim.dc3, lim.dc3)},
    w_missing_mass : { title -> limited_h2(title, 200, 200, lim.w, lim.missing_mass) },
    w_chi2 : {title -> limited_h2(title, 200, 200, lim.w, lim.chi2) },
    w_angle_ep : {title -> limited_h2(title, 200, 200, lim.w, lim.angle_ep) },
    w_theta_gamma : {title -> limited_h2(title, 200, 200, lim.w, lim.theta_gamma) },
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
              
		def hits1 = event.dc1.get(idx)
		def hits2 = event.dc2.get(idx)
		def hits3 = event.dc3.get(idx)
		def hit1 = hits1.find{it.layer==12}
		def hit2 = hits2.find{it.layer==24}
		def hit3 = hits3.find{it.layer==36}

                histos.computeIfAbsent('w_inclusive_', histoBuilders.w).fill(kin.w)

                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    pro.pindex = it
                    pro.sphi = shiftPhi(Math.toDegrees(pro.phi()))
                    def pkin = getPKin(beam, target, ele, pro)

		    def ctof = event.ctof_status.contains(it).findResult{ stat -> stat ? "CTOF" : "FTOF"}

                    histos.computeIfAbsent('w_' + ctof, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('angle_ep_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
                    histos.computeIfAbsent('theta_gamma_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
		    histos.computeIfAbsent('missing_mass_' + ctof, histoBuilders.missing_mass).fill(pkin.missing_mass)

		    def pass_theta_gamma = pkin.theta_gamma < cuts.theta_gamma[1]
		    def pass_angle_ep = pkin.angle > cuts.angle[0] && pkin.angle < cuts.angle[1]
		    def pass_w_elastic = pkin.w < cuts.w[1]
		    def pass_missing_mass = pkin.missing_mass > cuts.missing_mass[0] && pkin.missing_mass < cuts.missing_mass[1]
		    def pass_high_w = pkin.w > cuts.high_w[0]
		    def pass_missing_mass_ftof = pkin.missing_mass > cuts.missing_mass_ftof[0] && pkin.missing_mass < cuts.missing_mass_ftof[1]

		    histos.computeIfAbsent('w_angle_ep_' + ctof, histoBuilders2.w_angle_ep).fill(pkin.w, pkin.angle)

		    if (pass_w_elastic){
			histos.computeIfAbsent('angle_ep_pass_w_elastic_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if (pass_angle_ep){
			histos.computeIfAbsent('w_pass_angle_' + ctof, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('theta_gamma_pass_angle_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
 			histos.computeIfAbsent('theta_egamma_pass_angle_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_egamma)
			histos.computeIfAbsent('theta_e_theta_gamma_pass_angle_' + ctof, histoBuilders2.theta_theta).fill(
			    Math.toDegrees(ele.theta()), pkin.theta_gamma
			)
			histos.computeIfAbsent('missing_mass_pass_angle_' + ctof, 
					       histoBuilders.missing_mass).fill(pkin.missing_mass)

			if (pkin.w > 2){
			    histos.computeIfAbsent('missing_mass_high_w_pass_angle_' + ctof, 
						   histoBuilders.missing_mass).fill(pkin.missing_mass)
			}

			histos.computeIfAbsent('w_theta_sum_pass_angle_' + ctof, histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			)
			histos.computeIfAbsent('w_missing_mass_pass_angle_ep_elastic_' + ctof, histoBuilders2.w_missing_mass).fill(
			    pkin.w, pkin.missing_mass
			)
			histos.computeIfAbsent('chi2_pass_angle_' + ctof, histoBuilders.chi2).fill(event.chi2pid[it])
			histos.computeIfAbsent('w_chi2_pass_angle_' + ctof, histoBuilders2.w_chi2).fill(pkin.w, event.chi2pid[it])
		    }

		    if (pass_theta_gamma){
			histos.computeIfAbsent('w_pass_theta_gamma_' + ctof, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('angle_ep_pass_theta_gamma_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if (pass_angle_ep && pass_w_elastic){
			histos.computeIfAbsent('p_ele_theta_ele_elastic_' + ctof, histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta()))
			if (hit1 && hit2 && hit3){
			    histos.computeIfAbsent('dc1_elastic_' + ctof, histoBuilders2.dc1).fill(hit1.x, hit1.y)
			    histos.computeIfAbsent('dc2_elastic_' + ctof, histoBuilders2.dc2).fill(hit2.x, hit2.y)
			    histos.computeIfAbsent('dc3_elastic_' + ctof, histoBuilders2.dc3).fill(hit3.x, hit3.y)
			}
		    }

		    if (pass_theta_gamma && pass_angle_ep && pass_missing_mass_ftof && pass_high_w){
			histos.computeIfAbsent('w_pass_all_angles_' + ctof, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('p_ele_theta_ele_isr_' + ctof, histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta()))
			histos.computeIfAbsent('w_theta_sum_isr_' + ctof, histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			)
			if (hit1 && hit2 && hit3){
			    histos.computeIfAbsent('dc1_isr_' + ctof, histoBuilders2.dc1).fill(hit1.x, hit1.y)
			    histos.computeIfAbsent('dc2_isr_' + ctof, histoBuilders2.dc2).fill(hit2.x, hit2.y)
			    histos.computeIfAbsent('dc3_isr_' + ctof, histoBuilders2.dc3).fill(hit3.x, hit3.y)
			}

		    }

		    if (pass_missing_mass_ftof){			
			histos.computeIfAbsent('w_theta_sum_pass_missing_mass_' + ctof, histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			)
			if (pass_angle_ep){
			    histos.computeIfAbsent('w_theta_sum_pass_missing_mass_angles_' + ctof, histoBuilders2.w_theta_sum).fill(
				pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			    )
			    histos.computeIfAbsent('w_theta_gamma_' + ctof, histoBuilders2.w_theta_gamma).fill(pkin.w,pkin.theta_gamma)
			}
		    }
 
		    // Everything Else Passed (eep)
		    if (pass_angle_ep && pass_missing_mass_ftof && pass_theta_gamma){
			histos.computeIfAbsent('w_eep_' + ctof, histoBuilders.w).fill(pkin.w)
		    }

		    if (pass_missing_mass_ftof && pass_high_w && pass_theta_gamma){
			histos.computeIfAbsent('angle_ep_eep_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if (pass_missing_mass_ftof && pass_angle_ep && pass_high_w){
			histos.computeIfAbsent('theta_gamma_eep_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
		    }
		    
		    if (pass_high_w && pass_theta_gamma && pass_angle_ep){
			histos.computeIfAbsent('missing_mass_eep_' + ctof, histoBuilders.missing_mass).fill(pkin.missing_mass)
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
out.writeFile("event-selection.hipo")
