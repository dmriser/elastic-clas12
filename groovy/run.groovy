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
    w: [0.8, 1.35],
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
        w           : [0.7, cuts.w[1]],
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
    missing_mass: [-0.25, 0.25],
    chi2:[0, 10],
    dc1:[-150,150],
    dc2:[-250,250],
    dc3:[-350,350],
    cots:[0,20]
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
    w: { title -> limited_h1(title, 200, lim.w) },
    missing_mass: { title -> limited_h1(title, 200, lim.missing_mass) },
]

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

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

	def run = filename.split("/")[-1].split('\\.')[0][-4..-1].toInteger()

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
              
                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    pro.pindex = it
                    def pkin = getPKin(beam, target, ele, pro)

		    def ctof = event.ctof_status.contains(it).findResult{ stat -> stat ? "CTOF" : "FTOF"}
		    
		    def pass_angle = pkin.angle < cuts.angle[1] && pkin.angle > cuts.angle[0]
		    def pass_w = pkin.w < cuts.w[1]

		    if (pass_angle && pkin.missing_mass.abs() < 0.05){
			histos.computeIfAbsent("missing_mass_" + ctof + '_' + run, histoBuilders.missing_mass).fill(pkin.missing_mass)
			histos.computeIfAbsent("w_" + ctof  + '_' + run, histoBuilders.w).fill(pkin.w)
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
out.writeFile("out.hipo")
