#Equivalent SDOF of archetype M31y: nonlinear time history analysis under the floor motions
# N, mm, s
import openseespy.opensees as op
import os
import sys
import numpy as np
from sdof_from_2d import load_sdof_params
from sdof_nltha_results import save_peak_displacements

# Floor motions are read from the repo's motions/ folder (motions/motion_sets.yaml)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, REPO_DIR)
from piperom.motions import motion_set


# Output in Pushover_SDOF, wherever the script is run from
direc = os.path.join(SCRIPT_DIR, 'ResultsM31y_SDOF')

if(os.path.isdir(direc)==False):
    os.mkdir(direc)

# Peak absolute displacements of all IMs and records of this archetype
peakFile = os.path.join(direc, 'M31y_NLTHA_peak_displacements.csv')

# Envelope recorder file of the record being analysed (read, then deleted)
envFile = os.path.join(direc, 'DispEnvelope_tmp.out')

motions = motion_set('S4_IM')   # 44 records (22 bidirectional pairs), IM levels 1-12
IMs = range(1, 11)              # intensity levels analysed
floor = 4                       # floor where the system is installed

# Equivalent SDOF from the 2D pushover (Pushover2D) at the step closest to dc_target
sdof = load_sdof_params('M31y', dc_target=12.0)

# Viscous damping: mass-proportional (a constant dashpot), XI_MASS of critical at the initial
# (elastic) frequency of the SDOF. The trapezes dissipate energy through their Pinching4 hysteresis.
XI_MASS = 0.01
Gamma = sdof['Gamma']

# Restraint failure and collapse as in the 3D model (3D_models/*_biron.py): each spring is wrapped in
# a MinMax at its 3D cap (mm, at the spring's DOF), i.e. cap/(Gamma*phi) in SDOF displacement. The
# analysis stops when the largest braced-node displacement, Gamma*max(phi)*u, exceeds COLLAPSE_DISP.
CAP_LONG = 60.0       # longitudinal restraints (orthogonals)
CAP_TRANS = 35.0      # transverse restraints (stiff hangers)
COLLAPSE_DISP = 60.0
CHECK_STEPS = 100     # time steps between checks
braced_factor = Gamma*max(phi for phi, _ in sdof['trans'] + sdof['long'])
mass = [sdof['mass'], sdof['mass'], 0]


#SteelTnd Stiff PR
#parameters for C-TPS-L
LePf1 = 600 #floating point values defining force points on the positive response envelope
LePf2  = 7500.0 #floating point values defining force points on the positive response envelope
LePf3 = 10000.0 #floating point values defining force points on the positive response envelope
LePf4 = 11500.0 #floating point values defining force points on the positive response envelope
LePd1 = 0.1 #floating point values defining deformation points on the positive response envelope
LePd2 = 12.0 #floating point values defining deformation points on the positive response envelope
LePd3 = 24.0 #floating point values defining deformation points on the positive response envelope
LePd4 = 61.0 #floating point values defining deformation points on the positive response envelope
LeNf1 = -600.0 #floating point values defining force points on the negative response envelope
LeNf2 = -7500.0 #floating point values defining force points on the negative response envelope
LeNf3 = -10000.0 #floating point values defining force points on the negative response envelope
LeNf4 = -11500.0 #floating point values defining force points on the negative response envelope
LeNd1 = -0.1 #floating point values defining deformation points on the negative response envelope
LeNd2 = -12 #floating point values defining deformation points on the negative response envelope
LeNd3 = -24 #floating point values defining deformation points on the negative response envelope
LeNd4 = -61.0 #floating point values defining deformation points on the negative response envelope
LrDispP = 0.1 #floating point value defining the ratio of the deformationTt which reloading occurs to the maximum historic deformation demand
LrForceP = 0.45 #floating point value defining the ratio of the forceTt which reloading begins to force corresponding to the maximum historic deformation demand
LuForceP = -0.4 #floating point value defining the ratio of strength developed upon unloading from negative load to the maximum strength developed under monotonic loading
LrDispN = 0.1 #floating point value defining the ratio of the deformationTt which reloading occurs to the minimum historic deformation demand
LrForceN = 0.45 #floating point value defining the ratio of the forceTt which reloading begins to force corresponding to the minimum historic deformation demand
LuForceN = -0.4 #floating point value defining the ratio of strength developed upon unloading from negative load to the minimum strength developed under monotonic loading
LgK1 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
LgK2 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
LgK3 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
LgK4 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
LgKLim = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
LgD1 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
LgD2 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
LgD3 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
LgD4 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
LgDLim = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
LgF1 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
LgF2 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
LgF3 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
LgF4 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
LgFLim = 0.0 #floating point values controlling cyclic degradation model for strength degradation
LgE = 10.0 #floating point value used to define maximum energy dissipation under cyclic loading. Total energy dissipation capacity is definedTs this factor multiplied by the energy dissipated under monotonic loading.
LdmgType = "cycle" #string to indicate type of damage (option: "cycle", "energy"

#parameters for C-TPS-T
TePf1 = 600 #floating point values defining force points on the positive response envelope
TePf2 = 6000 #floating point values defining force points on the positive response envelope
TePf3 = 9000 #floating point values defining force points on the positive response envelope
TePf4 = 9100 #floating point values defining force points on the positive response envelope
TePd1 = 0.1 #floating point values defining deformation points on the positive response envelope
TePd2 = 10.0 #floating point values defining deformation points on the positive response envelope
TePd3 = 17.0 #floating point values defining deformation points on the positive response envelope
TePd4 = 36.0 #floating point values defining deformation points on the positive response envelope
TeNf1 = -600 #floating point values defining force points on the negative response envelope
TeNf2 = -6000 #floating point values defining force points on the negative response envelope
TeNf3 = -9000 #floating point values defining force points on the negative response envelope
TeNf4 = -9100 #floating point values defining force points on the negative response envelope
TeNd1 = -0.1 #floating point values defining deformation points on the negative response envelope
TeNd2 = -10.0 #floating point values defining deformation points on the negative response envelope
TeNd3 = -17.0 #floating point values defining deformation points on the negative response envelope
TeNd4 = -36.0 #floating point values defining deformation points on the negative response envelope
TrDispP = 0.1 #floating point value defining the ratio of the deformationTt which reloading occurs to the maximum historic deformation demand
TrForceP = 0.4 #floating point value defining the ratio of the forceTt which reloading begins to force corresponding to the maximum historic deformation demand
TuForceP = -0.3 #floating point value defining the ratio of strength developed upon unloading from negative load to the maximum strength developed under monotonic loading
TrDispN = 0.1 #floating point value defining the ratio of the deformationTt which reloading occurs to the minimum historic deformation demand
TrForceN = 0.4 #floating point value defining the ratio of the forceTt which reloading begins to force corresponding to the minimum historic deformation demand
TuForceN = -0.3 #floating point value defining the ratio of strength developed upon unloading from negative load to the minimum strength developed under monotonic loading
TgK1 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
TgK2 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
TgK3 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
TgK4 = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
TgKLim = 0.0 #floating point values controlling cyclic degradation model for unloading stiffness degradation
TgD1 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
TgD2 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
TgD3 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
TgD4 = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
TgDLim = 0.0 #floating point values controlling cyclic degradation model for reloading stiffness degradation
TgF1 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
TgF2 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
TgF3 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
TgF4 = 0.0 #floating point values controlling cyclic degradation model for strength degradation
TgFLim = 0.0 #floating point values controlling cyclic degradation model for strength degradation
TgE = 10.0 #floating point value used to define maximum energy dissipation under cyclic loading. Total energy dissipation capacity is definedTs this factor multiplied by the energy dissipated under monotonic loading.
TdmgType = "cycle" #string to indicate type of damage (option: "cycle", "energy"


def build_sdof():
    """SDOF model: one zeroLength element with one Pinching4 spring per support of the 2D model."""
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3)

    op.node(nodefix, 0.0, 0.0)
    op.node(nodefree, 0.0, 0.0, '-mass', *mass)

    op.fix(nodefix, 1, 1, 1)

    # Forces scaled by the number of identical springs it represents (n),
    # deformations by 1/(Gamma*phi)
    matTags = []

    for i, (phi, n) in enumerate(sdof['long']):
        matTag = 100 + i
        capTag = 2000 + matTag
        op.uniaxialMaterial('Pinching4', matTag, n*LePf1, LePd1/(Gamma*phi), n*LePf2, LePd2/(Gamma*phi), n*LePf3, LePd3/(Gamma*phi), n*LePf4, LePd4/(Gamma*phi),
                            n*LeNf1, LeNd1/(Gamma*phi), n*LeNf2, LeNd2/(Gamma*phi), n*LeNf3, LeNd3/(Gamma*phi), n*LeNf4, LeNd4/(Gamma*phi),
                            LrDispP, LrForceP, LuForceP, LrDispN, LrForceN, LuForceN, LgK1, LgK2, LgK3, LgK4, LgKLim, LgD1, LgD2, LgD3, LgD4, LgDLim, LgF1, LgF2, LgF3, LgF4, LgFLim, LgE, LdmgType)
        op.uniaxialMaterial('MinMax', capTag, matTag, '-min', -CAP_LONG/(Gamma*phi), '-max', CAP_LONG/(Gamma*phi))
        matTags.append(capTag)

    for i, (phi, n) in enumerate(sdof['trans']):
        matTag = 300 + i
        capTag = 2000 + matTag
        op.uniaxialMaterial('Pinching4', matTag, n*TePf1, TePd1/(Gamma*phi), n*TePf2, TePd2/(Gamma*phi), n*TePf3, TePd3/(Gamma*phi), n*TePf4, TePd4/(Gamma*phi),
                            n*TeNf1, TeNd1/(Gamma*phi), n*TeNf2, TeNd2/(Gamma*phi), n*TeNf3, TeNd3/(Gamma*phi), n*TeNf4, TeNd4/(Gamma*phi),
                            TrDispP, TrForceP, TuForceP, TrDispN, TrForceN, TuForceN, TgK1, TgK2, TgK3, TgK4, TgKLim, TgD1, TgD2, TgD3, TgD4, TgDLim, TgF1, TgF2, TgF3, TgF4, TgFLim, TgE, TdmgType)
        op.uniaxialMaterial('MinMax', capTag, matTag, '-min', -CAP_TRANS/(Gamma*phi), '-max', CAP_TRANS/(Gamma*phi))
        matTags.append(capTag)

    matTot = 1000
    op.uniaxialMaterial('Parallel', matTot, *matTags)

    matRig = 4
    op.uniaxialMaterial('Elastic', 	 matRig,	10e12)

    eleID1 = 1
    op.element('zeroLength', eleID1, nodefix, nodefree, '-mat', matTot, matRig, matRig, '-dir', 1, 2, 3)


nodefix = 1
nodefree = 2

print('Model generation started...')

for im in IMs:
    print('IM = '+str(im)+':')
    peakDisp = []      # peak absolute SDOF displacement (mm) of each record
    completed = []     # True if the analysis reached the end of the record
    collapsedAll = []  # True if a braced node exceeded COLLAPSE_DISP (analysis stopped)
    for Name in motions.records:

        build_sdof()

        # 10th Step: Eigenvalue Analysis
        omega = []
        freq =  []
        T = []

        lamb = op.eigen(1)

        for lam in lamb:
            omega.append((lam)**0.5)
            freq.append((lam)**0.5/(2*np.pi))
            T.append((2*np.pi)/(lam)**0.5)

        for t in range(len(T)):
            print('T'+str(t+1)+' = '+str(T[t])+' s')


        #11th step: Damping: mass-proportional, XI_MASS of critical at the initial frequency of the SDOF
        alphaM = 2*XI_MASS*omega[0]
        op.rayleigh(alphaM, 0.0, 0.0, 0.0)


        #12th step: Run Dynamic Nonlinear time history analysis

        print('Running NLTH Analysis, GM '+Name+'...')

        # Floor acceleration (mm/s^2) at the floor of the system
        motion = motions.load(Name, level=im, floor=floor)
        acc = motion.acc
        dts = motion.dt

        GMfact = 1

        Npts = len(acc)

        dta = dts
        TmaxAnalysis = dts*Npts

        GMfatt =  GMfact

        IDLoadTag = 100
        IDTimeSeries = 1000


        #Define time series
        op.timeSeries('Path',IDTimeSeries,'-dt',dts,'-values',*acc,'-factor',GMfatt, '-prependZero')

        #Define load pattern
        op.pattern('UniformExcitation',IDLoadTag,1,'-accel',IDTimeSeries)

        # Displacement envelope only (rows: min, max, abs max)
        op.recorder('EnvelopeNode', '-file', envFile, '-node', nodefree, '-dof', 1, 'disp')


        testParams = [1.e-4, 500]                        # convergence tolerance for test
        op.wipeAnalysis()
        op.integrator('Newmark', 0.5, 0.25) # determine the next time step for an analysis
        op.numberer('RCM')                  # renumber dof's to minimize band-width (optimization), if you want to
        op.system('FullGeneral')          # how to store and solve the system of equations in the analysis
        op.constraints('Transformation')             # how it handles boundary conditions
        op.test('EnergyIncr',*testParams)      # determine if convergence has been achieved at the end of an iteration step
        op.algorithm('Newton')              # use Newton
        op.analysis('Transient') # define type of analysis static or transient

        # Analysis in chunks of CHECK_STEPS steps, stopped if a braced node exceeds COLLAPSE_DISP
        def over_cap():
            return braced_factor*abs(op.nodeDisp(nodefree, 1)) > COLLAPSE_DISP
        collapsed = False
        done, ok = 0, 0
        while done < Npts and ok == 0:
            nChunk = min(CHECK_STEPS, Npts - done)
            ok = op.analyze(nChunk, dta)     # returns ok=0 if analysis was successful
            done += nChunk
            if ok == 0 and over_cap():
                collapsed = True
                break

        if(ok != 0 and not collapsed):  # analysis was not successful.
            # --------------------------------------------------------------------------------------------------
            # change some analysis parameters to achieve convergence
            # performance is slower inside this loop
            #    Time-controlled analysis
            ok = 0
            controlTime = op.getTime()
            while(controlTime < TmaxAnalysis and ok == 0):
                controlTime = op.getTime()
                if over_cap():
                    collapsed = True
                    break
                ok = op.analyze(1, dta)
                if(ok != 0):
                    testParams2 = [1.e-3, 3000]
                    op.test('EnergyIncr',*testParams2)
                    op.algorithm('Newton','-initial')
                    print("Trying Newton with Initial Tangent ..")
                    ok = op.analyze(1,dta/2)
                    op.test('EnergyIncr',*testParams)
                    op.algorithm('Newton')
                if(ok != 0):
                    op.algorithm('Broyden',50)
                    op.test('EnergyIncr',*testParams2)
                    print("Trying Broyden ..")
                    ok = op.analyze(1,dta/2)
                    op.test('EnergyIncr',*testParams)
                    op.algorithm('Newton')
                if(ok != 0):
                    op.algorithm('NewtonLineSearch')
                    op.test('EnergyIncr',*testParams2)
                    print("Trying NewtonWithLineSearch ..")
                    ok = op.analyze(1,dta/10)
                    op.test('EnergyIncr',*testParams)
                    op.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
                if(ok != 0):
                    print("Trying KrylovNewton ..")
                    op.algorithm('KrylovNewton')
                    op.test('EnergyIncr',*testParams2)
                    ok = op.analyze(1,dta/20)
                    op.test('EnergyIncr',*testParams)
                    op.algorithm('Newton')

        collapsed = collapsed or over_cap()     # also after the last step
        if collapsed:
            print('Collapse: a braced node exceeded %g mm at t = %g s' % (COLLAPSE_DISP, op.getTime()))
        endTime = op.getTime()
        print("Ground Motion Done. End Time:"+str(endTime))
        op.wipe()      # closes the recorder → envelope written

        peakDisp.append(float(np.loadtxt(envFile)[2]))
        completed.append(endTime >= TmaxAnalysis - dta/2)
        collapsedAll.append(collapsed)
        os.remove(envFile)

    # All records of this IM to the archetype's file (written after each IM)
    save_peak_displacements(peakFile, im, motions.records, peakDisp, completed, collapsedAll)

# No recorder files left once all the peaks are in the results file
if os.path.exists(envFile):
    os.remove(envFile)
print('All analyses done. Peak displacements in '+peakFile)
