
"""
netParams.py 

High-level specifications for M1 network model using NetPyNE

Contributors: salvadordura@gmail.com
"""

from netpyne import specs
import pickle, json

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

netParams.version = 56

from cfg_cell import cfg

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------
netParams.scale = cfg.scale # Scale factor for number of cells
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0 # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightModels = {'HH_simple': 1.0, 'HH_reduced': 1.0, 'HH_full': 1.0} #scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0 #0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
netParams.defineCellShapes = True  # convert stylized geoms to 3d points


# special condition to change Kgbar together with ih when running batch
# note min Kgbar is assumed to be 0.5, so this is set here as an offset 
if cfg.makeKgbarFactorEqualToNewFactor:
    cfg.KgbarFactor = 0.5 + cfg.modifyMechs['newFactor']

#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------
cellModels = ['HH_simple', 'HH_reduced', 'HH_full']
layer = {'1':[0.0, 0.1], '2': [0.1,0.29], '4': [0.29,0.37], '5A': [0.37,0.47], '24':[0.1,0.37], '5B': [0.47,0.8], '6': [0.8,1.0], 
'longTPO': [2.0,2.1], 'longTVL': [2.1,2.2], 'longS1': [2.2,2.3], 'longS2': [2.3,2.4], 'longcM1': [2.4,2.5], 'longM2': [2.5,2.6], 'longOC': [2.6,2.7]}  # normalized layer boundaries

netParams.correctBorder = False
    #{'threshold': [cfg.correctBorderThreshold, cfg.correctBorderThreshold, cfg.correctBorderThreshold],
                        #'yborders': [layer['2'][0], layer['5A'][0], layer['6'][0], layer['6'][1]]}  # correct conn border effect

#------------------------------------------------------------------------------
## Load cell rules previously saved using netpyne format
cellParamLabels = ['PT5B_full'] #'NGF_simple', 'VIP_reduced'] # list of cell rules to load from file
loadCellParams = cellParamLabels
saveCellParams = False #True

for ruleLabel in loadCellParams:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/'+ruleLabel+'_cellParams.pkl')
    
    # Adapt K gbar
    if ruleLabel in ['IT2_reduced', 'IT4_reduced', 'IT5A_reduced', 'IT5B_reduced', 'IT6_reduced', 'CT6_reduced', 'IT5A_full']:
        cellRule = netParams.cellParams[ruleLabel]
        for secName in cellRule['secs']:
            for kmech in [k for k in cellRule['secs'][secName]['mechs'].keys() if k.startswith('k') and k!='kBK']:
                cellRule['secs'][secName]['mechs'][kmech]['gbar'] *= cfg.KgbarFactor 

#------------------------------------------------------------------------------
# Specification of cell rules not previously loaded
# Includes importing from hoc template or python class, and setting additional params

#------------------------------------------------------------------------------
# Reduced cell model params (6-comp) deleted

#------------------------------------------------------------------------------
## PT5B full cell model params (700+ comps)
#UC Davis PT Cell

if 'PT5B_full' not in loadCellParams:
    ihMod2str = {'harnett': 1, 'kole': 2, 'migliore': 3}

    netParams.loadCellParams('PT5B_full', 'Na1216TF.pkl')
    netParams.renameCellParamsSec(label='PT5B_full', oldSec ='soma_0', newSec ='soma')
    cellRule = netParams.cellParams['PT5B_full']

    cellRule['secs']['axon_0']['geom']['pt3d'] = [[1e30, 1e30, 1e30]]
    cellRule['secs']['axon_1']['geom']['pt3d'] = [[1e30, 1e30, 1e30]]

    nonSpiny = ['apic_0', 'apic_1']
    netParams.addCellParamsSecList(label='PT5B_full', secListName='perisom', somaDist=[0, 50])  # sections within 50 um of soma
    netParams.addCellParamsSecList(label='PT5B_full', secListName='below_soma', somaDistY=[-600, 0])  # sections within 0-300 um of soma
    for sec in nonSpiny: # N.B. apic_1 not in `perisom` . `apic_0` and `apic_114` are
        if sec in cellRule['secLists']['perisom']: # fixed logic
            cellRule['secLists']['perisom'].remove(sec)
    cellRule['secLists']['alldend'] = [sec for sec in cellRule.secs if ('dend' in sec or 'apic' in sec)] # basal+apical
    cellRule['secLists']['apicdend'] = [sec for sec in cellRule.secs if ('apic' in sec)] # apical
    cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['alldend'] if sec not in nonSpiny]

#     netParams.addCellParamsWeightNorm('PT5B_full', 'conn/PT5B_full_weightNorm.pkl', threshold=cfg.weightNormThreshold)  # load weight norm
#     if saveCellParams: netParams.saveCellParamsRule(label='PT5B_full', fileName='cells/PT5B_full_cellParams.pkl')


#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
## load densities
with open('cells/cellDensity.pkl', 'rb') as fileObj: density = pickle.load(fileObj)['density']

## Local populations


netParams.popParams['PT5B'] = {'cellModel': cfg.cellmod['PT5B'], 'cellType': 'PT', 'ynormRange': layer['5B'], 'density': 0.5*density[('M1','E')][3]}

if cfg.singleCellPops:
    for pop in netParams.popParams.values(): pop['numCells'] = 1

#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3*cfg.AMPATau2Factor, 'e': 0}
netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93} 
netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
netParams.synMechParams['GABAASlow'] = {'mod': 'MyExp2SynBB','tau1': 2, 'tau2': 100, 'e': -80}
netParams.synMechParams['GABAASlowSlow'] = {'mod': 'MyExp2SynBB', 'tau1': 200, 'tau2': 400, 'e': -80}

ESynMech = ['AMPA', 'NMDA']
SOMESynMech = ['GABAASlow','GABAB']
SOMISynMech = ['GABAASlow']
PVSynMech = ['GABAA']


#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
    for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}

#------------------------------------------------------------------------------
# NetStim inputs
#------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, ynorm, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
        [params[s] for s in ['pop', 'ynorm', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']] 

        # cfg.analysis['plotTraces']['include'] = [(pop,0)]

        if synMech == ESynMech:
            wfrac = cfg.synWeightFractionEE
        elif synMech == SOMESynMech:
            wfrac = cfg.synWeightFractionSOME
        else:
            wfrac = [1.0]

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        # connect stim source to target
        # for i, syn in enumerate(synMech):
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop, 'ynorm': ynorm},
            'sec': sec, 
            'loc': loc,
            'synMech': synMech,
            'weight': weight,
            'synMechWeightFactor': synMechWeightFactor,
            'delay': delay}

#------------------------------------------------------------------------------
# Local connectivity parameters
#------------------------------------------------------------------------------
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
pmat = connData['pmat']
wmat = connData['wmat']
bins = connData['bins']
#import conn_param
# pmat = conn_param.pmat
# wmat = conn_param.wmat
# bins = conn_param.bins

#------------------------------------------------------------------------------
## E -> E
if cfg.addConn and cfg.EEGain > 0.0:
    labelsConns = [('W+AS_norm', 'IT', 'L2/3,4'), ('W+AS_norm', 'IT', 'L5A,5B'), 
                   ('W+AS_norm', 'PT', 'L5B'), ('W+AS_norm', 'IT', 'L6'), ('W+AS_norm', 'CT', 'L6')]
    labelPostBins = [('W+AS', 'IT', 'L2/3,4'), ('W+AS', 'IT', 'L5A,5B'), ('W+AS', 'PT', 'L5B'), 
                    ('W+AS', 'IT', 'L6'), ('W+AS', 'CT', 'L6')]
    labelPreBins = ['W', 'AS', 'AS', 'W', 'W']
    preTypes = [['IT'], ['IT'], ['IT', 'PT'], ['IT','CT'], ['IT','CT']] 
    postTypes = ['IT', 'IT', 'PT', 'IT','CT']
    ESynMech = ['AMPA','NMDA']

    for i,(label, preBinLabel, postBinLabel) in enumerate(zip(labelsConns,labelPreBins, labelPostBins)):
        for ipre, preBin in enumerate(bins[preBinLabel]):
            for ipost, postBin in enumerate(bins[postBinLabel]):
                for cellModel in cellModels:
                    ruleLabel = 'EE_'+cellModel+'_'+str(i)+'_'+str(ipre)+'_'+str(ipost)
                    netParams.connParams[ruleLabel] = { 
                        'preConds': {'cellType': preTypes[i], 'ynorm': list(preBin)}, 
                        'postConds': {'cellModel': cellModel, 'cellType': postTypes[i], 'ynorm': list(postBin)},
                        'synMech': ESynMech,
                        'probability': pmat[label][ipost,ipre],
                        'weight': wmat[label][ipost,ipre] * cfg.EEGain / cfg.synsperconn[cellModel], 
                        'synMechWeightFactor': cfg.synWeightFractionEE,
                        'delay': 'defaultDelay+dist_3D/propVelocity',
                        'synsPerConn': cfg.synsperconn[cellModel],
                        'sec': 'spiny'}
            

#------------------------------------------------------------------------------
## E -> I
if cfg.EIGain: # Use IEGain if value set
    cfg.EPVGain = cfg.EIGain
    cfg.ESOMGain = cfg.EIGain
else: 
    cfg.EIGain = (cfg.EPVGain+cfg.ESOMGain)/2.0

if cfg.addConn and (cfg.EPVGain > 0.0 or cfg.ESOMGain > 0.0):
    labelsConns = ['FS', 'LTS']
    labelPostBins = ['FS/LTS', 'FS/LTS']
    labelPreBins = ['FS/LTS', 'FS/LTS']
    preTypes = ['IT', 'PT', 'CT']
    postTypes = ['PV', 'SOM']
    ESynMech = ['AMPA','NMDA']
    lGain = [cfg.EPVGain, cfg.ESOMGain] # E -> PV or E -> SOM
    for i,(label, preBinLabel, postBinLabel) in enumerate(zip(labelsConns,labelPreBins, labelPostBins)):
        for ipre, preBin in enumerate(bins[preBinLabel]):
            for ipost, postBin in enumerate(bins[postBinLabel]):
                ruleLabel = 'EI_'+str(i)+'_'+str(ipre)+'_'+str(ipost)
                netParams.connParams[ruleLabel] = {
                    'preConds': {'cellType': preTypes, 'ynorm': list(preBin)},
                    'postConds': {'cellType': postTypes[i], 'ynorm': list(postBin)},
                    'synMech': ESynMech,
                    'probability': pmat[label][ipost,ipre],
                    'weight': wmat[label][ipost,ipre] * lGain[i],
                    'synMechWeightFactor': cfg.synWeightFractionEI,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'sec': 'soma'} # simple I cells used right now only have soma



#------------------------------------------------------------------------------
## I -> all
if cfg.IEGain: # Use IEGain if value set
    cfg.PVEGain = cfg.IEGain
    cfg.SOMEGain = cfg.IEGain
else: 
    cfg.IEGain = (cfg.PVEGain+cfg.SOMEGain)/2.0

if cfg.IIGain:  # Use IIGain if value set
    cfg.SOMPVGain = cfg.IIGain
    cfg.PVSOMGain = cfg.IIGain
    cfg.SOMSOMGain = cfg.IIGain
    cfg.PVPVGain = cfg.IIGain
else:
    cfg.IIGain = (cfg.PVSOMGain+cfg.SOMPVGain+cfg.SOMSOMGain+cfg.PVPVGain)/4.0

if cfg.addConn and (cfg.IEGain > 0.0 or cfg.IIGain > 0.0):
    # Local, intralaminar only; all-to-all but distance-based; high weights; L5A/B->L5A/B
    preCellTypes = ['SOM', 'SOM', 'SOM', 'PV', 'PV', 'PV']
    ynorms = [(layer['2'][0],layer['4'][1]), (layer['5A'][0],layer['5B'][1]), (layer['6'][0],layer['6'][1])]*2
    IEweights = cfg.IEweights * 2  # [I->E2/3+4, I->E5, I->E6] weights (Note * 2 is repeat list operator)
    IIweights = cfg.IIweights * 2  # [I->I2/3+4, I->I5, I->I6] weights (Note * 2 is repeat list operator)
    postCellTypes = ['PT', ['IT','CT'], 'PV', 'SOM']
    IEdisynBiases = [None, cfg.IEdisynapticBias, cfg.IEdisynapticBias, None, cfg.IEdisynapticBias, cfg.IEdisynapticBias]
    disynapticBias = None  # default, used for I->I

    # preCellTypes = ['SOM', 'PV']
    # ynorms = [[0,1]]*2
    # IEweights = cfg.IEweights * 2  # [I->E2/3+4, I->E5, I->E6] weights (Note * 2 is repeat list operator)
    # IIweights = cfg.IIweights * 2  # [I->I2/3+4, I->I5, I->I6] weights (Note * 2 is repeat list operator)
    # postCellTypes = ['PT', ['IT','CT'], 'PV', 'SOM']
    # IEdisynBiases = [cfg.IEdisynapticBias, cfg.IEdisynapticBias]
    # disynapticBias = None  # default, used for I->I

    for i,(preCellType, ynorm, IEweight, IIweight, IEdisynBias) in enumerate(zip(preCellTypes, ynorms, IEweights, IIweights, IEdisynBiases)):
        for ipost, postCellType in enumerate(postCellTypes):
            for cellModel in cellModels:
                if postCellType == 'PV':    # postsynaptic I cell
                    sec = 'soma'
                    synWeightFraction = [1]
                    if preCellType == 'PV':             # PV->PV
                        weight = IIweight * cfg.PVPVGain
                        synMech = PVSynMech
                    else:                           # SOM->PV
                        weight = IIweight * cfg.SOMPVGain
                        synMech = SOMISynMech
                elif postCellType == 'SOM': # postsynaptic I cell
                    sec = 'soma'
                    synWeightFraction = [1]
                    if preCellType == 'PV':             # PV->SOM
                        weight = IIweight * cfg.PVSOMGain
                        synMech = PVSynMech
                    else:                           # SOM->SOM
                        weight = IIweight * cfg.SOMSOMGain
                        synMech = SOMISynMech
                elif postCellType == ['IT','CT']: # postsynaptic IT,CT cell
                    disynapticBias = IEdisynBias
                    if preCellType == 'PV':             # PV->E
                        weight = IEweight * cfg.PVEGain
                        synMech = PVSynMech
                        sec = 'perisom'
                    else:                           # SOM->E
                        weight = IEweight * cfg.SOMEGain
                        synMech = SOMESynMech
                        sec = 'spiny'
                        synWeightFraction = cfg.synWeightFractionSOME
                elif postCellType == 'PT': # postsynaptic PT cell
                    disynapticBias = IEdisynBias
                    if preCellType == 'PV':             # PV->E
                        weight = IEweight * cfg.IPTGain * cfg.PVEGain
                        synMech = PVSynMech
                        sec = 'perisom'
                    else:                           # SOM->E
                        weight = IEweight * cfg.IPTGain * cfg.SOMEGain
                        synMech = SOMESynMech
                        sec = 'spiny'
                        synWeightFraction = cfg.synWeightFractionSOME
                if cellModel == 'HH_full':
                    weight = weight * cfg.IFullGain


                ruleLabel = 'I_'+cellModel+'_'+str(i)+'_'+str(ipost)
                netParams.connParams[ruleLabel] = {
                    'preConds': {'cellType': preCellType, 'ynorm': ynorm},
                    'postConds': {'cellModel': cellModel, 'cellType': postCellType, 'ynorm': ynorm},
                    'synMech': synMech,
                    'probability': '1.0 * exp(-dist_3D_border/probLambda)',
                    'weight': weight / cfg.synsperconn[cellModel],
                    'delay': 'defaultDelay+dist_3D_border/propVelocity',
                    'synsPerConn': cfg.synsperconn[cellModel],
                    'synMechWeightFactor': synWeightFraction,
                    'sec': sec,
                    'disynapticBias': disynapticBias}


#------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
#------------------------------------------------------------------------------         
if cfg.addSubConn:
    with open('conn/conn_dend_PT.json', 'r') as fileObj: connDendPTData = json.load(fileObj)
    with open('conn/conn_dend_IT.json', 'r') as fileObj: connDendITData = json.load(fileObj)
    
    #------------------------------------------------------------------------------
    # L2/3,TVL,S2,cM1,M2 -> PT (Suter, 2015)
    lenY = 30 
    spacing = 50
    gridY = range(0, -spacing*lenY, -spacing)
    synDens, _, fixedSomaY = connDendPTData['synDens'], connDendPTData['gridY'], connDendPTData['fixedSomaY']
    for k in synDens.keys():
        prePop,postType = k.split('_')  # eg. split 'M2_PT'
        if prePop == 'L2': prePop = 'IT2'  # include conns from layer 2/3 and 4
        netParams.subConnParams[k] = {
        'preConds': {'pop': prePop}, 
        'postConds': {'cellType': postType},  
        'sec': 'spiny',
        'groupSynMechs': ESynMech, 
        'density': {'type': '1Dmap', 'gridX': None, 'gridY': gridY, 'gridValues': synDens[k], 'fixedSomaY': fixedSomaY}} 


    #------------------------------------------------------------------------------
    # TPO, TVL, M2, OC  -> E (L2/3, L5A, L5B, L6) (Hooks 2013)
    lenY = 26
    spacing = 50
    gridY = range(0, -spacing*lenY, -spacing)
    synDens, _, fixedSomaY = connDendITData['synDens'], connDendITData['gridY'], connDendITData['fixedSomaY']
    for k in synDens.keys():
        prePop,post = k.split('_')  # eg. split 'M2_L2'
        postCellTypes = ['IT','PT','CT'] if prePop in ['OC','TPO'] else ['IT','CT']  # only OC,TPO include PT cells
        postyRange = list(layer[post.split('L')[1]]) # get layer yfrac range 
        if post == 'L2': postyRange[1] = layer['4'][1]  # apply L2 rule also to L4 
        netParams.subConnParams[k] = {
        'preConds': {'pop': prePop}, 
        'postConds': {'ynorm': postyRange , 'cellType': postCellTypes},  
        'sec': 'spiny',
        'groupSynMechs': ESynMech, 
        'density': {'type': '1Dmap', 'gridX': None, 'gridY': gridY, 'gridValues': synDens[k], 'fixedSomaY': fixedSomaY}} 


    #------------------------------------------------------------------------------
    # S1, S2, cM1 -> E IT/CT; no data, assume uniform over spiny
    netParams.subConnParams['S1,S2,cM1->IT,CT'] = {
        'preConds': {'pop': ['S1','S2','cM1']}, 
        'postConds': {'cellType': ['IT','CT']},
        'sec': 'spiny',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'} 


    #------------------------------------------------------------------------------
    # rest of local E->E (exclude IT2->PT); uniform distribution over spiny
    netParams.subConnParams['IT2->non-PT'] = {
        'preConds': {'pop': ['IT2']}, 
        'postConds': {'cellType': ['IT','CT']},
        'sec': 'spiny',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'} 
        
    netParams.subConnParams['non-IT2->E'] = {
        'preConds': {'pop': ['IT4','IT5A','IT5B','PT5B','IT6','CT6']}, 
        'postConds': {'cellType': ['IT','PT','CT']},
        'sec': 'spiny',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'} 


    #------------------------------------------------------------------------------
    # PV->E; perisomatic (no sCRACM)
    netParams.subConnParams['PV->E'] = {
        'preConds': {'cellType': 'PV'}, 
        'postConds': {'cellType': ['IT', 'CT', 'PT']},  
        'sec': 'perisom', 
        'density': 'uniform'} 


    #------------------------------------------------------------------------------
    # SOM->E; apical dendrites (no sCRACM)
    netParams.subConnParams['SOM->E'] = {
        'preConds': {'cellType': 'SOM'}, 
        'postConds': {'cellType': ['IT', 'CT', 'PT']},  
        'sec': 'apicdend',
        'groupSynMechs': SOMESynMech,
        'density': 'uniform'} 


    #------------------------------------------------------------------------------
    # All->I; apical dendrites (no sCRACM)
    netParams.subConnParams['All->I'] = {
        'preConds': {'cellType': ['IT', 'CT', 'PT', 'SOM', 'PV']},
        'postConds': {'cellType': ['SOM', 'PV']},  
        'sec': 'spiny',
        'groupSynMechs': ESynMech,
        'density': 'uniform'} 


#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
netParams.description = """ 
- M1 net, 6 layers, 7 cell types 
- NCD-based connectivity from  Weiler et al. 2008; Anderson et al. 2010; Kiritani et al. 2012; 
  Yamawaki & Shepherd 2015; Apicella et al. 2012
- Parametrized version based on Sam's code
- Updated cell models and mod files
- Added parametrized current inputs
- Fixed bug: prev was using cell models in /usr/site/nrniv/local/python/ instead of cells 
- Use 5 synsperconn for 5-comp cells (HH_reduced); and 1 for 1-comp cells (HH_simple)
- Fixed bug: made global h params separate for each cell model
- Fixed v_init for different cell models
- New IT cell with same geom as PT
- Cleaned cfg and moved background inputs here
- Set EIGain and IEGain for each inh cell type
- Added secLists for PT full
- Fixed reduced CT (wrong vinit and file)
- Added subcellular conn rules to distribute synapses
- PT full model soma centered at 0,0,0 
- Set cfg seeds here to ensure they get updated
- Added PVSOMGain and SOMPVGain
- PT subcellular distribution as a cfg param
- Cylindrical volume
- DefaultDelay (for local conns) = 2ms
- Added long range connections based on Yamawaki 2015a,b; Suter 2015; Hooks 2013; Meyer 2011
- Updated cell densities based on Tsai 2009; Lefort 2009; Katz 2011; Wall 2016; 
- Separated PV and SOM of L5A vs L5B
- Fixed bugs in local conn (PT, PV5, SOM5, L6)
- Added perisom secList including all sections 50um from soma
- Added subcellular conn rules (for both full and reduced models)
- Improved cell models, including PV and SOM fI curves
- Improved subcell conn rules based on data from Suter15, Hooks13 and others
- Adapted Bdend L of reduced cell models
- Made long pop rates a cfg param
- Set threshold to 0.0 mV
- Parametrized I->E/I layer weights
- Added missing subconn rules (IT6->PT; S1,S2,cM1->IT/CT; long->SOM/PV)
- Added threshold to weightNorm (PT threshold=10x)
- weightNorm threshold as a cfg parameter
- Separate PV->SOM, SOM->PV, SOM->SOM, PV->PV gains 
- Conn changes: reduced IT2->IT4, IT5B->CT6, IT5B,6->IT2,4,5A, IT2,4,5A,6->IT5B; increased CT->PV6+SOM6
- Parametrized PT ih gbar
- Added IFullGain parameter: I->E gain for full detailed cell models
- Replace PT ih with Migliore 2012
- Parametrized ihGbar, ihGbarBasal, dendNa, axonNa, axonRa, removeNa
- Replaced cfg list params with dicts
- Parametrized ihLkcBasal and AMPATau2Factor
- Fixed synMechWeightFactor
- Parametrized PT ih slope
- Added disynapticBias to I->E (Yamawaki&Shepherd,2015)
- Fixed E->CT bin 0.9-1.0
- Replaced GABAB with exp2syn and adapted synMech ratios
- Parametrized somaNa
- Added ynorm condition to NetStims
- Added option to play back recorded spikes into long-range inputs
- Fixed Bdend pt3d y location
- Added netParams.convertCellShapes = True to convert stylized geoms to 3d points
- New layer boundaries, cell densities, conn, FS+SOM L4 grouped with L2/3, low cortical input to L4
- Increased exc->L4 based on Yamawaki 2015 fig 5
- v54: Moved from NetPyNE v0.7.9 to v0.9.1 (v54_batch1-6)
- v54: Moved to NetPyNE v0.9.1 and py3 (v54_batch7 onwards)
- v56: Reduced dt from 0.05 to 0.025 (note this version follows from v54, i.e. without new cell types; branch 'paper2019_py3')
- v56: (included in prev version): Added cfg.KgbarFactor
"""
