import numpy as np
import matplotlib.pyplot as plt
import uproot
import ROOT
import sys
import os
import csv
import sys
import pandas as pd

#######################################################################################################################
#Define Useful Functions
#######################################################################################################################
print('Started script')

if sys.argv[1] == "write":
    write = 1
else:
    write = 0

#Defining fiducial cuts
polyParams = [6.92312336e+02, -2.43739552e+02,  2.84329329e+03, -2.09128802e+04, 9.67156033e+04, -2.84628041e+05,  5.31021633e+05, -6.06995609e+05, 3.87230616e+05, -1.05486503e+05]

def Fiducial(r,dt_us):
    if abs(dt_us) > 79 and abs(dt_us) < 789:
        r_max = 0
        for i in range(len(polyParams)):
            r_max += polyParams[i]*pow(dt_us/1000., i)
        if r < r_max/10.:
            return True
        else:
            return False

def get_files_from_list(path):
    with open(path, 'r') as file:
        file_paths = file.read().split('\n')
    return [path for path in file_paths if path]

#ER/NR Bands
def get_mdc3_band(dat_path):
    #Bin Center, Bin Actual, Gaus Mean, Mean Error, Gaus 90CL, 90CL Error, X^2/DOF
    df = pd.read_csv(dat_path, sep='\t', usecols= [0,1,2,3,4], names= ['bin_center', 'S1_bin_actual','log10_S2_mean','mean_error','logS2_90CL'])
    df['90CLupper'] = df['log10_S2_mean'] + df['logS2_90CL']
    df['90CLlower'] = df['log10_S2_mean'] - df['logS2_90CL']
    return df

def plot_mdc3_band(df,label,**kwargs):
    kwargs.update({'linestyle': '-'})
    kwargs.update({'alpha': 0.5})
    linewidth = kwargs.setdefault('linewidth',4)
    plt.plot(df['S1_bin_actual'], df['log10_S2_mean'], label='{} mean'.format(label), **kwargs)
    kwargs.update({'linestyle': '--', 'linewidth': linewidth-2})
    plt.plot(df['S1_bin_actual'], df['90CLupper'], label='{} 90 CL Upper'.format(label), **kwargs)
    plt.plot(df['S1_bin_actual'], df['90CLlower'], label='{} 90 CL lower'.format(label), **kwargs)
    return

def NR_Cut(s1, s2, df):
    CL90u = df['90CLupper']
    CL90l = df['90CLlower']
    i = 0
    if 0 < s1 < 999:
        while s1 > i:
            i += 1
        else:
            if s2 < CL90u[i] and s2 > CL90l[i]:
                return True
            else:
                return False
    else:
        return False


ER_Band = get_mdc3_band('../Input_Files/MDC3_ER_band.dat')
NR_Band = get_mdc3_band('../Input_Files/MDC3_NR_band.dat')
    
#######################################################################################################################
#Load in LZAP Files
#######################################################################################################################
file_prefix = 'root://gfe02.grid.hep.ph.ic.ac.uk'
lzap_files = []
lzap_file_list = get_files_from_list('../Input_Files/lzap_files.txt')
for i in range(len(lzap_file_list)):
    lzap_files.append(file_prefix+lzap_file_list[i])

if sys.argv[2] == 'full':
    num_files = len(lzap_files)
else:
    num_files = int(sys.argv[2])
print(type(num_files),num_files)

lzap_files = lzap_files[0:num_files]

for n in range (len(lzap_files)):
    rootfile = lzap_files[n]
    eventTree = uproot.open(rootfile)['Events']
    scattersTree = uproot.open(rootfile)['Scatters']

    eventBranches = ['eventHeader.runID', 'eventHeader.eventID', 'eventHeader.rawFileName']
    tpcBranches = ['pulsesTPC.nPulses', 'pulsesTPC.pulseStartTime_ns', 'pulsesTPC.pulseEndTime_ns', 'pulsesTPC.pulseArea_phd']
    ssBranches = ['ss.nSingleScatters', 'ss.s1Area_phd', 'ss.s2Area_phd', 'ss.driftTime_ns', 'ss.x_cm', 'ss.y_cm', 'ss.correctedX_cm', 'ss.correctedY_cm', 'ss.correctedS1Area_phd', 'ss.correctedS2Area_phd', 'ss.s1PulseID', 'ss.s2PulseID', 'ss.skinTotalArea', 'ss.odPromptArea']

    if n == 0:
        
        eventDataLZAP=eventTree.arrays(eventBranches)
        tpcDataLZAP=eventTree.arrays(tpcBranches)
        ssDataLZAP=scattersTree.arrays(ssBranches)
        
        runID_LZAP_list = list(eventDataLZAP[b'eventHeader.runID'])
        eventID_LZAP_list = list(eventDataLZAP[b'eventHeader.eventID'])
        fileName_LZAP_list = list(eventDataLZAP[b'eventHeader.rawFileName'])
        
        nTPCpulses_LZAP_list = list(tpcDataLZAP[b'pulsesTPC.nPulses'])
        pulseArea_LZAP_list = list(tpcDataLZAP[b'pulsesTPC.pulseArea_phd'])
        
        nSS_LZAP_list = list(ssDataLZAP[b'ss.nSingleScatters'])
        s1Area_LZAP_list = list(ssDataLZAP[b'ss.s1Area_phd'])
        s2Area_LZAP_list = list(ssDataLZAP[b'ss.s2Area_phd'])
        driftTime_LZAP_list = list(ssDataLZAP[b'ss.driftTime_ns'])
        x_LZAP_list = list(ssDataLZAP[b'ss.x_cm'])
        y_LZAP_list = list(ssDataLZAP[b'ss.y_cm'])
        s1cArea_LZAP_list = list(ssDataLZAP[b'ss.correctedS1Area_phd'])
        s2cArea_LZAP_list = list(ssDataLZAP[b'ss.correctedS2Area_phd'])

    else:
        
        eventDataLZAP=eventTree.arrays(eventBranches)
        tpcDataLZAP=eventTree.arrays(tpcBranches)
        ssDataLZAP=scattersTree.arrays(ssBranches)
        
        runID_LZAP_list += list(eventDataLZAP[b'eventHeader.runID'])
        eventID_LZAP_list += list(eventDataLZAP[b'eventHeader.eventID'])
        fileName_LZAP_list += list(eventDataLZAP[b'eventHeader.rawFileName'])
        
        nTPCpulses_LZAP_list += list(tpcDataLZAP[b'pulsesTPC.nPulses'])
        pulseArea_LZAP_list += list(tpcDataLZAP[b'pulsesTPC.pulseArea_phd'])
        x_LZAP_list += list(ssDataLZAP[b'ss.x_cm'])
        y_LZAP_list += list(ssDataLZAP[b'ss.y_cm'])        
        nSS_LZAP_list += list(ssDataLZAP[b'ss.nSingleScatters'])
        s1Area_LZAP_list += list(ssDataLZAP[b'ss.s1Area_phd'])
        s2Area_LZAP_list += list(ssDataLZAP[b'ss.s2Area_phd'])
        driftTime_LZAP_list += list(ssDataLZAP[b'ss.driftTime_ns'])

        s1cArea_LZAP_list += list(ssDataLZAP[b'ss.correctedS1Area_phd'])
        s2cArea_LZAP_list += list(ssDataLZAP[b'ss.correctedS2Area_phd'])

        
#        print(len(fileName_LZAP_list),len(runID_LZAP_list))
print('formed header lists')


#######################################################################################################################
#Load MCTruth Files
#######################################################################################################################
mctruth_files = []
mctruth_file_list = get_files_from_list('../Input_Files/mctruth_files.txt')
for i in range(len(mctruth_file_list)):
    mctruth_files.append(file_prefix+mctruth_file_list[i])
mctruth_files=mctruth_files[0:num_files]
    
for n in range(len(mctruth_files)):

    mctruth_rootfile = mctruth_files[n]
    mcTruthTree=uproot.open(mctruth_rootfile)['RQMCTruth']

    shortEventBranches = [ 'mcTruthEvent.runNumber', 'mcTruthEvent.derEvent', 'mcTruthEvent.baccEvent', 'mcTruthEvent.parentParticle', 'mcTruthEvent.parentVolume']
    eventBranches = [ 'mcTruthEvent.runNumber', 'mcTruthEvent.derEvent', 'mcTruthEvent.baccEvent', 'mcTruthEvent.parentParticle', 'mcTruthEvent.parentVolume', 'mcTruthEvent.eventID', 'mcTruthEvent.runID']
    pulseBranches = ['mcTruthPulses.nRQMCTruthPulses', 'mcTruthPulses.pulseIdentifier', 'mcTruthPulses.vertexNumber', 'mcTruthPulses.pheCount', 'mcTruthPulses.firstPheTime_ns', 'mcTruthPulses.lastPheTime_ns']
    vertexBranches = ['mcTruthVertices.nRQMCTruthVertices', 'mcTruthVertices.volumeName', 'mcTruthVertices.particleName', 'mcTruthVertices.time_ns', 'mcTruthVertices.positionX_mm', 'mcTruthVertices.positionY_mm', 'mcTruthVertices.positionZ_mm', 'mcTruthVertices.energyDep_keV', 'mcTruthVertices.detectedS1Photons', 'mcTruthVertices.s1PhotonHits', 'mcTruthVertices.s2PhotonHits', 'mcTruthVertices.rawS1Photons', 'mcTruthVertices.rawS2Photons', 'mcTruthVertices.detectedS2Photons', 'mcTruthVertices.s1PulseIndex', 'mcTruthVertices.s2PulseIndex', 'mcTruthVertices.detectedScintPhotons']

    if n == 0:
        fileName_list =[]
        try:
            eventData=mcTruthTree.arrays(eventBranches)
        except:
            eventData=mcTruthTree.arrays(shortEventBranches)
        pulseData=mcTruthTree.arrays(pulseBranches)
        vertexData=mcTruthTree.arrays(vertexBranches)
        
        try:
            eventID_list = list(eventData[b'mcTruthEvent.eventID'])
            runID_list = list(eventData[b'mcTruthEvent.runID'])
        except:
            eventID_list = [0]*len(list(eventData[b'mcTruthEvent.runNumber']))
            runID_list = [0]*len(list(eventData[b'mcTruthEvent.runID']))
        parentParticle_list = list(eventData[b'mcTruthEvent.parentParticle'])
        parentVolume_list = list(eventData[b'mcTruthEvent.parentVolume'])
    
        nPulses_list = list(pulseData[b'mcTruthPulses.nRQMCTruthPulses'])
        
        nVertices_list = list(vertexData[b'mcTruthVertices.nRQMCTruthVertices'])
        volumeName_list = list(vertexData[b'mcTruthVertices.volumeName'])
        particleName_list = list(vertexData[b'mcTruthVertices.particleName'])
        vertexTime_list = list(vertexData[b'mcTruthVertices.time_ns'])
        posX_list = list(vertexData[b'mcTruthVertices.positionX_mm'])
        posY_list = list(vertexData[b'mcTruthVertices.positionY_mm'])
        posZ_list = list(vertexData[b'mcTruthVertices.positionZ_mm'])
        energyDep_list = list(vertexData[b'mcTruthVertices.energyDep_keV'])
        detectedS1Photons_list = list(vertexData[b'mcTruthVertices.detectedS1Photons'])
        detectedS2Photons_list = list(vertexData[b'mcTruthVertices.detectedS2Photons'])
        detectedScintPhotons_list = list(vertexData[b'mcTruthVertices.detectedScintPhotons'] )
        for j in range(len(vertexData[b'mcTruthVertices.nRQMCTruthVertices'])):
            fileName_list.append(lzap_files[n]) 
    else:
        try:
            eventData=mcTruthTree.arrays(eventBranches)
        except:
            eventData=mcTruthTree.arrays(shortEventBranches)
        pulseData=mcTruthTree.arrays(pulseBranches)
        vertexData=mcTruthTree.arrays(vertexBranches)
        
        try:
            eventID_list += list(eventData[b'mcTruthEvent.eventID'])
            runID_list += list(eventData[b'mcTruthEvent.runID'])
        except:
            eventID_list += [0]*len(list(eventData[b'mcTruthEvent.runNumber']))
            runID_list += [0]*len(list(eventData[b'mcTruthEvent.runID']))
        parentParticle_list += list(eventData[b'mcTruthEvent.parentParticle'])
        parentVolume_list += list(eventData[b'mcTruthEvent.parentVolume'])
    
        nPulses_list += list(pulseData[b'mcTruthPulses.nRQMCTruthPulses'])
        
        nVertices_list += list(vertexData[b'mcTruthVertices.nRQMCTruthVertices'])
        volumeName_list += list(vertexData[b'mcTruthVertices.volumeName'])
        particleName_list += list(vertexData[b'mcTruthVertices.particleName'])
        vertexTime_list += list(vertexData[b'mcTruthVertices.time_ns'])
        posX_list += list(vertexData[b'mcTruthVertices.positionX_mm'])
        posY_list += list(vertexData[b'mcTruthVertices.positionY_mm'])
        posZ_list += list(vertexData[b'mcTruthVertices.positionZ_mm'])
        energyDep_list += list(vertexData[b'mcTruthVertices.energyDep_keV'])
        detectedS1Photons_list += list(vertexData[b'mcTruthVertices.detectedS1Photons'])
        detectedS2Photons_list += list(vertexData[b'mcTruthVertices.detectedS2Photons'])
        detectedScintPhotons_list += list(vertexData[b'mcTruthVertices.detectedScintPhotons'] )
        for j in range(len(vertexData[b'mcTruthVertices.nRQMCTruthVertices'])):
            fileName_list.append(lzap_files[n])

#    print(n)
print('loaded mctruth lists')

#######################################################################################################################
#Define Entry Specific Variables
#######################################################################################################################
s1, s2, volumeName, particleName, scintPhotons = 0,0,0,0,0
s1raw, s2raw, s1hits, s2hits = 0,0,0,0,
x,y,z = 0,0,0
energy = 0
nVertices = 0
runID, eventID, parentParticle, parentVolume = 0,0,0,0
baccEvent, derEvent, runNumber = 0,0,0
fileName = 0
nPulses, pulseType, pulseVertexNumber, pulsePheCount = 0,0,0,0

correctedS1, correctedS2 = 0., 0.
skinTotalArea, odPromptArea = 0., 0.
driftTime = 0.
x_ss, y_ss, xc_ss, yc_ss = 0., 0., 0., 0.
nSS = 0
lzapRun, lzapEvent, lzapRawFile = 0, 0, 0

def GetEntry(i):
    global s1, s2, scintPhotons
    global s1raw, s2raw, s1hits, s2hits
    global particleName, volumeName
    global x, y, z
    global energy
    global nVertices
    global eventID, runID, parentParticle, parentVolume
    global baccEvent, derEvent, runNumber
    global fileName, lzapRun, lzapEvent, lzapRawFile
    global nPulses, pulseType, pulseVertexNumber, pulsePheCount
    global correctedS1, correctedS2, driftTime, x_ss, y_ss, xc_ss, yc_ss, nSS  #LZAP Quantities
    global skinTotalArea, odPromptArea
    runID = runID_list[i]
    eventID = eventID_list[i]
    parentParticle = parentParticle_list[i]
    parentVolume = parentVolume_list[i]
    s1, s2 = detectedS1Photons_list[i], detectedS2Photons_list[i]    #array of ints/longs
    particleName = particleName_list[i] #array of ints/longs
    volumeName = volumeName_list[i]     #array of ints/longs
    scintPhotons = detectedScintPhotons_list[i] #array of ints/longs
    x = posX_list[i]
    y = posY_list[i]
    z = posZ_list[i]
    energy = energyDep_list[i]
    nVertices = nVertices_list[i] #integer
    fileName = fileName_list[i]
    nPulses = nPulses_list[i]
    #LZAP info -- not truth
    nSS, correctedS1, correctedS2 = nSS_LZAP_list[i], s1cArea_LZAP_list[i], s2cArea_LZAP_list[i]
    x_ss, y_ss, driftTime = x_LZAP_list[i], y_LZAP_list[i], driftTime_LZAP_list[i]
    lzapRun, lzapEvent, lzapRawFile = runID_LZAP_list[i], eventID_LZAP_list[i], fileName_LZAP_list[i]

numEvents = len(nVertices_list)
print(numEvents)

#######################################################################################################################
#Count Single Scatters
#######################################################################################################################
S1c, logS2c = [], []
r2, drift = [], []
rc2, drift = [], []
nSS_observed = 0
for n in range(numEvents):
    GetEntry(n)
    if nSS == 1:
        r = (np.sqrt(x_ss*x_ss + y_ss*y_ss))
        dt = (driftTime/1000.)
        S1c.append(correctedS1)
        logS2c.append(np.log10(correctedS2))
        r2.append((x_ss*x_ss + y_ss*y_ss))
        rc2.append((xc_ss*xc_ss + yc_ss*yc_ss))
        drift.append(driftTime/1000.)
        nSS_observed += 1
print("%i Single Scatters Seen! %.3f %% of events!" % (nSS_observed, nSS_observed*100./numEvents))

#######################################################################################################################
#Define Vertex Specific variables
#######################################################################################################################
vert_S1_sum, vert_logS2_sum = [], []
r2_true, z_true, r2_trueA, r2_trueB, z_trueA, z_trueB = [], [], [], [], [], []
ind_r, ind_z, ind_r_print, ind_z_print = [], [], [], []
r_true = []
mssi_S1c, mssi_logS2c = [], []
mssi_r2, mssi_drift = [], []
ss_S1c, ss_logS2c = [], []
ss_r2, ss_drift = [], []
allr, allz, allS1, allS2, allR2, allDrift = [], [], [], [], [], []
nMSSI, ss_events, NR_events = 0, 0, 0
cut_s1, cut_s2, cut_R2, cut_drift = [], [], [], []
ss_particle, mssi_particle = [], []

#######################################################################################################################
#Iterate over Verticies looking for energy dep in liquid xenon and RFR
#######################################################################################################################
print('NumEvents=', numEvents)
for n in range(numEvents):
    GetEntry(n)
    rad2_cm = (x_ss*x_ss + y_ss*y_ss)
    if nVertices > 0 and nSS == 1 and Fiducial(np.sqrt(rad2_cm), driftTime/1000) == True:
        allS1.append(correctedS1)
        allS2.append(np.log10(correctedS2))
        allR2.append(rad2_cm)
        allDrift.append(driftTime/1000)
        if NR_Cut(correctedS1, np.log10(correctedS2), NR_Band) == True:
            NR_events += 1
            cut_s1.append(correctedS1)
            cut_s2.append(np.log10(correctedS2))
            cut_R2.append(rad2_cm)
            cut_drift.append(driftTime/1000)
            vetoScint = 0
            thisS1, thisS2 = 0, 0
            nFFRs, nRFRs = 0, 0
            particles = []
            theseRadii, theseZ, theseZ2 = [], [], []
            for i in range(nVertices):
                if any([scint > 0 for scint in scintPhotons]):
                    vetoScint += scintPhotons[i]
                if energy[i] > 0 and 'LiquidXenonTarget' in str(volumeName[i]):
                    nFFRs += 1
                    thisS1 += s1[i]
                    thisS2 += s2[i]
                    theseRadii.append(np.sqrt(x[i]*x[i] + y[i]*y[i]))
                    theseZ.append(z[i])
                    particles.append(str(particleName[i]))
                if energy[i] > 0 and 'ReverseFieldRegion' in str(volumeName[i]):
                    thisS1 += s1[i]
                    thisS2 += s2[i]
                    nRFRs += 1
                    allr.append(np.sqrt(x[i]*x[i] + y[i]*y[i]))
                    allz.append(z[i])
                    theseRadii.append(np.sqrt(x[i]*x[i] + y[i]*y[i]))
                    theseZ.append(z[i])
            if (nFFRs >= 1 and nRFRs >= 1):
                nMSSI += 1 
#                print( "MSSI:", correctedS1, np.log10(correctedS2), theseRadii, theseZ, vetoScint, particleName)
                x, y = 0, 0
                for i in range(len(theseZ)):
                    if theseZ[i] > 0:
                        x = 1
                    if theseZ[i] < 0:
                        y = 1
                    if x == 1 and y == 1:
                        del theseRadii[i+1:]
                        del theseZ[i+1:]
                        break
                ind_r.append(theseRadii)
                ind_z.append(theseZ)
                mssi_S1c.append(correctedS1)
                mssi_logS2c.append(np.log10(correctedS2))
                mssi_r2.append(x_ss*x_ss + y_ss*y_ss)
                mssi_drift.append(driftTime/1000)
                mssi_particle.append(particleName[-1])
            else:
                ss_events += 1
                ss_S1c.append(correctedS1)
                ss_logS2c.append(np.log10(correctedS2))
                ss_r2.append(x_ss*x_ss + y_ss*y_ss)
                ss_drift.append(driftTime/1000)
                ss_particle.append(particleName[-1])
                #break
print(ss_events)
print(nMSSI)
print(NR_events)
print("%i MSSI Events seen! %.3f %% of SingleScatters, %.3f %% of events!" % (nMSSI, nMSSI*100./nSS_observed, nMSSI*100./numEvents))

#######################################################################################################################
#Write Output to Files
#######################################################################################################################
if write == 1:
    Output_Folder = '../Output_Files/NR_Band_MSSI_Events_v5'
    print('writing files')
    with open(Output_Folder+'/Single_Scatters.txt', 'w') as filehandle:
        filehandle.write('ss_S1c, ss_logS2c, ss_R2, ss_DriftTime_us\n')
        for i in range(len(allS1)):
            filehandle.write('%f, %f, %f, %f\n' % (allS1[i], allS2[i], allR2[i], allDrift[i]))
        filehandle.close()

    with open(Output_Folder+'/Single_Scatters_w_NR_cut.txt', 'w') as filehandle:
        filehandle.write('ss_S1c, ss_logS2c, ss_R2, ss_DriftTime_us\n')
        for i in range(len(cut_s1)):
            filehandle.write('%f, %f, %f, %f\n' % (cut_s1[i], cut_s2[i], cut_R2[i], cut_drift[i]))
        filehandle.close()

    with open(Output_Folder+'/MSSI_vertices.txt', 'w') as filehandle:
        for x in range(len(ind_r)):
            for i in range(len(ind_r[x])):
                filehandle.write('%f, %f,' % (ind_r[x][i], ind_z[x][i]))
            filehandle.write('\n')
        filehandle.close()

    with open(Output_Folder+'/Single_Scatters_w_NR_MSSI_cuts.txt', 'w') as filehandle2:
        filehandle2.write('mssi_S1c, mssi_logS2c, mssi_r2, mssi_drift, mssi_particle\n')
        for i in range(len(mssi_S1c)):
            filehandle2.write('%f, %f, %f, %f, %s\n' % (mssi_S1c[i], mssi_logS2c[i], mssi_r2[i], mssi_drift[i], mssi_particle[i]))
        filehandle2.close()

    with open(Output_Folder+'/Single_Scatters_w_NR_not_MSSI_cuts.txt', 'w') as filehandle2:
        filehandle2.write('ss_S1c, ss_logS2c, ss_r2, ss_drift, ss_particle\n')
        for i in range(len(ss_S1c)):
            filehandle2.write('%f, %f, %f, %f, %s\n' % (ss_S1c[i], ss_logS2c[i], ss_r2[i], ss_drift[i], ss_particle[i]))
        filehandle2.close()
