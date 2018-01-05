"""

Fitness inference code for bulk fitness assay, by Atish Agarwala. Latest version of fitness inference algorithm as
described in: http://dx.doi.org/10.1016/j.cell.2016.08.002.

"""

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import norm


def inferFitness(barcodes,cycleTimes,allReads,outputFolder=None,experimentName=None,neutralBarcodes=None,
                 multNoiseThresh=3e3,zCutoff=2.7,multNoiseBase=0.1,lowCoverageThresh=5e5,sparsityThresh=None):
    """
    inferFitness - main fitness inference function. Expects barcodes, cycle times, and reads; returns dictionary of
    fitness inference data for all replicates. See Venkataram et. al. Cell 2016 for details of fitness assay algorithm.
    All fitnesses reported per cycle.

    Required parameters:

    :param barcodes: N x 1 list of all barcodes
    :param cycleTimes: 1 x q list of cycle times
    :param allReads: dictionary of length r, where key is replicate name, and value is N x q matrix of reads

    Optional parameters:
    :param outputFolder: folder to save fitness data output. Saved in tab separated column formats.
    :param experimentName: name of experiment, used in saving fitness files
    :param neutralBarcodes: list of putatively neutral (ie reference strain) barcodes
    :param multNoiseThresh: read threshold used to compute multiplicative noise coefficient
    :param zCutoff: cutoff used to establish neutrals
    :param multNoiseBase: default value of multiplicative noise, used if only 1 replicate exists
    :param lowCoverageThresh: threshold by which low coverage points are tossed out
    :param sparsityThresh: Remove timepoints where >sparsityThresh lineages have zero reads. Optional.

    :return repFitnessData: dictionary with same keys as allReads. values are dictionaries, with key value pairs:
    neutralBarcodes - N_{neut} x 1 list of barcodes which were assumed to be neutral.
    timePointsUsed - 1 x q list of timepoints used in inference
    multNoiseParams - 1 x q-1 list of multiplicative noise parameters
    kappas - 1 x q-1 list of additive noise parameters
    meanFitnesses - 1 x q-1 list of mean fitnesses
    allTimeFitness - N x q-1 list of fitnesses
    allTimeErrors - N x q-1 list of errors
    aveFitness - N x 1 list of overall fitness estimates
    aveError - N x 1 list of overall error estimates
    """

    repNames = list(allReads) # names of replicates

    # Make sure reads are a numpy array
    for repName in repNames:
        allReads[repName] = np.asarray(allReads[repName],dtype=np.dtype(float))
    cycleTimes = np.asarray(cycleTimes)

    # sort barcodes, rearrange data appropriately
    sortedBarcodeIdx = np.argsort(barcodes)
    barcodes = np.sort(barcodes)
    if neutralBarcodes is not None:
        neutralBarcodes = np.sort(neutralBarcodes)
    for repName in repNames:
        allReads[repName] = allReads[repName][sortedBarcodeIdx,:]

    # filter out low coverage/sparse timepoints
    filteredCycleTimes,filteredReads = filterTimepoints(cycleTimes,allReads,lowCoverageThresh,sparsityThresh)

    # set neutral indices (boolean array)

    if neutralBarcodes is None:
        # if no neutrals defined, use double filtering to figure them out
        neutralIndices = np.array(len(barcodes)*[True])
        dummyMean,dummyKappa,neutralIndices,zScores = meanVarAndNeutrals(neutralIndices,filteredReads[repNames[0]],zCutoff,filteredCycleTimes[repNames[0]],firstPass=True)
        dummyMean,dummyKappa,neutralIndices,zScores = meanVarAndNeutrals(neutralIndices,filteredReads[repNames[0]],zCutoff,filteredCycleTimes[repNames[0]])
        dummyMean,dummyKappa,neutralIndices,zScores = meanVarAndNeutrals(neutralIndices,filteredReads[repNames[0]],zCutoff,filteredCycleTimes[repNames[0]])
    else:
        neutralIndices = returnNeutralIndices(barcodes,neutralBarcodes)

    # estimate multiplicative noise
    multNoiseParams = inferMultNoise(filteredReads,filteredCycleTimes,multNoiseThresh,multNoiseBase)

    # infer fitnesses

    repFitnessData = inferFitnessAndError(filteredReads,filteredCycleTimes,multNoiseParams,neutralIndices,zCutoff,barcodes)

    # output consistency checks on multiplicative noise estimation
    print('Multiplicative noise consistency checks')
    for rep in repFitnessData:
        print('')
        print(rep,' inconsistent times:')
        print('')
        inconsistentTimes = ''
        for kappa,multNoiseParam,time in zip(repFitnessData[rep]['kappas'],repFitnessData[rep]['multNoiseParams'],
                                             repFitnessData[rep]['timePointsUsed'][:-1]):
            if kappa/multNoiseThresh>multNoiseParam:
                inconsistentTimes = inconsistentTimes+str(time)+', '
        if len(inconsistentTimes)==0:
            print('No clear inconsistencies')
        else:
            inconsistentTimes = inconsistentTimes[:-2]
            print(inconsistentTimes)
        print('')
    # save data

    if outputFolder!=None:
        if experimentName == None:
            experimentName = 'Exp0'
        outputFitnessData(repFitnessData,outputFolder,experimentName)

    return repFitnessData


def filterTimepoints(cycleTimes,allReads,lowCoverageThresh,sparsityThresh=None):

    """
    filterTimepoints - filters out low coverage timepoints, returns only high coverage data
    :param cycleTimes: 1 x q vector of cycle times for all replicates
    :param allReads: dictionary of length r, where key is replicate name, and value is N x q matrix of reads
    :param lowCoverageThresh: threshold by which low coverage points are tossed out
    :param sparsityThresh: Toss out timepoints where >sparsityThresh lineages have zero reads. Optional.
    :return filteredCycleTimes: dictionary length r with cycle times for each replicate
    :return filteredReads: dictionary length r with only high coverage timepoints
    """
    repNames = list(allReads)
    filteredCycleTimes = {}
    filteredReads = {}
    if sparsityThresh==None:
        sparsityCheck = False
    else:
        sparsityCheck = True
    for repName in repNames:
        currentCycleTimes = np.zeros(len(cycleTimes))
        currentCycleIdx = 0
        for t in range(0,len(cycleTimes)):
            totReads = np.sum(allReads[repName][:,t])
            sparsity = np.sum(allReads[repName][:,t]==0)/np.size(allReads[repName][:,t])
            if totReads>lowCoverageThresh and (not sparsityCheck or sparsity<sparsityThresh):
                currentCycleTimes[currentCycleIdx] = t
                currentCycleIdx = currentCycleIdx+1
        currentCycleTimes = currentCycleTimes[0:currentCycleIdx].astype('int')
        filteredCycleTimes[repName] = cycleTimes[currentCycleTimes]
        filteredReads[repName] = allReads[repName][:,currentCycleTimes]

    return filteredCycleTimes,filteredReads

def returnNeutralIndices(barcodes,neutralBarcodes):
    """
    returnNeutralIndices - matches neutral (relative to ancestor) barcodes with barcode list, returns indices of all
    matches. Assumes both barcode lists are previously sorted.
    :param barcodes: N x 1 list of barcodes
    :param neutralBarcodes: list of putatively neutral barcodes
    :return neutralIndices: boolean array of neutral barcodes
    """

    neutralIndices = np.array([False]*len(barcodes))
    neutralBarcodeIdx = 0
    barcodeIdx = 0

    while barcodeIdx<len(barcodes) and neutralBarcodeIdx<len(neutralBarcodes):
        if barcodes[barcodeIdx]==neutralBarcodes[neutralBarcodeIdx]:
            neutralIndices[barcodeIdx] = True
            neutralBarcodeIdx = neutralBarcodeIdx+1
        barcodeIdx = barcodeIdx+1

    return neutralIndices

def meanVarAndNeutrals(neutralIndices,replicateReads,zCutoff,cycleTimes,firstPass = False):
    """
    meanVarAndNeutrals - uses neutral indices to compute mean fitness, additive noise parameter and neutral indices
    :param neutralIndices: putatively neutral indices (logical index)
    :param replicateReads: N x q matrix of reads from a single replicate
    :param zCutoff: Cutoff of Z value used to identify neutrals
    :param cycleTimes: 1 x q vector of cycle times for replicate
    :param firstPass: if True, harsher conditioning in order to find neutrals
    :return meanFitness: 1 x q-1 vector of mean fitnesses
    :return kappas: 1 x q-1 vector of additive noise parameters
    :return newNeutralIndices: list of indices of lineages which pass new neutrality filter.
    """



    meanFitness = np.zeros(len(cycleTimes)-1)
    newNeutralIndices = np.zeros((len(neutralIndices),len(cycleTimes)-1))
    kappas = np.zeros(len(cycleTimes)-1)
    for initIdx,initTime in enumerate(cycleTimes[0:-1]):
        Rinit = np.sum(replicateReads[:,initIdx])
        Rfinal = np.sum(replicateReads[:,initIdx+1])
        skipInference = False # Skip inference if data not good


        # Use median for mean fitness estimation

        initReads =  np.median(replicateReads[neutralIndices, initIdx])
        finalReads = np.median(replicateReads[neutralIndices,initIdx+1])

        # If median is low due to sparseness, then use the mean instead.
        # If first pass, don't use this timepoint for inference.
        medianLowCutoff = 10  # if median below this value, use mean instead

        if initReads<medianLowCutoff or finalReads<medianLowCutoff:
            initReads = np.mean(replicateReads[neutralIndices, initIdx])
            finalReads = np.mean(replicateReads[neutralIndices, initIdx + 1])
            if firstPass:
                skipInference = True

        meanFitness[initIdx] = np.log(initReads/Rinit)\
                                -np.log(finalReads/Rfinal)
        expectedReads = (Rfinal/Rinit)*np.exp(-meanFitness[initIdx])*replicateReads[:,initIdx]
        zScores = replicateReads[:,initIdx+1]-expectedReads
        zScores = zScores*np.power(expectedReads,-0.5)

        # compute kappa
        kappas[initIdx] = np.var(zScores[neutralIndices*np.isfinite(zScores)])
        # new neutrals

        if firstPass:
            # use quartiles of distribution to get cutoff scale for first pass
            zMed = np.median(zScores[np.isfinite(zScores)])
            # quartile to use for inference
            pWidth = 0.25
            sortedZ = np.sort(zScores)
            minIdx = int(round((0.5-pWidth)*len(sortedZ)))
            maxIdx = int(round((0.5+pWidth)*len(sortedZ)))
            deltaZ = sortedZ[maxIdx]-sortedZ[minIdx]
            estWidth = deltaZ/(norm.ppf(0.5+pWidth)-norm.ppf(0.5-pWidth))
            for idx,isNeutral in enumerate(neutralIndices):
                if isNeutral:
                    if np.abs(zScores[idx]-zMed)<zCutoff*estWidth or skipInference:
                        newNeutralIndices[idx,initIdx] = 1
        else:
            # If information about neutrals known, use variance to estimate widths
            for idx,isNeutral in enumerate(neutralIndices):
                if isNeutral:
                    if np.abs(zScores[idx])<np.sqrt(kappas[initIdx])*zCutoff:
                        newNeutralIndices[idx,initIdx] = 1
        binVec = np.concatenate(([-1e8],np.linspace(-15,15,90),[1e8]))

        # normalize to mean fitness per cycle
        meanFitness[initIdx] = meanFitness[initIdx]/(cycleTimes[initIdx+1]-initTime)

    newNeutralIndices = np.prod(newNeutralIndices,1)
    newNeutralIndices = (newNeutralIndices == 1)



    return meanFitness,kappas,newNeutralIndices,zScores

def inferMultNoise(allReads,allCycleTimes,multNoiseThresh,multNoiseBase):
    """
    inferMultNoise - infer multiplicative noise parameter using high frequency lineages

    :param allReads: dictionary length r containing N x q read matrix
    :param allCycleTimes: dictionary length r containing cycle times of each replicate
    :param multNoiseThresh: lower read threshold for multiplicative noise inference
    :param multNoiseBase: value of multiplicative noise to use if only one replicate present
    :return multNoiseParams: dictionary length r containing 1 x q-1 vector of multiplicative noise parameter
    """

    # compute list of cycle times. Also compute vector to map calculated multiplicative noise parameters to cycle times.

    repNames = list(allReads)
    numLineages = np.shape(allReads[repNames[0]])[0]

    multNoiseParams = {}
    if len(repNames)>1:
        # find minimal matching set
        totalTimes = allCycleTimes[repNames[0]]
        for repName in repNames:
            totalTimes = np.intersect1d(totalTimes,allCycleTimes[repName])
        if len(totalTimes)>1:
            # find time indices and mapping from calculated parameters to cycle noise
            timeIndices = {}
            mappedTimePoints = {}
            for repName in repNames:
                timeIndices[repName] = np.zeros(len(totalTimes))
                mappedTimePoints[repName] = np.zeros(len(allCycleTimes[repName]))
                totalTimeIdx = 0
                specificTimeIdx = 0
                while totalTimeIdx<len(totalTimes):
                    specificTimeFound = False
                    while not specificTimeFound:
                        specificTimeFound = allCycleTimes[repName][specificTimeIdx]==totalTimes[totalTimeIdx]
                        if specificTimeFound:
                            timeIndices[repName][totalTimeIdx] = specificTimeIdx
                        if totalTimeIdx<=len(totalTimes)-2:
                            mappedTimePoints[repName][specificTimeIdx] = totalTimeIdx
                        else:
                            mappedTimePoints[repName][specificTimeIdx] = len(totalTimes)-2
                        specificTimeIdx = specificTimeIdx+1
                    totalTimeIdx = totalTimeIdx+1
                mappedTimePoints[repName][specificTimeIdx:] = mappedTimePoints[repName][specificTimeIdx-1]
                mappedTimePoints[repName] = mappedTimePoints[repName][0:-1]
            # compute noise parameters
            calculatedMultNoise = np.zeros(len(totalTimes)-1)
            for initIdx,initTime in enumerate(totalTimes[0:-1]):
                initReads = np.zeros((numLineages,len(repNames)))
                finalReads = np.zeros((numLineages,len(repNames)))
                repIdx = 0
                for repName in repNames:
                    initReads[:,repIdx] = allReads[repName][:,int(timeIndices[repName][initIdx])]
                    finalReads[:,repIdx] = allReads[repName][:,int(timeIndices[repName][initIdx])+1]
                    repIdx = repIdx+1
                initR = np.sum(initReads,axis=0)
                finalR = np.sum(finalReads,axis=0)

                filteredIdx = np.mean(initReads,axis=1)>multNoiseThresh
                # if nothing passes filter, use default value
                if sum(filteredIdx)==0:
                    calculatedMultNoise[initIdx] = multNoiseBase
                else:
                    # filtered frequencies
                    initFreq = initReads[filteredIdx,:]/initR
                    finalFreq = finalReads[filteredIdx,:]/finalR
                    # Variance of log slopes
                    logSlopes = np.log(finalFreq)-np.log(initFreq)
                    calculatedMultNoise[initIdx] = np.sqrt(np.mean(np.var(logSlopes,axis=1,ddof=1)))
            for repName in repNames:
                multNoiseParams[repName] = calculatedMultNoise[mappedTimePoints[repName].astype('int')]
        else:
            for repName in repNames:
                multNoiseParams[repName] = multNoiseBase*np.ones(len(allCycleTimes[repName])-1)
    else:

        multNoiseParams[repNames[0]] = multNoiseBase*np.ones(len(allCycleTimes[repNames[0]])-1)

    return multNoiseParams

def inferFitnessAndError(allReads,allCycleTimes,multNoiseParams,neutralIndices,zCutoff,barcodes):
    """
    inferFitnessAndError - Carries out fitness inference once neutrals identified and multiplicative noise parameter
    is estimated
    :param allReads: dictionary length r containing N x q read matrix
    :param allCycleTimes: dictionary length r containing cycle times of each replicate
    :param multNoiseParams: dictionary length r containing 1 x q-1 vector of multiplicative noise parameter
    :param neutralIndices: putatively neutral indices
    :param barcodes: - N x 1 list of barcodes
    :return repFitnessData: dictionary with same keys as allReads. values are dictionaries, with key value pairs:
    neutralBarcodes - N_{neut} x 1 list of barcodes which were assumed to be neutral.
    timePointsUsed - 1 x q list of timepoints used in inference
    multNoiseParams - 1 x q-1 list of multiplicative noise parameters
    kappas - 1 x q-1 list of additive noise parameters
    meanFitnesses - 1 x q-1 list of mean fitnesses
    allTimeFitness - N x q-1 list of fitnesses
    allTimeErrors - N x q-1 list of errors
    aveFitness - N x 1 list of overall fitness estimates
    aveError - N x 1 list of overall error estimates
    """
    repNames = list(allReads)
    repFitnessData = {}

    for repName in repNames:
        tempDataDict = {}


        # find mean fitness and kappa values

        meanFitness,kappas,newNeutralIndices,zScores = meanVarAndNeutrals(neutralIndices,allReads[repName],zCutoff
                                                                  ,allCycleTimes[repName])
        # Store neutrals
        tempDataDict['neutralBarcodes'] = barcodes[newNeutralIndices]
        # second pass, with cleaned up neutrals
        meanFitness,kappas,newNeutralIndices,zScores = meanVarAndNeutrals(newNeutralIndices,allReads[repName],zCutoff
                                                                  ,allCycleTimes[repName])

        # calculate fitness and error
        deltaTs = allCycleTimes[repName][1:]-allCycleTimes[repName][0:-1]
        totReads = np.sum(allReads[repName],0)

        allTimeFitness = np.log(allReads[repName][:,1:]/totReads[1:])-np.log(allReads[repName][:,0:-1]/totReads[0:-1])
        allTimeFitness = allTimeFitness/deltaTs
        allTimeFitness = allTimeFitness + meanFitness
        

        allTimeErrors = np.sqrt(np.power(allReads[repName][:,1:]/kappas,-1)
                                +np.power(multNoiseParams[repName],2))/(allCycleTimes[repName][1:]-allCycleTimes[repName][0:-1])

        # calculate ave fitness and error
        aveFitness,aveError = inverseVarAve(allTimeFitness,allTimeErrors)

        # save data
        tempDataDict['barcodes'] = barcodes
        tempDataDict['timePointsUsed'] = allCycleTimes[repName]
        tempDataDict['kappas'] = kappas
        tempDataDict['multNoiseParams'] = multNoiseParams[repName]
        tempDataDict['meanFitness'] = meanFitness
        tempDataDict['allTimeFitness'] = allTimeFitness
        tempDataDict['allTimeErrors'] = allTimeErrors
        tempDataDict['aveFitness'] = aveFitness
        tempDataDict['aveError'] = aveError
        tempDataDict['zScores'] = zScores
        repFitnessData[repName] = tempDataDict



    return repFitnessData

def inverseVarAve(meanVals,standardDevs):
    """
    inverseVarAve - take weighted average with inverse variances.
    :param meanVals: Values to be averaged, N x q. Averaged across second dimension
    :param standardDevs: Standard errors of each value, N x q
    :return weightedMeans: N x 1 vector of weighted average
    :return weightedStandardDevs: N x 1 vector of final standard error
    """

    weightedMeans = np.sum(meanVals*np.power(standardDevs,-2),axis=1)/np.sum(np.power(standardDevs,-2),axis=1)
    weightedStandardDevs = np.power(np.sum(np.power(standardDevs,-2),axis=1),-0.5)

    return weightedMeans,weightedStandardDevs

def outputFitnessData(repFitnessData,outputFolder,experimentName):

    """
    outputFitnessData - save fitness data to file
    :param repFitnessData: fitness inference data. See inferFitness for details
    :param outputFolder: folder to save fitness data output. Saved in tab separated column formats.
    :param experimentName: name of experiment, used in saving fitness files
    :return:
    """
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)
    fileHeader = outputFolder+'/'+experimentName
    for repName in repFitnessData.keys():
        # timepoints, mean fitness, multiplicative noise parameters, kappa

        repData = repFitnessData[repName]

        fileName = fileHeader+'inferenceStats'+str(repName)+'.txt'
        f = open(fileName,'w')
        f.write('Timepoints: '+' '.join("{:1.3f}".format(t) for t in repData['timePointsUsed'])+'\n')
        f.write('Mean fitness: '+' '.join("{:1.3f}".format(t) for t in repData['meanFitness'])+'\n')
        f.write('Mult noise (std dev/cycle): '+' '.join("{:1.3f}".format(t) for t in repData['multNoiseParams'])+'\n')
        f.write('Additive noise (var/cycle): '+' '.join("{:1.3f}".format(t) for t in repData['kappas'])+'\n')
        f.close()


        # all time fitness

        fileName = fileHeader+'allTimeFitness'+str(repName)+'.txt'
        np.savetxt(fileName,repData['allTimeFitness'],'%1.3f',' ')

        # all time errors

        fileName = fileHeader+'allTimeErrors'+str(repName)+'.txt'
        np.savetxt(fileName,repData['allTimeErrors'],'%1.3f',' ')

        # barcode fitness error

        fileName = fileHeader+'aveFitnessAndError'+str(repName)+'.txt'
        aveFitness = [ "{:1.3f}".format(x) for x in repData['aveFitness']]
        aveError = [ "{:1.3f}".format(x) for x in repData['aveError']]
        fitnessArray  = np.stack((repData['barcodes'],aveFitness,aveError))
        fitnessArray = np.transpose(fitnessArray)
        np.savetxt(fileName,fitnessArray,'%1s',delimiter=' ')

    return
