#everything here assumes initial scan is high to low


#physical constants
F      = 96485   # [=] C/mol, Faraday's constant
R      = 8.3145  # [=] J/mol-K, ideal gas constant

molarMasses = {"H": 1.008, "B": 10.81, "LI": 6.94, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "NA": 22.990, "MG": 24.305, "P": 30.974, "S":32.06, "CL": 35.45, "K": 39.098, "CA": 40.078, "FE": 55.845, "ZN": 65.38, "BR": 79.904}


def normalizedF(temp):
    return F/(R * temp)


def getMolarMass(formula): # enter formula as C,6;H,12;O,6
    elements = formula.split(";")
    weight = 0
    for element in elements:
        atoms = element.split(",")
        weight += int(atoms[1]) * self.molarMasses[atoms[0]]
    return weight

def importCV(foldername, datestrings, electrodes, speeds, analyte, dataArray):
    for date in datestrings:
        for m in range(len(electrodes)):
            electrode = electrodes[m]
            for n in range(len(speeds)):
                speed = speeds[n]
                str = date + '-' + electrode + '-' + analyte + '-' + speed + 'mvs.txt'
                newObject = CV(electrode, analyte, date, speed)
                potentials = []
                currents = []
                try:
                    fileObject = open(foldername+'/'+str, 'r')
                    lines = fileObject.readlines()
                    i = 0
                    while (i < len(lines)):
                        line = lines[i]
                        if line[0:6] == 'Init E':
                            initE = line.split(" = ")[1]
                            initE = float(initE.strip())
                            newObject.setInitE(initE)
                        if line[0:4] == 'High':
                            highE = line.split(' = ')[1]
                            highE = float(highE.strip())
                            newObject.setHighE(highE)
                        if line[0:5] == 'Low E':
                            lowE = line.split(' = ')[1]
                            lowE = float(lowE.strip())
                            newObject.setLowE(lowE)
                        if line[0:9] == 'Potential':
                            i += 2
                            while (i < len(lines)):
                                point = lines[i]
                                points = point.split(', ')
                                potentials.append(float(points[0]))
                                currents.append(float(points[1]))
                                i +=1
                        i += 1      
                    newObject.inputCVPoints(potentials, currents)
                    dataArray[n, m] = newObject
                except:
                    print("no " + str)

def importCVtoArray(foldername, datestrings, electrodes, speeds, analyte):
    import numpy
    dataArray = numpy.ndarray((len(speeds),len(electrodes)), dtype = CV)
    importCV(foldername, datestrings, electrodes, speeds, analyte, dataArray)
    return dataArray

class CV:
    def __init__(self, electrode, analyte, date, scanRate):
        self.analyte = analyte
        self.data = date
        self.electrode = electrode
        self.initE = None
        self.highE = None
        self.lowE = None
        self.scanRate = scanRate #note that scan rate is in mv/sec
        self.Ep1 = None
        self.Ep2 = None
        self.ip1 = None
        self.ip2 = None
        self.potentials = None
        self.currents = None
        self.halfLength = None
        self.upperPotential = None
        self.lowerPotential = None
        self.upperCurrent = None
        self.lowerCurrent = None
    
    def represent(self):
        return self.electrode + " electrode in " + self.analyte + " at " + self.scanRate + " mV/s"
    
    def setInitE(self, initE):
        self.initE = initE
    
    def setHighE(self, highE):
        self.highE = highE
    
    def setLowE(self, lowE):
        self.lowE = lowE
    
    
    def inputCVPoints(self, potential, current):
        import numpy
        self.potentials = potential
        self.currents = current
        self.halfLength = int(len(current)/2)
        self.lowerPotential = potential[:self.halfLength]
        self.upperPotential = potential[self.halfLength:]
        self.lowerCurrent = current[:self.halfLength]
        self.upperCurrent = current[self.halfLength:]
        try:
            self.Ep1 = self.correctedEp1()
        except:
            self.Ep1 = None
        try:
            self.Ep2 = self.correctedEp2()
        except:
            self.Ep2 = None
        try:
            self.ip1 = self.correctedMinCurrent()
        except:
            self.ip1 = None
        try:
            self.ip2 = self.correctedMaxCurrent()
        except:
            self.ip2 = None
        
    def Ehalf(self):
        if self.Ep1 == None or self.Ep2 == None:
            raise Exception("inadequate data")
        return (self.Ep1 + self.Ep2)/2
        
    def graph(self, graph):
        import matplotlib.pyplot as plt
        if self.potentials == None or self.currents == None:
            raise Exception("no scan data")  
        graph.plot(self.potentials, self.currents)
        graph.set_xlabel('potential (V)')
        graph.set_ylabel('current (Amp)')
    
    def normalize(self, array):
        range = max(array) - min(array)
        return [item/range for item in array]

    def differentiate(self, xarray, yarray):
        if len(xarray) != len(yarray):
            raise Exception("x and y arrays must be the same length")
        derivativesArray = []
        for i in range(len(xarray)-1):
            derivativesArray.append((yarray[i+1] - yarray[i])/(xarray[i+1] - xarray[i]))
        return derivativesArray
    
    def smootherSG(self, y, window_size, order, deriv = 0, rate = 1):
        import numpy
        from math import factorial
        y = numpy.array(y)
        try:
             window_size = numpy.abs(int(window_size))
             order = numpy.abs(int(order))
        except ValueError as msg:
             raise ValueError("window_size and order have to be of type int")
        if (window_size % 2 != 1 or window_size < 1):
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = numpy.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = numpy.concatenate((firstvals, y, lastvals))
        return numpy.convolve( m[::-1], y, mode='valid')
    
    def findZero(self, array):
        outputArray = []
        for i in range(len(array) -1):
            if (array[i+1]>0 and array[i]<0) or (array[i+1]<0 and array[i]>0):
                outputArray.append(i)
        return outputArray
    
    def upperSegmentInflectionIndex(self):
        import numpy
        xnorm = self.normalize(self.upperPotential)
        ynorm = self.normalize(self.upperCurrent)
        firstDeriv = self.differentiate(xnorm, ynorm)
        smoothFirstDeriv = self.smootherSG(firstDeriv, 101, 3)
        truncateSmoothFirstDerive= smoothFirstDeriv[100:] 
        #zeros = self.findZero(truncateSmoothFirstDerive)
        extreme = numpy.argmax(ynorm)
        truncateSmoothFirstDeriveTwo = truncateSmoothFirstDerive[:extreme]
        firstDerivMinusOne = [item -1 for item in truncateSmoothFirstDeriveTwo]
        return self.findZero(firstDerivMinusOne)[0] +100
      
    def upperSegmentInflectionPoint(self):
        i = self.upperSegmentInflectionIndex()
        return self.upperPotential[i], self.upperCurrent[i]
    
    def lowerSegmentInflectionIndex(self):
        import numpy
        xnorm = self.normalize(self.lowerPotential)
        ynorm = self.normalize(self.lowerCurrent)
        firstDeriv = self.differentiate(xnorm, ynorm)
        smoothFirstDeriv = self.smootherSG(firstDeriv, 101, 3)
        truncateSmoothFirstDerive= smoothFirstDeriv[100:] 
        #zeros = self.findZero(truncateSmoothFirstDerive)
        extreme = numpy.argmin(ynorm)
        truncateSmoothFirstDeriveTwo = truncateSmoothFirstDerive[:extreme]
        firstDerivMinusOne = [item -1 for item in truncateSmoothFirstDeriveTwo]
        return self.findZero(firstDerivMinusOne)[0] + 100
    
    def lowerSegmentInflectionPoint(self):
        i = self.lowerSegmentInflectionIndex()
        return self.lowerPotential[i], self.lowerCurrent[i]
    
    def lowerCorrectionLine(self, x):
        x2, y2 = self.lowerSegmentInflectionPoint()
        x1 = self.lowerPotential[-1]
        y1 = self.lowerCurrent[-1]
        m = (y2 - y1)/(x2 - x1)
        return m*(x - x1) + y1
    
    def upperCorrectionLine(self, x):
        x1 = self.upperPotential[-1]
        y1 = self.upperCurrent[-1]
        x2, y2 = self.upperSegmentInflectionPoint()
        m = (y2 - y1)/(x2 - x1)
        return m*(x - x1) + y1
    
    def correctedLowerCurrent(self):
        import numpy
        return numpy.array(self.lowerCurrent) - self.lowerCorrectionLine(numpy.array(self.lowerPotential))
    
    def correctedEp1(self):
        import numpy
        i = self.lowerSegmentInflectionIndex()
        array = self.correctedLowerCurrent()[i:]
        j = numpy.argmin(array) + i
        return self.lowerPotential[j]
    
    def correctedMinCurrent(self):
        return min(self.correctedLowerCurrent()[self.lowerSegmentInflectionIndex():])
    
    def correctedUpperCurrent(self):
        import numpy
        return numpy.array(self.upperCurrent) - self.upperCorrectionLine(numpy.array(self.upperPotential))
    
    def correctedEp2(self):
        import numpy
        i = self.upperSegmentInflectionIndex()
        array = self.correctedUpperCurrent()[i:]
        j = numpy.argmax(array) + i
        return self.upperPotential[j]
        
    def correctedMaxCurrent(self):
        return max(self.correctedUpperCurrent()[self.upperSegmentInflectionIndex():])
        
    def graphInflectionPoints(self, graph, col):
        inflectLowerX, inflectLowerY = self.lowerSegmentInflectionPoint()
        inflectUpperX, inflectUpperY = self.upperSegmentInflectionPoint()
        graph.plot(inflectUpperX, inflectUpperY, '.', ms = 10, color = col)
        graph.plot(inflectLowerX, inflectLowerY, '.', ms = 10, color = col)
    
    def graphLowerLine(self, graph):
        import numpy
        graph.plot(self.lowerPotential, self.lowerCorrectionLine(numpy.array(self.lowerPotential)), 'k-')
    
    def graphUpperLine(self, graph):
        import numpy
        graph.plot(self.upperPotential, self.upperCorrectionLine(numpy.array(self.upperPotential)), 'k-')

#making python simulation of code
#k0 = electrochemical rate constant, units of cm/s
#kc = chemical rate constant, units of 1/s
#O refers to the oxidizing agent
#assume that O & R, the reducing agent, have the same diffusion constant
#units of diffusion constant are cm^2/s
#this code was adapted from Peter M Attia
#scan rate here is in V/s
#concentrations are in molar
#L is number of iterations
#DM is model diffusion coefficient
        
class simulateCV:
    def __init__(self, initialConcO = 1, ORDiffusionConstant = 1e-5, initialOverpotenital = 0.2, finalOverpotential = -0.2, scanRate = 1e-3,  numElectrons = 1, alpha = 0.5, k0 = 1e-2, kc = 1e-3, temp = 298.15, iterations = 500, DM = 0.45):
        import numpy
        self.initialConcO = initialConcO
        self.ORDiffusionConstant = ORDiffusionConstant
        self.initialOverpotential = initialOverpotenital
        self.finalOverpotential = finalOverpotential
        self.scanRate = scanRate
        self.numElectrons = numElectrons
        self.alpha = alpha
        self.k0 = k0
        self.kc = kc
        self.temp = temp
        self.L = iterations
        self.DM = DM
        self.updateConstants()
        self.simulate()
        
    def changeInitialConcO(self, newConc):
        self.initialConcO = newConc
        self.updateConstants()
        self.simulate()
        
    def changeORDiffusionConstant(self, diffusionConstant):
        self.ORDiffusionConstant = diffusionConstant
        self.updateConstants()
        self.simulate()
    
    def changeInitialOverpotential(self, startV):
        self.initialOverpotential = startV
        self.updateConstants()
        self.simulate()
        
    def changeFinalOverpotential(self, endV):
        self.finalOverpotential = endV
        self.updateConstants()
        self.simulate()
        
    def changeScanRate(self, v):
        self.scanRate = v
        self.updateConstants()
        self.simulate()
    
    def changeNumElectrons(self, n):
        self.numElectrons = n
        self.updateConstants()
        self.simulate()
        
    def changeAlpha(self, a):
        self.alpha = a
        self.updateConstants()
        self.simulate()
        
    def changek0(self, k):
        self.k0 = k
        self.updateConstants()
        self.simulate()
    
    def changekc(self, k):
        self.kc = k
        self.updateConstants()
        self.simulate()
        
    def changeTemp(self, T):
        self.temp = T
        self.updateConstants()
        self.simulate()
    
    def changeNumIterations(self, l):
        self.L = l
        self.updateConstants()
        self.simulate()
    
    def changeDM(self, dm):
        self.DM = dm
        self.updateConstants()
        self.simulate()
        
    def updateConstants(self):
        import numpy
        self.tk = 2 + (self.initialOverpotential - self.finalOverpotential)/self.scanRate
        self.Dt = self.tk/self.L
        self.Dx = numpy.sqrt(self.ORDiffusionConstant * self.Dt/self.DM) 
        self.j = numpy.ceil(4.2 * self.L**0.5) +5
        self.f = F/(R *self.temp)
        self.km = (self.kc * self.tk)/self.L
        if self.km > 0.1:
            raise Exception ('kc * tk /L = ' + str(round(self.km, 5)) + ", which exceeds the upper limit of 0.1")
    
    def simulate(self):
        import numpy
        C = self.initialConcO/1000
        k = numpy.array(range(self.L+1))
        t = self.Dt *k
        overpotentialVectorNeg = self.initialOverpotential - self.scanRate*t
        overpotentialVectorPos = self.finalOverpotential + self.scanRate * t
        eta1 = [item for item in overpotentialVectorNeg if item > self.finalOverpotential]
        eta2 = [item for item in overpotentialVectorPos if item <= self.initialOverpotential]
        overpotentialScan = numpy.array(eta1 + eta2)
        overpotentialNorm = self.f * overpotentialScan
        kf = self.k0 * numpy.exp(-1 * self.alpha * self.numElectrons * overpotentialNorm)
        kb = self.k0 * numpy.exp((1-self.alpha)* self.numElectrons *overpotentialNorm)
        intJ = int(self.j)
        newL = len(overpotentialScan)
        O = C * numpy.ones([newL, intJ])
        R = numpy.zeros([newL, intJ])
        JO = numpy.zeros((1, newL))
        for i in range(0, newL-1):
            for l in range(1, intJ-1):
                O[i+1, l] = O[i, l] + self.DM*(O[i, l+1] + O[i, l-1] - 2*O[i,l])
                R[i+1, l] = R[i, l] + self.DM*(R[i, l+1] + R[i, l-1] - 2*R[i, l]) - self.km*R[i,l]
            JO[0, i+1] = (kf[i+1] * O[i+1, 1] - kb[i+1] * R[i+1, 1])/(1 + self.Dx/self.ORDiffusionConstant * (kf[i+1] + kb[i+1]))
            O[i+1, 0] = O[i+1, 1] - JO[0, i+1] * (self.Dx/self.ORDiffusionConstant)
            R[i+1, 0] = R[i+1, 1] + JO[0, i+1] *(self.Dx/self.ORDiffusionConstant) - self.km*R[i+1, 1]
        Z = -self.numElectrons * F * JO.transpose()
        if len(overpotentialScan) > len(Z):
            overpotentialScan = overpotentialScan[:-1]
        self.overpotentialScan = overpotentialScan
        self.currentDensity = Z
        
    def graph(self, graph):
        from matplotlib import pyplot as plt             
        graph.plot(self.overpotentialScan, self.currentDensity)
        graph.set_xlabel("Overpotential (E - E0') (V)")
        graph.set_ylabel("Current Density (A/cm^2)")


class simulationComparison: #assume in water
    
    def __init__(self, CVobject, concO, diffusionConstant, electrodeArea, temp = 298.15, numElectrons = 1):
        import numpy
        self.initialConcO = concO
        self.CV = CVobject
        self.scanRate = float(self.CV.scanRate)/1000
        self.ORDiffusionConstant = diffusionConstant
        self.temp = temp
        self.numElectrons = numElectrons
        self.redoxPotential = self.CV.Ehalf()
        self.initialOverpotential = self.CV.highE - self.redoxPotential
        self.finalOverepotential = self.CV.lowE - self.redoxPotential
        self.electrodeArea = electrodeArea
        self.simulation = simulateCV(self.initialConcO, self.ORDiffusionConstant, self.initialOverpotential, self.finalOverepotential, self.scanRate, self.numElectrons, temp = self.temp)
        self.expPotentials = self.CV.potentials        
        self.expCurrents = self.CV.currents
        self.uploadFromSim()
    
    def uploadFromSim(self):
        self.simulationPotentials = self.simulation.overpotentialScan + self.redoxPotential
        self.simulationCurrentDensity = self.simulation.currentDensity
        self.simulationCurrent = self.simulationCurrentDensity * self.electrodeArea
    
    def changeAlpha(self, alpha):
        self.simulation.changeAlpha(alpha)
        self.uploadFromSim()
    
    def changek0(self, k):
        self.simulation.changek0(k)
        self.uploadFromSim()
    
    def changekc(self, k):
        self.simulation.changekc(k)
        self.uploadFromSim()
    
    def normalizeExpCurrents(self):
        import numpy
        range = numpy.max(self.expCurrents) - numpy.min(self.expCurrents)
        return self.expCurrents/range
    
    def normalizeSimCurrents(self):
        import numpy
        range = numpy.max(self.simulationCurrentDensity) - numpy.min(self.simulationCurrentDensity)
        return self.simulationCurrentDensity/range
    
    def graphSimulationCurrentDensity(self, graph):
        self.simulation.graph(graph)
    
    def graphExperiment(self, graph):
        self.CV.graph(graph)
        
    def graphSimulationCurrent(self, graph):
        import matplotlib.pyplot as plt
        graph.plot(self.simulationPotentials, self.simulationCurrent)
        graph.set_xlabel('potential (V)')
        graph.set_ylabel('current (AU)')
        
    def graphNormalizedSimulation(self, graph):
        import matplotlib.pyplot as plt
        graph.plot(self.simulationPotentials, self.normalizeSimCurrents())
        graph.set_xlabel('potential (V)')
        graph.set_ylabel('current (AU)')
    
    def graphNormalizedExp(self, graph):
        import matplotlib.pyplot as plt
        if self.expPotentials == None or self.expCurrents == None:
            raise Exception("no scan data")  
        graph.plot(self.expPotentials, self.normalizeExpCurrents())
        graph.set_xlabel('potential (V)')
        graph.set_ylabel('current (AU)')
        
    def getExperimentalParamters(self):
        return "CV of " + self.CV.analyte + " with " + self.CV.electrode + " electrode in water at " + str(self.scanRate) + " V/s from " + str(round(self.CV.highE, 1)) + " to " + str(round(self.CV.lowE, 1)) + " V"
    
    def getSimulationParameters(self):
        return "initial conc O: " + str(self.initialConcO) + "\nOR diffusion constant: " + str(self.ORDiffusionConstant) + "\nalpha: " + str(self.simulation.alpha) + "\nk0: " + str(self.simulation.k0) + "\nkc: " + str(self.simulation.kc) + "\ninterations: " + str(self.simulation.L) + "\nDM: " + str(self.simulation.DM)