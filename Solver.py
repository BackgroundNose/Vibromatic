import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import sys

def convertUto2ReducedAU(mu, Ulist, deinf):
    for i in range(0,len(Ulist)):
        Ulist[i] *= 2*mu

    deinf *= 2*mu
    return Ulist, deinf

def revertFrom2ReducedAU(mu, Ulist, deinf):
    for i in range(0,len(Ulist)):
        Ulist[i] /= 2*mu

    deinf /= 2*mu
    return Ulist, deinf

def enToEh(mu, en):
    return en/(2*mu)

def getFloat(inp, lineNo):
    try:
        out = float(inp)
    except ValueError:
        print "INPUT ERROR: line ",lineNo
        print "Unable to read expected float."
        sys.exit()
    return out

def readInput(inpname):
    try:
        finput = open(inpname,'r')
    except IOError:
        print "ERROR: Failed to find file " + inpname
        sys.exit()

    Rlist = []
    Ulist = []
    re = None
    de = None

    M1 = None
    M2 = None

    lineNo = 0

    readingData = False

    for line in finput:
        lineNo += 1

        if (not readingData):
            if (line.find("[") != -1):
                readingData = True
                if (len(line.strip()) > 1):
                    print "INPUT ERROR: New line before data please! line", lineNo
                    sys.exit()
            else:
                line = line.upper()
                parts = line.strip().split('=')

                if (len(parts) != 2):
                    print "INPUT ERROR: I don't understand line", lineNo
                    sys.exit()

                field = parts[0].strip()
                val = parts[1].strip()


                if (field == "D0"):
                    de = getFloat(val, lineNo)

                #elif (field == "RE"):
                #    re = getFloat(val, lineNo)

                elif (field == "M1"):
                    M1 = getFloat(val, lineNo)

                elif (field == "M2"):
                    M2 = getFloat(val, lineNo)

                else:
                    print "WARNING: I don't understand the field " + field + " on line ", lineNo

        else:
            if (line.find("]") != -1):
                readingData = False
                if (len(line.strip()) > 1):
                    print "INPUT ERROR: New line before more input please! line", lineNo
                    sys.exit()
            else:
                line = line.strip().split(',')
                line = filter(None, line)

                if (len(line) != 2):
                    print "INPUT ERROR: Bad data line ", lineNo
                    print "Comma separated values for data please. e.g."
                    print "\tseparation, energy"
                    sys.exit()

                line[0] = line[0].strip()
                line[1] = line[1].strip()

                # Compensate for mathematica's annoying lack of 0. e.g. "1." for "1.0"
                if (line[0][-1] == '.'):
                    line[0] = line[0] + '0'
                if (line[1][-1] == '.'):
                    line[0] = line[0] + '0'

                Rlist.append(getFloat(line[0], lineNo))
                Ulist.append(getFloat(line[1], lineNo))

    print "Potential consists of ", len(Ulist)," points."

    if re is None:
        re = Ulist.index(min(Ulist))

    if de is None:
        print "WARNING: Guessing at D0 as no input value found."
        print "         This may be inaccurate if potential is not long ranged enough."
        de = Ulist[-1]

    if M1 is None:
        print "INPUT ERROR: No mass found for atom 1.  e.g. M1 = 1.0"
        sys.exit()
    if M2 is None:
        print "INPUT ERROR: No mass found for atom 2.  e.g. M2 = 1.0"
        sys.exit()

    # Convert potential to be 0 at Re. (Currently always does this as no fitting is done).
    setUtoRe(Ulist)


    # Determine a.u. reduced mass for potential scaling
    mu = reducedMass(convertMassToAU(M1), convertMassToAU(M2))

    # REMEMBER! THE MORSE POTENTIAL FOR COOLEY'S REFERENCE HAS THIS BUILT IN. SET MU to 0.5
    # deinf = 0.16454007507   # AND SET THIS TO 188.4355

    # Convert potential and de to a.u. scaled reduced mass
    Ulist, de = convertUto2ReducedAU(mu, Ulist, de)

    return Rlist, Ulist, mu, M1, M2, abs(de)

def convertMassToAU(M):
    return M*1822.88839

def reducedMass(M1, M2):
    return M1*M2/(M1+M2)

def rotationPotential(Rlist, Ulist, J, lam):
    for i in range(0,len(Rlist)):
        if Rlist[i] != 0:
            Ulist[i] += (J*(J+1) - lam**2)/Rlist[i]**2
        else:
            Ulist[i] = 1E12
    return Ulist

def setUtoRe(Ulist):
    smallest = min(Ulist)

    for i in range(0,len(Ulist)):
        Ulist[i] -= smallest

def setUtoDe(Ulist, deinf):
    for i in range(0,len(Ulist)):
        Ulist[i] -= deinf

def morse(alpha, re, de, r):
    return de*(1.0 - math.exp(-alpha*(r - re)))**2

def morseResidual(Rlist, Ulist, re, de, alpha):
    sumOSqrs = 0
    for i in range(0, len(Rlist)):
        if (Ulist[i] <= de):
            sumOSqrs += (morse(alpha, re, de, Rlist[i]) - Ulist[i])**2
    return sumOSqrs

def plotFitting(calc, alpha):
    # Generate list of morse points
    calc.morseList = []
    for r in calc.Rlist:
        calc.morseList.append(morse(alpha, calc.Re, calc.D0, r))
    calc._toInterest_fired()
    calc.plotRU()
    #fig = plt.figure()
    #Uplot = fig.add_subplot(111)
    #Uplot.set_xlim([min(Rlist), max(Rlist)])
    #Uplot.set_ylim([min(Ulist), 1.1*de])
    #
    #Uplot.plot(Rlist, Ulist, color='green', label="Input Potential")
    #Uplot.plot(Rlist, morseList, color='red', label="Fit Morse", linestyle='--')
    #
    #Uplot.set_xlabel("Separation ($r_0$)")
    #Uplot.set_ylabel("Reduced Mass Scaled Potential (2*Eh*mu)")
    #Uplot.legend()
    #
    #Uplot.set_title("Fit \nAlpha: "+str(alpha))

    #fig.show()


def assessFit( alpha, Rlist, Ulist, re, de):
    return morseResidual(Rlist, Ulist, re, de, alpha)

def fitMorse(calc):

    result =  minimize_scalar(assessFit, args=(calc.Rlist, calc.Ulist, calc.Re, calc.D0))

    alpha = result.x
    calc.add_line("Morse potential fit for starting guess")

    calc.add_line("Alpha: "+str(alpha))
    calc.add_line("Sum of squares: "+str(assessFit(alpha, calc.Rlist, calc.Ulist, calc.Re, calc.D0)))


    # TODO:  Convert this to use calc plotting!
    plotFitting(calc, alpha)
    return alpha

def morseEigs(calc, alpha):
    # mass enters as reduced mass in Daltons for agreement with Cooley
    eig = 0
    calc.guessEigList = []
    k = int((4*math.pi*math.sqrt(2*calc.muma*calc.D0)/(2*math.pi*alpha)-1)/2.0)

    calc.add_line("# Eigenstates: " + str(k))

    calc.add_line("Morse Eigenvalues", True)

    n = 0
    while  n < min(k,50):
        omega0 = (alpha/(2*math.pi))*math.sqrt(2*calc.D0/calc.muma)
        eig = -calc.D0 + 2*math.pi*omega0*(n+0.5) - ((2*math.pi)**2*omega0**2) / (4.0*calc.D0) * (n+0.5)**2
        calc.guessEigList.append(eig)

        calc.add_line("v"+str(n)+": "+str(eig))
        n += 1
    return calc.guessEigList

def generateStartingE(calc):
    alpha = fitMorse(calc)
    return morseEigs(calc, alpha)

def nextPout(i, Rlist, Ulist, E, h, Plist):
    if False:   # Simple Integration formula
        return h**2*(Ulist[i] - E)* Plist[i] + 2.0*Plist[i] - Plist[i-1]
    else:
        Yim1 = (1 - (h**2/12.0) * (Ulist[i-1] - E))*Plist[i-1]
        Yi = (1 - (h**2/12.0) * (Ulist[i] - E))*Plist[i]
        num = h**2*(Ulist[i] - E) * Plist[i] + 2*Yi - Yim1
        return num / (1 - (h**2/12.0)*(Ulist[i+1]-E))

def nextPin(i, Rlist, Ulist, E, h, Plist):
    if False:   # Simple Integration formula
        return h**2*(Ulist[i]-E)*Plist[i] + 2.0*Plist[i] - Plist[i+1]
    else:
        Yip1 = (1 - (h**2/12.0) * (Ulist[i+1] - E))*Plist[i+1]
        Yi = (1 - (h**2/12.0) * (Ulist[i] - E))*Plist[i]
        num = h**2*(Ulist[i] - E) * Plist[i] + 2*Yi - Yip1
        return num / (1 - (h**2/12.0)*(Ulist[i-1]-E))

def integrate(calc, Rlist, Ulist, E, h):
    if (E > Ulist[-1]):
        calc.add_line( "Unbound state!")
        return "UNBOUND"

    # Integrate inwards
    Pin = [0.0]*(len(Ulist))


    Pin[-1] = 1.0
    Pin[-2] = Pin[-1]*math.exp( Rlist[-1]*math.sqrt( Ulist[-1] - E ) - Rlist[-2]*math.sqrt(Ulist[-2] - E) )

    m = -1

    for i in reversed(range(0,len(Pin)-2)):
        Pin[i] = nextPin(i+1, Rlist, Ulist, E, h, Pin)

        if (Pin[i] < Pin[i+1]):
            Pin[i] = 0.0
            m = i+1
            break

    if (m == -1):
        calc.add_line("ERROR in integration!")
        calc.add_line("Inward integration never hits maximum.")
        return "ERROR"

    Pout = [0.0]*len(Rlist)
    Pout[0] = 0
    Pout[1] = 1.0
    for i in range(2,m+1):
        Pout[i] = nextPout(i-1, Rlist, Ulist, E, h, Pout)

    return Pout, Pin, m

def normaliseIntegration(Pout, Pin, m):
    for i in reversed(range(0,len(Pin))):
        Pin[i] /= Pin[m]

    for i in range(0,len(Pout)):
        Pout[i] /= Pout[m]

    return Pout, Pin

def getDerivs(Pout, Pin, Rlist, m):
    dPout = (Pout[m-1] - Pout[m]) / (Rlist[m-1] - Rlist[m])
    dPin = (Pin[m] - Pin[m+1]) / (Rlist[m] - Rlist[m+1])
    return dPout, dPin

def integrateP2(Pout, Pin, Rlist, m):
    total = 0.0
    for i in range(0, len(Rlist)-1):
        if i < m:
            P2 = ((Pout[i] + Pout[i+1])/2.0)**2
        else:
            P2 = ((Pin[i] + Pin[i+1])/2.0)**2
        total += P2 * (Rlist[i+1] - Rlist[i])
    return total

def simpleCorrectE( Pout, Pin, Rlist, m):
    dPout, dPin = getDerivs(Pout, Pin, Rlist, m)

    return (dPout - dPin) / integrateP2(Pout, Pin, Rlist, m)

def correctE( Pout, Pin, Ulist, E, h, m):
    sumP2 = 0.0
    for i in range(0, len(Pout)):
        if i < m:
            sumP2 += Pout[i]**2
        else:
            sumP2 += Pin[i]**2

    Ymm1 = (1 - (h**2/12.0) * (Ulist[m-1] - E))*Pout[m-1]
    Ym = (1 - (h**2/12.0) * (Ulist[m] - E))*Pout[m]
    Ymp1 = (1 - (h**2/12.0) * (Ulist[m+1] - E))*Pin[m+1]
    yBit = ((-Ymm1 + 2*Ym - Ymp1)/h**2)
    return (yBit + (Ulist[m] - E)*Pout[m])/sumP2


def optimizeEigenvalue(calc, Rlist, Ulist, E, h):
    calc.add_line("Initial Guess: "+str(E))

    dE = 1E12
    tol = 1E-6

    it = 0
    while (abs(dE) >= tol):
        try:
            Pout, Pin, m = integrate(calc, Rlist, Ulist, E, h)
            Pout, Pin = normaliseIntegration(Pout, Pin, m)

            #if False:
            #    Pfig = plt.figure()
            #    Pplot = Pfig.add_subplot(111)
            #    Pplot.set_xlim([min(Rlist), max(Rlist)])
            #    Pplot.set_ylim([min([min(Pin),min(Pout)]), max([max(Pin),max(Pout)])])
            #
            #    Pplot.set_title("Iteration: " +str(it))
            #    Pplot.plot(Rlist[:m], Pout[:m], color='green')
            #    Pplot.plot(Rlist[m:], Pin[m:], color='red')
            #
            #    Pfig.show()

            dE = correctE( Pout, Pin, Ulist, E, h, m)

            E += dE

            calc.add_line("\n-----------\n"+"Iteration "+str(it)+"\n-----------")
            calc.add_line("Correction: "+str(dE))
            calc.add_line("New E: "+str(E))

            it += 1
        except ValueError as e:
            return ("UNBOUND", "UNBOUND")

    calc.add_line( "\n\nConverged! At: "+str(E))
    calc.add_line("With dE: "+str(dE))
    return E, Pout[:m]+Pin[m:]
