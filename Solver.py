import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import sys

"""
This module actually solves the nuclear Schrodinger equation for a diatomic.
"""

class Details():
    """Structure to hold the values read from input."""
    Rlist = []
    Ulist = []
    bkpUlist = []
    re = None
    mu = None
    muma = None
    M1 = None
    M2 = None
    D0 = None
    h = None
    levels = []
    guessEigVals = []
    convergedValues = []
    levelsToFind = []
    morseList = []
    Plists = []

def driver(data, calc):
    print "Beginning Quantum wobbler..."

    data.bkpUlist = data.Ulist[:]

    fig = generateStartingE(data, calc)

    setUtoDe(data.Ulist, data.D0)



    if len(data.guessEigVals) == 0:
        if calc is not None:
            calc.add_line("No bound Eigenvalues were found for guess potential.")
            calc.add_line("For the hell of it, lets try optimizing something very nearly unbound...")
        else:
            print "No bound Eigenvalues were found for guess potential."
            print "For the hell of it, lets try optimizing something very nearly unbound..."

        # try looking for one right at the dissociation limit.
        data.guessEigVals.append(-data.D0*0.95)

    data.h = data.Rlist[-1] - data.Rlist[-2]      # This is the separation between R points

    if calc is not None:
        calc.add_line("\nR Points Separation: "+str(data.h))
    else:
        print "\nR Points Separation: "+str(data.h)

    data.convergedValues = [0.0]*len(data.guessEigVals)
    for val in data.levelsToFind:
        if val >= len(data.guessEigVals):
            if calc is not None:
                calc.add_line("\nERROR: No guess state for level "+str(val))
                calc.add_line("Possbily state is unbound, manual guesses to be implimented eventually...")
            else:
                print "\nERROR: No guess state for level "+str(val)
                print "Possbily state is unbound, manual guesses to be implimented eventually..."
        else:
            if calc is not None:
                calc.add_line("Optimizing level "+str(val),True)
            else:
                print "======================"
                print " Optimizing level "+str(val)
                print "======================"
            data.convergedValues[val], tlist = optimizeEigenvalue(calc, data, data.guessEigVals[val])
            data.Plists.append(tlist)


    data.Ulist, tmp = revertFrom2ReducedAU(data.mu, data.Ulist, data.D0)
    data.morseList, data.D0 = revertFrom2ReducedAU(data.mu, data.morseList, data.D0)
    setUtoDe(data.morseList, data.D0)

    for val in data.levelsToFind:
        if val < len(data.convergedValues):
            if data.convergedValues[val] != "UNBOUND":
                data.convergedValues[val] = enToEh(data.mu, data.convergedValues[val])
    
    if calc is not None:
        calc.add_line("Summary",True)
    else:
        print "==========="
        print "  Summary"
        print "==========="

    for val in data.levelsToFind:
        if val < len(data.convergedValues):
            if data.convergedValues[val] != "UNBOUND":
                if calc is not None:
                    calc.add_line("v"+str(val)+": " + str(data.convergedValues[val])+" Eh\t"+str(abs(data.convergedValues[val])*219474.63)+" cm-1")
                else:
                    print "v"+str(val)+": " + str(data.convergedValues[val])+" Eh\t"+str(abs(data.convergedValues[val])*219474.63)+" cm-1"
            else:
                if calc is not None:
                    calc.add_line("v"+str(val)+" is not a bound state!")
                else:
                    print "v"+str(val)+" is not a bound state!"

    return data

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

def parseLevels(levels, calc=None):
    levelsToFind = []
    raw = levels.split(',')
    for item in raw:
        item = item.strip()
        if item.find('-') != -1:
            temp = item.split('-')
            try:
                for i in range(int(temp[0]), int(temp[1])+1):
                    levelsToFind.append(i)
            except ValueError:
                if calc is not None:
                    calc.add_line("ERROR: Problem in level choice!")
                    calc.add_line("Choose individual levels with comma separation")
                    calc.add_line("or a range with start-stop")
                else:
                    print "ERROR: Problem in level choice!"
                    print "Choose individual levels with comma separation"
                    print "or a range with start-stop"
                return []
        else:
            try:
                levelsToFind.append(int(item))
            except ValueError:
                if calc is not None:
                    calc.add_line("ERROR: Problem in level choice!")
                    calc.add_line("Choose individual levels with comma separation")
                    calc.add_line("or a range with start-stop")
                else:
                    print "ERROR: Problem in level choice!"
                    print "Choose individual levels with comma separation"
                    print "or a range with start-stop"
                return []

    return levelsToFind

def readInput(inpname):
    try:
        finput = open(inpname,'r')
    except IOError:
        print "ERROR: Failed to find file " + inpname
        sys.exit()

    out = Details()

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
                    out.D0 = getFloat(val, lineNo)

                elif (field == "M1"):
                    out.M1 = getFloat(val, lineNo)

                elif (field == "M2"):
                    out.M2 = getFloat(val, lineNo)

                elif (field == "LEVELS"):
                    out.levelsToFind = parseLevels(val)
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

                out.Rlist.append(getFloat(line[0], lineNo))
                out.Ulist.append(getFloat(line[1], lineNo))

    print "Potential consists of ", len(out.Ulist)," points."

    if out.re is None:
        out.re = out.Rlist[out.Ulist.index(min(out.Ulist))]

    if out.D0 is None:
        print "WARNING: Guessing at D0 as no input value found."
        print "         This may be inaccurate if potential is not long ranged enough."
        out.D0 = out.Ulist[-1]

    if out.M1 is None:
        print "INPUT ERROR: No mass found for atom 1.  e.g. M1 = 1.0"
        sys.exit()
    if out.M2 is None:
        print "INPUT ERROR: No mass found for atom 2.  e.g. M2 = 1.0"
        sys.exit()

    # Convert potential to be 0 at Re. (Currently always does this as no fitting is done).
    setUtoRe(out.Ulist)


    # Determine a.u. reduced mass for potential scaling
    out.mu = reducedMass(convertMassToAU(out.M1), convertMassToAU(out.M2))
    out.muma = reducedMass(out.M1, out.M2)

    # WHEN TESTING REMEMBER! THE MORSE POTENTIAL FOR COOLEY'S REFERENCE HAS THIS BUILT IN. SET MU to 0.5
    # deinf = 0.16454007507   # AND SET THIS TO 188.4355

    # Convert potential and de to a.u. scaled reduced mass
    out.Ulist, out.de = convertUto2ReducedAU(out.mu, out.Ulist, out.de)

    out.de = abs(out.de)

    return out

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

def plotFitting(calc, data, alpha):
    # Generate list of morse points
    data.morseList = []
    for r in data.Rlist:
        data.morseList.append(morse(alpha, data.re, data.D0, r))

    # else:
    #     fig = plt.figure()
    #     Uplot = fig.add_subplot(111)
    #     Uplot.set_xlim([min(data.Rlist), max(data.Rlist)])
    #     Uplot.set_ylim([min(data.Ulist), 1.1*data.de])
        
    #     Uplot.plot(data.Rlist, data.Ulist, color='green', label="Input Potential")
    #     Uplot.plot(data.Rlist, data.morseList, color='red', label="Fit Morse", linestyle='--')
        
    #     Uplot.set_xlabel("Separation ($r_0$)")
    #     Uplot.set_ylabel("Reduced Mass Scaled Potential (2*Eh*mu)")
    #     Uplot.legend()
        
    #     Uplot.set_title("Fit \nAlpha: "+str(alpha))

    #     return fig    

def assessFit( alpha, Rlist, Ulist, re, de):
    return morseResidual(Rlist, Ulist, re, de, alpha)

def fitMorse(data, calc):

    result =  minimize_scalar(assessFit, args=(data.Rlist, data.Ulist, data.re, data.D0))

    alpha = result.x
    if calc is not None:
        calc.add_line("Morse potential fit for starting guess")
        calc.add_line("Alpha: "+str(alpha))
        calc.add_line("Sum of squares: "+str(assessFit(alpha, data.Rlist, data.Ulist, data.re, data.D0)))
    else:
        print "Morse potential fit for starting guess"
        print "Alpha: "+str(alpha)
        print "Sum of squares: "+str(assessFit(alpha, data.Rlist, data.Ulist, data.re, data.D0))

    # None if GUI is being used. Otherwise contains the matplotlib figure
    fig = plotFitting(calc, data, alpha)

    return alpha, fig

def morseEigs(data, calc, alpha):
    # mass enters as reduced mass in Daltons for agreement with Cooley
    eig = 0

    k = int((4*math.pi*math.sqrt(2*data.muma*data.D0)/(2*math.pi*alpha)-1)/2.0)

    if calc is not None:
        calc.add_line("# Eigenstates: " + str(k))
        calc.add_line("Morse Eigenvalues", True)
    else:
        print "# Eigenstates: " + str(k)
        print "\n================="
        print "Morse Eigenvalues"
        print "=================\n"

    n = 0
    while  n < min(k,50):
        omega0 = (alpha/(2*math.pi))*math.sqrt(2*data.D0/data.muma)
        eig = -data.D0 + 2*math.pi*omega0*(n+0.5) - ((2*math.pi)**2*omega0**2) / (4.0*data.D0) * (n+0.5)**2

        data.guessEigVals.append(eig)
        if calc is not None:
            calc.add_line("v"+str(n)+": "+str(eig))
        else:
            print "v"+str(n)+": "+str(eig)
        n += 1
    
    if calc is not None:
        calc.guessEigList = data.guessEigVals[:]

    return 

def generateStartingE(data, calc=None):
    alpha, fig = fitMorse(data, calc)
    morseEigs(data, calc, alpha)
    return fig

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

def integrate(calc, data, E):
    if (E > data.Ulist[-1]):
        calc.add_line( "Unbound state!")
        return "UNBOUND"

    # Integrate inwards
    Pin = [0.0]*(len(data.Ulist))


    Pin[-1] = 1.0
    Pin[-2] = Pin[-1]*math.exp( data.Rlist[-1]*math.sqrt( data.Ulist[-1] - E ) - data.Rlist[-2]*math.sqrt(data.Ulist[-2] - E) )

    m = -1

    for i in reversed(range(0,len(Pin)-2)):
        Pin[i] = nextPin(i+1, data.Rlist, data.Ulist, E, data.h, Pin)

        if (Pin[i] < Pin[i+1]):
            Pin[i] = 0.0
            m = i+1
            break

    if (m == -1):
        if calc is not None:
            calc.add_line("ERROR in integration!")
            calc.add_line("Inward integration never hits maximum.")
        else:
            print "ERROR in integration!"
            print "Inward integration never hits maximum."
        return "ERROR"

    Pout = [0.0]*len(data.Rlist)
    Pout[0] = 0
    Pout[1] = 1.0
    for i in range(2,m+1):
        Pout[i] = nextPout(i-1, data.Rlist, data.Ulist, E, data.h, Pout)

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


def optimizeEigenvalue(calc, data, E):
    if calc is not None:
        calc.add_line("Initial Guess: "+str(E))
    else:
        print "Initial Guess: "+str(E)

    dE = 1E12
    tol = 1E-6

    it = 0
    while (abs(dE) >= tol):
        try:
            Pout, Pin, m = integrate(calc, data, E)
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

            dE = correctE( Pout, Pin, data.Ulist, E, data.h, m)

            E += dE
            if calc is not None:
                calc.add_line("\n-----------\n"+"Iteration "+str(it)+"\n-----------")
                calc.add_line("Correction: "+str(dE))
                calc.add_line("New E: "+str(E))
            else:
                print "\n-----------\n"+"Iteration "+str(it)+"\n-----------"
                print "Correction: "+str(dE)
                print "New E: "+str(E)
            it += 1
        except ValueError as e:
            return ("UNBOUND", "UNBOUND")
    if calc is not None:
        calc.add_line("\n\nConverged! At: "+str(E))
        calc.add_line("With dE: "+str(dE))
    else:
        print "\n\nConverged! At: "+str(E)
        print "With dE: "+str(dE)

    return E, Pout[:m]+Pin[m:]
