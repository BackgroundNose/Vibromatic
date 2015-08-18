# -*- coding: utf-8 -*-
from enthought.traits.api import *
from enthought.traits.ui.api import *
from chaco.api import ArrayPlotData, Plot, PlotGraphicsContext
from enable.component_editor import ComponentEditor
from pyface.api import FileDialog, confirm, YES

from Solver import *

class Calculation(HasTraits):
    M1 = Float()
    M2 = Float()
    D0 = Float()
    Re = Float()
    levels = String("0,1,2-6")
    levelsToFind = []
    loadPot = Button()

    muma = 0.0
    mu = 0.0

    solve = Button()

    plot = Instance(Plot)
    toInterest = Button()

    results = String()
    saveLog = Button()
    clear = Button()

    savePlot = Button()

    scaled = False

    Rlist = []
    Ulist = []

    plotRangeX = []
    plotRangeY = []

    self.h = 0  # The separation between R points

    morseList = []
    guessEigList = []
    convergedValues = []
    Plists = []

    view = View(HSplit(
                    Group(
                        Item('loadPot',show_label=False, label="\nLoad Potential\n"),
                        Item(name='M1'),
                        Item(name='M2'),
                        Item(name='D0'),
                        Item(name='Re'),
                        Item(name='levels'),
                        Item('solve', show_label=False, label="\nSolve the System!\n")
                    ),
                    Group(
                        Item(name='results', show_label=False,style='custom'),
                        Item(name='clear', show_label=False,label="Clear output"),
                        Item(name='saveLog', show_label=False,label="Save output")
                    ),
                    Group(
                        Item('plot', editor = ComponentEditor(), show_label=False),
                        Item('toInterest',show_label=False, label="Scale Plot"),
                        Item('savePlot',show_label=False, label="Save Plot")
                    )
                ),
                title = "Quantum Wobbler",
                width = 800,
                height = 400,
                resizable = True)


    def __init__(self):
        self.greeting()

    def _loadPot_fired(self):
        self.Rlist = []
        self.Ulist = []
        self.morseList = []
        self.guessEigList = []
        self.convergedValues = []
        self.Plists = []
        self.muma = 0.0
        self.mu = 0.0
        try:
            dialog = FileDialog(title='Select a potential',action='open')
            if not dialog.open():
                self.add_line("I can't open the file \'" + str(dialog.path)+"\'")
                return
            inp = open(dialog.path,'r')
            self.potential = []
            for line in inp:
                raw = line.strip().split(',')
                raw = filter(None, raw)

                try:
                    R = float(raw[0])
                    U = float(raw[1])
                    self.Rlist.append(R)
                    self.Ulist.append(U)

                except ValueError:
                    self.add_line("I cant convert this to two floats:")
                    self.add_line("\""+str(line)+"\"")

            self.plotRU()
            self.add_line("Loading new potential!",True)
            self.add_line("Loading potential from: "+str(inp.name))
        except IOError:
            self.add_line("I can't find the file \'" + str(dialog.path)+"\'")

    def _toInterest_fired(self):
        self.D0=abs(self.D0)
        if len(self.Rlist) == 0:
            self.add_line("Load a potential first!")
            return
        if self.D0 == 0:
            self.add_line("Set a D0 first!")
            return
        if (min(self.Ulist) < 0):
            if len(self.morseList) > 0 and min(self.morseList) >= 0:
                setUtoDe(self.morseList, self.D0)
            self.plotRangeY = (-self.D0, self.D0/10.0)
        else:
            self.plotRangeY = (0, self.D0*1.1)
        self.plotRU()

    def _savePlot_fired(self):
        dialog = FileDialog(title='Save plot as...', action='save as',
                            wildcard = FileDialog.create_wildcard('Portable Network Graphics', '*.png'))
        if not dialog.open():
            self.add_line("ERROR opening save file.")
            return
        self.plotRU(save=True,filename=dialog.path)

    def _clear_fired(self):
        if confirm(None, "Are you sure you want to clear the output?") == YES:
            self.results = ""
            self.greeting()

    def _saveLog_fired(self):
        dialog = FileDialog(title='Save plot as...', action='save as',
                            wildcard = FileDialog.create_wildcard('Text file', '*.txt'))
        if not dialog.open():
            self.add_line("ERROR opening save file.")
            return
        try:
            out = open(dialog.path,'w')
            out.write(self.results)
            out.close()
        except IOError:
            self.add_line("I failed to open "+str(dialog.path))
            return

    def _solve_fired(self):
        # Parse the levels asked for. Return False on error.
        if not self.parseLevels():
            return
        self.add_line("Beginning solver",title=True)

        if self.dataMissing():
            return

        if self.Re == 0:
            self.add_line("WARNING! : No Re set. Guessing from potential!")
            self.Re = self.Rlist[self.Ulist.index(min(self.Ulist))]
            self.add_line("Re set as: "+str(self.Re))

        setUtoRe(self.Ulist)
        #self._toInterest_fired()

        self.muma = reducedMass(self.M1, self.M2)
        self.mu = reducedMass(convertMassToAU(self.M1), convertMassToAU(self.M2))
        self.add_line("Converting Potential and D0 to 2 x reduced mass scale.")
        self.add_line("Reduced mass (au) = "+str(self.mu))
        self.add_line("Reduced mass (ma) = "+str(self.muma))
        self.Ulist, self.D0 = convertUto2ReducedAU(self.mu, self.Ulist, self.D0)
        self.scaled = True

        generateStartingE(self)

        setUtoDe(self.Ulist, self.D0)

        if len(self.guessEigList) == 0:
            self.add_line("No bound Eigenvalues were found for guess potential.")
            self.add_line("For the hell of it, lets try optimizing something very nearly unbound...")

            # try looking for one right at the dissociation limit.
            self.guessEigList.append(-self.D0*0.95)

        self.h = self.Rlist[-1] - self.Rlist[-2]      # This is the separation between R points

        self.add_line("\nR Points Separation: "+str(self.h))

        self.convergedValues = [0.0]*len(self.guessEigList)
        for val in self.levelsToFind:
            if val >= len(self.guessEigList):
                self.add_line("\nERROR: No guess state for level "+str(val))
                self.add_line("Possbily state is unbound, manual guesses to be implimented eventually...")
            else:
                self.add_line("Optimizing level "+str(val),True)
                self.convergedValues[val], tlist = optimizeEigenvalue(self, self.Rlist, self.Ulist, self.guessEigList[val], self.h)
                self.Plists.append(tlist)

        self.Ulist, self.D0 = revertFrom2ReducedAU(self.mu, self.Ulist, self.D0)
        self.morseList, tmp = revertFrom2ReducedAU(self.mu, self.morseList, 0.0)
        self.scaled = False

        for val in self.levelsToFind:
            if val < len(self.guessEigList):
                if self.convergedValues[val] != "UNBOUND":
                    self.convergedValues[val] = enToEh(self.mu, self.convergedValues[val])

        self._toInterest_fired()

        self.add_line("Summary",True)
        # Final Print Roundup
        for val in self.levelsToFind:
            if val < len(self.guessEigList):
                if self.convergedValues[val] != "UNBOUND":
                    self.add_line("v"+str(val)+": " + str(self.convergedValues[val])+" Eh\t"+str(abs(self.convergedValues[val])*219474.63)+" cm-1")
                else:
                    self.add_line("v"+str(val)+" is not a bound state!")


    def add_line(self, string, title=False):
        if title:
            pad = "-"*len(string)
            self.results += ("\n\n"+pad+"\n"+ string.upper() + "\n"+pad+"\n\n")
        else:
            self.results += string + "\n"

    def plotRU(self,rangeX=None, rangeY=None, save=False, filename=""):
        if save and filename == "":
            self.add_line("ERROR: I need a valid file name")
            return

        if save and filename.split('.')[-1] != "png":
            self.add_line("ERROR: File must end in .png")
            return

        if len(self.morseList) > 0:
            plotData = ArrayPlotData(x=self.Rlist, y=self.Ulist, morse=self.morseList, eigX=[self.Rlist[0], self.Rlist[-1]])
        else:
            plotData = ArrayPlotData(x=self.Rlist, y=self.Ulist)

        for val in self.levelsToFind:
            if val < len(self.guessEigList):
                plotData.set_data("eig"+str(val), [self.convergedValues[val], self.convergedValues[val]])

        plot = Plot(plotData)

        if len(self.morseList) > 0:
            plot.plot(("x","morse"), type = "line", color = "red")

        for val in self.levelsToFind:
            if val < len(self.guessEigList):
                plot.plot(("eigX","eig"+str(val)), type="line", color="green")

        plot.plot(("x","y"), type = "line", color = "blue")
        plot.plot(("x","y"), type = "scatter", marker_size = 1.0, color = "blue")
        #
        plot.index_axis.title = "Separation (r0)"
        if (self.scaled):
            plot.value_axis.title = "Potential (Eh * 2 * mu)"
        else:
            plot.value_axis.title = "Potential (Eh)"

        if len(self.plotRangeX) != 0:
            plot.x_axis.mapper.range.low = self.plotRangeX[0]
            plot.x_axis.mapper.range.high = self.plotRangeX[1]
        if len(self.plotRangeY) != 0:
            plot.y_axis.mapper.range.low = self.plotRangeY[0]
            plot.y_axis.mapper.range.high = self.plotRangeY[1]
        if not save:
            self.plot = plot
        else:
            plot.outer_bounds = [800,600]
            plot.do_layout(force=True)
            gc = PlotGraphicsContext((800,600), dpi = 72)
            gc.render_component(plot)
            gc.save(filename)

    def dataMissing(self):
        fail = False
        if len(self.Rlist) == 0:
            self.add_line("ERROR: Load a potential first!")
            fail = True
        if self.D0 == 0:
            self.add_line("ERROR: Set a D0!")
            fail = True
        if self.M1 == 0:
            self.add_line("ERROR: Set a mass for atom 1!")
            fail = True
        if self.M2 == 0:
            self.add_line("ERROR: Set a mass for atom 2!")
            fail = True
        return fail

    def parseLevels(self):
        self.levelsToFind = []
        raw = self.levels.split(',')
        for item in raw:
            item = item.strip()
            if item.find('-') != -1:
                temp = item.split('-')
                try:
                    for i in range(int(temp[0]), int(temp[1])+1):
                        self.levelsToFind.append(i)
                except ValueError:
                    self.add_line("ERROR: Problem in level choice!")
                    self.add_line("Choose individual levels with comma separation")
                    self.add_line("or a range with start-stop")
                    return False
            else:
                try:
                    self.levelsToFind.append(int(item))
                except ValueError:
                    self.add_line("ERROR: Problem in level choice!")
                    self.add_line("Choose individual levels with comma separation")
                    self.add_line("or a range with start-stop")
                    return False

        return True

    def greeting(self):
        self.add_line("=============================")
        self.add_line("       QUANTUM WOBBLER")
        self.add_line("=============================")
        self.add_line("\nA Diatomic vibrations solver.")
        self.add_line("By James Furness\n")

# Read the input file
inp = Calculation().configure_traits()
