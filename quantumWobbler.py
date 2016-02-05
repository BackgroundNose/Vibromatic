# -*- coding: utf-8 -*-
from enthought.traits.api import *
from enthought.traits.ui.api import *
from chaco.api import ArrayPlotData, Plot, PlotGraphicsContext
from enable.component_editor import ComponentEditor
from pyface.api import FileDialog, confirm, YES

from Solver import *

class Calculation(HasTraits):
    """
    GUI front end for the "Quantum Wobbler" Solver.py module.

    Written by James Furness (2016)
    """
    M1 = Float(0)
    M2 = Float(0)
    D0 = Float(0)
    Re = Float()
    levels = String("0,1,2-6")
    levelsToFind = []
    loadPot = Button()

    muma = 0.0
    mu = 0.0

    helpme = Button()
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
                        Item('helpme', show_label=False, label="\nHelp Me!\n"),
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

    def _helpme_fired(self):
        self.add_line("Help information for quantum wobbler:", True)
        self.add_line("Potential file should be comma separated values reading:\n\"separation,value\"\n\"separation,value\"\n...\n")
        self.add_line("The potential should be normalised so it tends to 0 at long range.")
        self.add_line(" i.e. the energy of the isolated atoms should be subtracted from the energy of the interacting pair")
        self.add_line("In this way D0 is the energy difference (in Hartree) between the equilibrium, and dissociation.")
        self.add_line("Re is the equilibrium distance of the diatomic. If this is ommitted the program will estimate it as the separation of the lowest energy data point. This can lead to inaccuracy.")
        self.add_line("M1 and M2 are the mass of atom 1 and 2 respectively. The atom numbering is arbitrary.")
        self.add_line("Levels gives the levels that should be optimized. These are supplied as a comma separated list. Ranges can be specified with a hyphen.")


    def _loadPot_fired(self):
        self.Rlist = []
        self.Ulist = []
        self.morseList = []
        self.guessEigList = []
        self.convergedValues = []
        self.Plists = []
        self.muma = 0.0
        self.mu = 0.0

        dialog = FileDialog(title='Select a potential',action='open')

        if not dialog.open():
            self.add_line("I can't open the file \'" + str(dialog.path)+"\'")
            return

        self.updateFromDetails(readInput(dialog.path))
        self.add_line("Loading new potential!",True)
        self.add_line("Loading potential from: "+str(dialog.path))
        self._toInterest_fired()


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

    def storeDataInDetails(self):
        out = Details()

        out.Rlist = self.Rlist[:]
        out.Ulist = self.Ulist[:]
        out.re = self.Re
        out.mu = self.mu
        out.muma = self.muma
        out.M1 = self.M1
        out.M2 = self.M2
        out.D0 = abs(self.D0)
        out.levels = ""+self.levels
        out.guessEigVals = self.guessEigList[:]
        out.convergedValues = self.convergedValues[:]
        out.levelsToFind = self.levelsToFind[:]
        return out

    def updateFromDetails(self, data):
        self.guessEigList = data.guessEigVals[:]
        self.convergedValues = data.convergedValues[:]
        self.levelsToFind = data.levelsToFind[:]
        self.levels = ""+data.levels
        self.morseList = data.morseList[:]
        self.Ulist = data.Ulist[:]
        self.D0 = data.D0
        self.M1 = data.M1
        self.M2 = data.M2
        self.Rlist = data.Rlist[:]
        self.levelsToFind = data.levelsToFind[:]
        self.mu = data.mu
        self.muma = data.muma



    def _solve_fired(self):
        
        self.add_line("Beginning solver",title=True)

        self.D0 = abs(self.D0)

        if self.dataMissing():
            return

        if self.Re == 0:
            self.add_line("WARNING! : No Re set. Guessing from potential!")
            self.Re = self.Rlist[self.Ulist.index(min(self.Ulist))]
            self.add_line("Re set as: "+str(self.Re))

        self.muma = reducedMass(self.M1, self.M2)
        self.mu = reducedMass(convertMassToAU(self.M1), convertMassToAU(self.M2))
        self.add_line("Converting Potential and D0 to 2 x reduced mass scale.")
        self.add_line("Reduced mass (au) = "+str(self.mu))
        self.add_line("Reduced mass (ma) = "+str(self.muma))

        self.scaled = True

        data = self.storeDataInDetails()
# Parse the levels asked for. Return False on error.
        data.levelsToFind = parseLevels(self.levels, self)
        if len(data.levelsToFind) == 0:
            return

        data = driver(data, self)

        self.updateFromDetails(data)

        self.scaled = False

        self._toInterest_fired()


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
            if val < len(self.convergedValues):
                plotData.set_data("eig"+str(val), [self.convergedValues[val], self.convergedValues[val]])

        plot = Plot(plotData)

        if len(self.morseList) > 0:
            plot.plot(("x","morse"), type = "line", color = "red")

        for val in self.levelsToFind:
            if val < len(self.convergedValues):

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

    

    def greeting(self):
        self.add_line("==========================")
        self.add_line("     QUANTUM WOBBLER")
        self.add_line("==========================")
        self.add_line("\nA Diatomic vibrations solver.")
        self.add_line("By James Furness\n")

# Read the input file
inp = Calculation().configure_traits()
