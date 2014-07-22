#  Este es el script principal. 
#De aquí se llama cada una de las clases definidas en adhesionComunicacionSteppables.py 
#  y se guardan los valores de las variables.

#1. Aquí se importan todas las librerías necesarias para el funcionamiento del programa.
import sys
from os import environ
from os import getcwd
import string
import CompuCellSetup
sys.path.append(environ["PYTHON_MODULE_PATH"])
sim,simthread = CompuCellSetup.getCoreSimulationObjects()
CompuCellSetup.initializeSimulationObjects(sim,simthread)
steppableRegistry=CompuCellSetup.getSteppableRegistry()

#2. La función coor_xy_tipos(X,Y,N, "Tipo1", "Tipo2") del programa piffaleatorio.py genera una distribución
#   aleatoria de N células de tipos celulares Tipo1 y Tipo2 en una rejilla de X por Y pixeles.
#   Esta función sólo se activa en el escenario de poblaciones celulares y de agregagación. 

from piffaleatorio import coor_xy_tipos
coor_xy_tipos(80, 80, 100, "Activada","Inhibida")


#3. Aquí se llaman las clases definidas en adhesionComunicacionSteppables.py.
#   Las clases ConstraintInitializerSteppable, SBMLsolver, ExtraFields y Plotting siempre están activadas. 
#   Las clases DirectaClass e IndirectaClass van a estar activadas o desactivadas según el escenario.                 

from adhesionComunicacionSteppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)

from adhesionComunicacionSteppables import DirectaClass
directaClass=DirectaClass(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(directaClass)

from adhesionComunicacionSteppables import IndirectaClass
secretion=Secretion(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(secretion)

from adhesionComunicacionSteppables import SBMLsolver
sbmlsolver=SBMLsolver(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(sbmlsolver)
        
from adhesionComunicacionSteppables import ExtraFields
extraFields=ExtraFields(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(extraFields)

from adhesionComunicacionSteppables import Plotting
plotting=Plotting(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(plotting)
  
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
