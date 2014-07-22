# Este script define cada una de las clases/steppables usados,
# (excepto la difusion externa que esta definida en el xml).

#AquÃ­ se importan todas las librerÃ­as necesarias para el funcionamiento del programa.

from PySteppables import *
import CompuCell
import sys
from random import *
import numpy as np
import math
import CompuCellSetup
from PySteppablesExamples import MitosisSteppableBase
            
#A continuaciÃ³n se muestran cada una de las siete clases definidas.


# 1. El creadorFilamento lo unico que hace es formar un filamento. 
#    Solo es llamado en el escenario de filamento (cuidando de haber apagado los 
#    steppables de forma en el XML).

class CreadorFilamento(SteppableBasePy):
	
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
   
 	def start(self):     
#Este bloque realiza un filamento con cÃ©lulas de 4 pixeles de longitud 
#  en una rejilla con 100 pixeles de longitud.
         y=1
         dic={}
         cell="cell"
         for i in np.arange(4, 98,4):
             dic[i]=cell+str(i)    
             print dic
         for i in np.arange(4,98,4):    
             dic[i]=self.potts.createCell()
             dic[i].type=self.ACTIVADA    
             self.cellField[i:i+4,y+5:y+6, 0]=dic[i]

    def step(self,mcs):
        pass

#2. Clase Directa. Esta es de los steppables  mÃ¡s importantes, ya que define la comunicaciÃ³n
#   entre cÃ©lulas y almacena las variables Ivec (suma de I en los vecinos) y 
#   Imul (concentraciÃ³n de I en la cÃ©lula multiplicada por el nÃºmero de vecinos).
#   Estas variables son almacendas para el resolvedor de  (EDO) de la clase SBML.

class Directa(SteppableBasePy):
	
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        pass

    def step(self,mcs):
#  Este loop es creado para evitar tomar las cÃ©lulas de la pared como vecinas.
        IdOfWall=[]
        for cell in self.cellListByType(self.WALL):
            IdOfWall.append(cell.id)
#En este loop, se busca cada una de las cÃ©lulas, y se les guarda el valor de I y de A.   
        for cell in self.cellListByType(self.ACTIVADA, self.INHIBIDA): 
            Ivec=0.0; nn=0; Imul=0.0   
            state=self.getSBMLState(_modelName='FI',_cell=cell)
            I= state['I']
            A = state['A']
# Para cada cÃ©lula,se buscan sus vecinos y se calcula cuÃ¡ntos son y cuÃ¡l es su suma de I 
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell): 
                if neighbor:
                    if neighbor.id in IdOfWall:
                        pass
                    else:
                        nn+=1
                        state=self.getSBMLState(_modelName='FI',_cell=neighbor)
                        Ivec+=state['I']    
# AquÃ­ se almacena en la variable Ivec: suma de I de los vecinos de la cÃ©lula
            state={} 
# Esta Ãºltima linea crea la variable Imul = I de la cÃ©lula X el nÃºmero de vecinos..       
            Imul= I*nn
# Las variables Imul e Ivec son guardadas en el diccionario de la cÃ©lula en la clase SBMLs       
            cellDict=self.getDictionaryAttribute(cell)
            cellDict["Ivec"]=Ivec
            cellDict["Imul"]=Imul


#3.Clase Indirecta. Esta clase:
#  a. Calcula la concentracion de I en el medio que estÃ¡ rodeando inmediatamente a la cÃ©lula 
#     y guarda la variable (Iext)     
#  b. Secreta una proporciÃ³n de I al medio, de manera distribuida a los pixeles
class Indirecta(SteppableBasePy): 
	
    pixelesFronteraMedio ={}
    boundaryPixelMedium= 0 
   #esta variable cuenta el nÃºmero de pixeles en la frontera de cada pixel. 

    def __init__(self,_simulator,_frequency): 
        SteppableBasePy.__init__(self,_simulator,_frequency) 
        self.boundaryStrategy=CompuCell.BoundaryStrategy.getInstance() 
        self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromNeighborOrder(1) 

    def start(self): 
        self.fieldI =self.getConcentrationField("Iext")
        self.fieldI[0,0,0]=2 
#aquÃ­ se crea un pequeÃ±o campo insignificante inicial de inhibidor 

    def step(self,mcs): 
# La funciÃ³n calcula el Iext de cada cÃ©lula (por medio de la funciÃ³n getIoutsideBoundary) 
# y lo guarda en el diccionaro de la cÃ©lula (para mandarlo al resolvedor, SBMLsolver).
        self.fieldI=self.getConcentrationField("Iext") 
        totalI={} 
        for cell in self.cellListByType(self.ACTIVADA, self.INHIBIDA):
            totalI[cell.id] = self.getIoutsideBoundary(cell)
 # concentraciÃ³n de I en el medio, en la frontera de cada cÃ©lula (Iext)
            cellDict=self.getDictionaryAttribute(cell)
 #Esta linea guarda los valores de Iext en el diccionario de la cÃ©lula
            cellDict["Iext"]=totalI[cell.id]

            state=self.getSBMLState(_modelName='FI',_cell=cell)
            I= state['I'] #concentraciÃ³n de I en cada cÃ©lula
            tasaSec= state['s'] # tasa de secreciÃ³n s
            S = I*tasaSec  #proporciÃ³n de I que sale de la celula 
            self.secreteOnMedium(cell,S) #la secreciÃ³n al medio.     


    def getIoutsideBoundary(self,cell): 
	
#esta funciÃ³n busca la concentraciÃ³n del campo (self.fieldI) en los pixeles de frontera.     
#Este pequeÃ±o loop crea una lista de de los pixeles vecinos, 
#es para evitar considerarlo en la bÃºsqueda de los pixeles que no son del medio

        neighPixelList=[]
        for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
            if neighbor:
                pixelNeigh= self.getCellPixelList(neighbor)
                for i in pixelNeigh:
                    a,b,c = i.pixel.x, i.pixel.y, i.pixel.z
                    neighPixelList.append([a,b,c])

        self.fieldI=self.getConcentrationField("Iext") 
        totalI=0 
        pixelList=self.getCellBoundaryPixelList(cell)
        visitedPixels=[]
        for bPixel in pixelList: 
            for i in xrange(self.maxNeighborIndex+1):
                pN=self.boundaryStrategy.getNeighborDirect(bPixel.pixel,i)
                x,y,z=pN.pt.x,pN.pt.y,pN.pt.z 

                if (([x,y,z] not in visitedPixels) and ([x,y,z] not in neighPixelList)):
                    cell2=self.cellField[x,y,z] 
                    if (not cell2 or cell2.id!=cell.id): 
                        totalI+=self.fieldI[x,y,z] 
                        visitedPixels.append([x,y,z])       
        iIndirecta.boundaryPixelMedium= len(visitedPixels)
        return totalI 

  #Esta funciÃ³n secreta la cantidad S de I al medio, de manera distribuida a los pixeles. 
    def secreteOnMedium(self,cell,S): 
        neighPixelList=[]
        for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
            if neighbor:
                pixelNeigh= self.getCellPixelList(neighbor)
                for i in pixelNeigh:
                    a,b,c = i.pixel.x, i.pixel.y, i.pixel.z
                    neighPixelList.append([a,b,c])
        pixelList=self.getCellBoundaryPixelList(cell)  
        visitedPixels=[] 
        for bPixel in pixelList: 
            for i in xrange(self.maxNeighborIndex+1): 
                pN=self.boundaryStrategy.getNeighborDirect(bPixel.pixel,i) 
                x,y,z=pN.pt.x,pN.pt.y,pN.pt.z 
                if (([x,y,z] not in visitedPixels) and ([x,y,z] not in neighPixelList)):
                    cell2=self.cellField[x,y,z] 
                    if (not cell2): 
                        Sporpixel = S/Indirecta.boundaryPixelMedium  

                        self.fieldI[x,y,z]+=Sporpixel
                        visitedPixels.append([x,y,z])

        if visitedPixels==[]:                
            Indirecta.pixelesFronteraMedio[cell.id]=0
        else:
            Indirecta.pixelesFronteraMedio[cell.id]=1
 
#4. Clase SBML. En esta clase se resuelve las EDO a cada paso montecarlo, 
#a partir de las variables generadas en las otras clases:
#	a. IVec, IMul de la clase DirectaClass
#	b. Iext de la clase Secretion
#	La EDO estÃ¡ definida en el archivo ComunicacionFinal.xml

class SBML(SteppableBasePy):    
    Aindi={}
    Iindi={}

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        import random 
        modelFile='Simulation/ComunicacionFinal.xml'  
        self.addSBMLToCellTypes(modelFile,'FI',[self.ACTIVADA, self.INHIBIDA],_stepSize=0.2) 
     #cargar los valores de SBML en cada celula
        state={} #diccionario para guardar los estados del SBML


  #En primer lugar,  se definen las condiciones iniciales, 
#  dependiendo del escenario (comunicaciÃ³n directa, indirecta o mixta) cambiamos los valores
        #en directa, s= 0, d= 0.009
        #en indirecta, s=0.009 y d = 0
        #en mixta, s=0.009, y d=0.009
        #condicion inicial aleatoria de activador (y por lo tanto de inhibidor).

        for cell in self.cellListByType(self.ACTIVADA, self.INHIBIDA):
            state['A'] = random.uniform(5,9) 
            state['I'] = 7
            state['X'] = 1   
            state['Ivec']=0
            state['Imul']=0
            state['Iext']=0
            state['e']= 0.000
            state['s']=0.000
            state['d']= 0.009
            self.setSBMLState(_modelName='FI',_cell=cell,_state=state)


    def step(self,mcs):
	    #aqui se definen los ids de las celulas que se va a seguir a lo largo del tiempo
        impulso = 7
        inhibited= [1,40]
        activated =[30,36]
        #aquÃ­ se sacan los valores de A y de I para cada cÃ©lula para poder graphicarla.  
        for cell in self.cellListByType(self.ACTIVADA, self.INHIBIDA):    
            state=self.getSBMLState(_modelName='FI',_cell=cell)           
            A = state['A']
            I = state['I'] 
            s = state['s']

            #este pequeÃ±o loop es para realizar una perturbaciÃ³n al paso 20000. 
            if mcs==20000:
                 if cell.id== 2 or cell.id == 42:
                     A = A+impulso  

           #esto es para seguir el inhibidor o activador de sÃ³lo 1 cÃ©lula, en el tiempo.
            if cell.id in activated or inhibited:
                SBML.Aindi[cell.id] = A
                SBML.Iindi[cell.id] = I

#En segundo lugar, se crea un diccionario para guardar las nuevas variables
            state={}
            cellDict=self.getDictionaryAttribute(cell)

#En tercer lugar, se llena el diccionario con los valores actualizados            
            try: 
                Ivec= cellDict['Ivec']
                Imul = cellDict['Imul']
                state['Ivec'] = Ivec
                state['Imul'] = Imul
            except KeyError:
                print "No Ivec o Imul"
                pass
            try: 
                Iext= cellDict['Iext']
                state['Iext'] = Iext
                s = Secretion.pixelesFronteraMedio[cell.id]*s
            except KeyError:
                print "no Iext"
                pass
            state['s'] = s
            state['A'] = A
            state['I'] = I

#En cuarto lugar, se calculan (es decir calcula la EDO) con las nuevas variables 
#de cada cÃ©lula.
            self.setSBMLState(_modelName='FI',_cell=cell,_state=state)

 #este Ãºltimo loop, cambia el tipo celular de cada cÃ©lula segÃºn su relaciÃ³n A/I.         

        for cell in self.cellListByType(self.ACTIVADA, self.INHIBIDA):
            activador = self.getSBMLValue(_modelName='FI',_valueName='A',_cell=cell)
            inhibidor = self.getSBMLValue(_modelName='FI',_valueName='I',_cell=cell)
            if inhibidor > activador:
                if cell.type== 1:
                    cell.type=2
            else:
                if cell.type==2:
                    cell.type=1  



#5. La clase contador TC actualiza cuantas celulas de cada tipo celular existen.
        
class ContadorTC(SteppableBasePy):
	ActivadasTotales=0
    InhibidasTotales=0

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

#Este comando genera la pared externa.
    def start(self):
		self.buildWall(self.WALL)

#Este bloque borra el numero de tipos celulares del paso montecarlo anterior y 
#agrega segun el tipo celular, 
#el nuevo numero de tipos celulares del paso montecarlo actual. 
    def step(self,mcs):
        contadorTC.ActivadasTotales = 0
        contadorTC.InhibidasTotales = 0
        for cell in self.cellListByType(self.INHIBIDA,self.ACTIVADA):
            if cell.type==1:
                contadorTC.ActivadasTotales += 1 
            elif cell.type==2:
                contadorTC.InhibidasTotales += 1

#6. Clase campoInhiActi. 
#Esta clase, lo Ãºnico que hace es permitir la visualizaciÃ³n de I y de A en cada cÃ©lula.

class CampoInhiActi(SteppableBasePy):

  def __init__(self,_simulator,_frequency=1):
      #aquÃ­ crea los campos de escalares.
    SteppableBasePy.__init__(self,_simulator,_frequency)
    self.scalarFieldActivador=CompuCellSetup.createScalarFieldCellLevelPy("Activador")
    self.scalarFieldInhibidor=CompuCellSetup.createScalarFieldCellLevelPy("Inhibidor")

  def step(self,mcs):  
      #AquÃ­ se borran los del paso anterior
    self.scalarFieldActivador.clear()
    self.scalarFieldInhibidor.clear()
    
    # Se busca en el diccionario de cada celula el valor de I y de A y se almacena
    for cell in self.cellListByType(self.ACTIVADA, self.INHIBIDA):
      if cell:
        cellDict=CompuCell.getPyAttrib(cell)
        activador =self.getSBMLValue(_modelName='FI',_valueName='A',_cell=cell) 
        inhibidor = self.getSBMLValue(_modelName='FI',_valueName='I',_cell=cell)
        
       #Se muestra el valor de I de A en la interfaz 
        self.scalarFieldActivador[cell]=activador
        self.scalarFieldInhibidor[cell]=inhibidor

                       
#7. Esta clase grafica todas las datos. 
#Se puede modificar segÃºn el tipo de grafo que se busque. 
                      
class Graficas(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
#La funciÃ³n start crea las grÃ¡ficas
 #estas dos listas cambian segÃºn el escenario. 
#Muestran que las cÃ©lulas que seran registradas en el tiempo.
        activated =[30,36]
        
#pW1 grafica los tipos celulares en el tiempo.
        self.pW1=self.addNewPlotWindow('Types of Cells','MonteCarlo Step (MCS)','Number')
        self.pW1.addPlot('Activada',_style='Lines',_color='blue',_size=2)
        self.pW1.addPlot('Inhibida',_style='Lines',_color='red',_size=2)
        
#pW3 grafica la A y I en una cÃ©lula cuyo estado estable es el Activado
        
        self.pW3=self.addNewPlotWindow('Activated Cells','MonteCarlo Step (MCS)','Number')
        for i in activated:
            self.pW3.addPlot("Act"+str(i),_style='Lines',_color='blue',_size=2)
            self.pW3.addPlot('Inh'+str(i),_style='Lines',_color='red',_size=2)
            
#pW4 grafica la A y I en una cÃ©lula cuyo estado estable es el inhibido

        self.pW4=self.addNewPlotWindow('Inhibited Cells','MonteCarlo Step (MCS)','Number')
        for i in inhibited:
            self.pW4.addPlot("Act"+str(i),_style='Lines',_color='blue',_size=2)
            self.pW4.addPlot('Inh'+str(i),_style='Lines',_color='red',_size=2)
                
        
    def step(self,mcs):
        #la funciÃ³n step va agregando los datos.
        inhibited= [1,40]
        activated =[30,36]
#         self.pW2.eraseAllData()

        self.pW1.addDataPoint("Activada",mcs,ConstraintInitializerSteppable.ActiveCells)
        self.pW1.addDataPoint("Inhibida",mcs,ConstraintInitializerSteppable.InhibitedCells)
        
        for i in activated:
            try:
                self.pW3.addDataPoint("Act"+str(i),mcs,SBMLsolver.Aindi[i] )
                self.pW3.addDataPoint("Inh"+str(i),mcs,SBMLsolver.Iindi[i] )
            except KeyError:
                pass         
        for i in inhibited:
            try:
                self.pW4.addDataPoint("Act"+str(i),mcs,SBMLsolver.Aindi[i] )
                self.pW4.addDataPoint("Inh"+str(i),mcs,SBMLsolver.Iindi[i] )
            except KeyError:
                pass            
        self.pW1.showAllPlots() 
        self.pW3.showAllPlots()
        self.pW4.showAllPlots() 
        
 #Este ultimo loop, guarda los datos de las grÃ¡ficas en un archivo de texto, 
#cada 100 pasos montecarlo.        
   
        if not mcs%100:
            self.pW1.savePlotAsData('ComIndiTipos.txt')
            self.pW3.savePlotAsData('IndiActivadas.txt')
            self.pW4.savePlotAsData('IndiInhibidas.txt')
