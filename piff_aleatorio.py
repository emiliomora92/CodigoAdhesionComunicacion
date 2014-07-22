#Este código simplemente crea una configuracion aleatoria de distribucion de celulas.

#Importamos las modulos importantes para el programa. 
import re
import sys
import random
from random import randrange
from random import choice

#Esta es la funcion que hace todo. Tiene cinco entradas.

#x se refiere al tamano en el eje x de la cuadrilla con la que se trabaja.
#y idem pero para el eje y.
# n es el numero de celulas que se quiere para la simulacion. las celulas miden un pixel, pero es importante considerar 

# ...que al primer paso montecarlo, creceran a su volumen optimo. Por lo tanto, antes de poner el numero es importante considerar el tamano maximo disponoble para las celulas.
# por ejemplo si x= y = 10. El tamano sera 100. Si el volumen maximo de cada celula es 5, solo se podran poner 20 celulas (de preferencia 19)
#t1 y t2 son los dos tipos celulares que se quieren (ACTIVADAS e INHIBIDAS). Si se quiere solo uno de los dos, solo se debe escribir dos veces su nombre, en t1 y t2



		
def coor_xy_tipos(x,y,n,t1,t2):
    
    
	ex = []
	rangox = range(x)
	
	rangoy = range(y)
	
	coor = {}
	num_of_cel = range(n+1)	
	num_of_cel.remove(0)
	# el siguiente loop asigna los diferentes tipos para cada celula
	for i in num_of_cel:
		#aqui se especifican las coordenadas aleatorias.
		a=choice(rangox)
		b=choice(rangoy)
	 # aqui verificamos si el lugar ya esta tomado, si si,se sortean otras coordenadas. 	
		while (a,b) in ex:
			a =choice(rangox)
			b= choice(rangoy)
		tipos=[t1,t2]	
		ch= choice (tipos)
    #llenamos el diccionario con el numero de celula como key, las coordenas y el tipo celular como valores.		
		coor[i]=(a,b,ch)
		ex.append((a,b))
	print coor
	
	#esta ultima parte escribe el piff en un documento. Hay que verificar en que archivo lo escribe, para poder encontrarlo posteriormente.
	
	ar_out = open ("/home/emilio/install_projects/3.7.0/EmilioTesis/adhesionComunication/Simulation/random_2D.pif",'w')
	for k in coor:
		ar_out.write("%s %s %s %s %s %s 0 0\t"%(k,coor[k][2],coor[k][0],coor[k][0],coor[k][1],coor[k][1]))
		ar_out.write("\n")
		
	ar_out.close()
