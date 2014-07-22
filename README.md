CodigoAdhesionComunicacion
==========================

Aquí está el código por si es necesario descargarlo. Hay seis archivos: 

1) El piff_aleatorio.py crea archivos pif (equivalentes a txt) para el estappable de inicio definido en adhesionComunicacion.xml, con células dispersas para los escenarios de agregacion y poblaciones celulares

2) El adhesionComunicacion.xml que genera las condiciones iniciales del programa

3) El adhesionComunicacionSteppables.py que define todas las clases

4) El adhesionComunicacion.py que corre todo el programa junto. 

5) El ComunicacionFinal.xml en donde están definidas las EDO bajo forma de SBML.

6) Un ejemplo del resultado del piff_aleatorio, denominado random2D.txt. Si el piffaleatorio tiene problemas, cargar directamente el random2D.txt en la misma carpeta que el resto.

7) El archivo adhesionComunicacion.cc3d. Este es el que es llamado al usar el programa CC3D (compucell3d.org) 

Para el correcto funcionamiento del programa, es importante generar las siguiente carpetas con los archivos adecuados:

     1. Carpeta X
     
        1.1 adhesionComunicacion.cc3d
     
        1.2 Simulation (carpeta)
     
            1.2.1 adhesionComunicacion.xml
     
            1.2.2 adhesionComunicacion.py
     
            1.2.3 adhesionComunicacionSteppables.py
     
            1.2.4 ComunicacionFinal.xml
     
            1.2.5 random_2D.txt
     
            1.2.6 piff_aleatorio.py (este es opcional).
     
     
Cualquier cosa que no funcione es probable que tenga que ver con la transcripción que hice de algunos términos. Cualquier duda avisenme y lo corrigo emiliomora92@gmail.com

Chido. 
           *
