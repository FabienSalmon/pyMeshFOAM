import numpy as np
from math import *
import linecache
import time
import os
import sys

# Rewriting the inp file to make it match with the input format for CalculiX

INPfile=open('Maillage.inp','r')
Mesh=INPfile.read()

Line_Aile=Mesh.find('ELSET=Aile')
Mesh2=Mesh[Line_Aile+10:]
k=0
Face_number_Aile=[]

while Mesh2[k]!='E':
    face=""
    while Mesh2[k]!=',':
        face+=Mesh2[k] 
        k+=1
    try:
        Face_number_Aile.append(int(face))
    except:
        break
    k+=2

Node_number_Aile=[]
Line_elem=Mesh.find('*******')
Mesh3=Mesh[Line_elem:]
for i in range(0,len(Face_number_Aile)):
    line_i=Mesh3.find(str(Face_number_Aile[i]))
    Mesh_i=Mesh3[line_i:]
    j=0
    while Mesh_i[j]!=',':
        j+=1
    j+=2
    n=0
    while Mesh_i[j+n]!='E':
        node=""
        while Mesh_i[j+n]!=',' and Mesh_i[j+n]!='*':
            node+=Mesh_i[j+n]
            n+=1
        try:
            if int(node) in Node_number_Aile:
                pass
            else:
                Node_number_Aile.append(int(node))
        except:
            break
        n+=2

Node_number_Aile.sort()

Line_fix=Mesh.find('ELSET=Fix')
Mesh2=Mesh[Line_fix+10:]
k=0
Face_number_Fix=[]

while Mesh2[k]!='E':
    face=""
    while Mesh2[k]!=',':
        face+=Mesh2[k] 
        k+=1
    try:
        Face_number_Fix.append(int(face))
    except:
        break
    k+=2

Node_number_Fix=[]
Line_elem=Mesh.find('*******')
Mesh3=Mesh[Line_elem:]
for i in range(0,len(Face_number_Fix)):
    line_i=Mesh3.find(str(Face_number_Fix[i]))
    Mesh_i=Mesh3[line_i:]
    j=0
    while Mesh_i[j]!=',':
        j+=1
    j+=2
    n=0
    while Mesh_i[j+n]!='E':
        node=""
        while Mesh_i[j+n]!=',' and Mesh_i[j+n]!='*':
            node+=Mesh_i[j+n]
            n+=1
        try:
            if int(node) in Node_number_Fix:
                pass
            else:
                Node_number_Fix.append(int(node))
        except:
            break
        n+=2

Node_number_Fix.sort()
INPfile.close()

INPfile2=open('Maillage2.inp','r')
Mesh=INPfile2.read()

Line_first_del=Mesh.find('*ELSET')

Fichier=open('All.msh','w')
for m in range(0,len(Mesh)):
    if m>=Line_first_del:
        break
    else: #if m>Line_fin_elem or m<Line_elem:
        Fichier.write(Mesh[m])

INPfile2.close()
Fichier.close()

Fichier=open('fix.nam','w')
Fichier.write("*NSET,NSET=Nfix\n")
line_write=""
for i in range(0,len(Node_number_Fix)):
    line_write+=str(Node_number_Fix[i])+", "
    if len(line_write)>50:
        Fichier.write(str(line_write)+"\n")
        line_write=""
Fichier.close()

Fichier=open('surface.nam','w')
Fichier.write("*NSET,NSET=Nsurface\n")
line_write=""
for i in range(0,len(Node_number_Aile)):
    line_write+=str(Node_number_Aile[i])+", "
    if len(line_write)>50:
        Fichier.write(str(line_write)+"\n")
        line_write=""
Fichier.close()
