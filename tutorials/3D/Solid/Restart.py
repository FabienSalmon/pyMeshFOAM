import numpy as np
import sys
import os
from math import *

File=sys.argv[1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Premiere partie : calcul de la condition limite en rotation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
Nom_frd=str(File)+'.frd'
fichier=open(str(Nom_frd),'r')
Nom_inp=str(File)+'.inp'
fichier2=open(str(Nom_inp),'r')
fichier3=open('Intermed.inp','w')

lines=fichier2.readlines()
n=0
for line in lines:
	n+=1
	if '*CLOAD' in line:
		break
	else:
		fichier3.write(line)

lignes=fichier.readlines()
k=0
for ligne in lignes:
	k+=1
	if "DISP" in ligne:
		Line_disp=k
	if "VELO" in ligne:
		Line_velo=k

j=0
disp=False
velo=False
for ligne in lignes:
	j+=1
	if len(ligne)<10:
		disp=False
		velo=False
	if j==Line_disp+5 and disp==False:
		fichier3.write('*INITIAL CONDITIONS, TYPE=DISPLACEMENT\n')
	if j==Line_disp+5 or disp==True:
		disp=True
		l=0
		node,disp_x,disp_y,disp_z="","","",""
		while ligne[4+l]==' ':
			l+=1
		while ligne[4+l]!=' ' and (ligne[3+l]=='E' or ligne[4+l]!='-'):
			node+=str(ligne[4+l])
			l+=1
		if ligne[4+l]=='-':
			disp_x+="-"
		l+=1
		while ligne[4+l]!=' ' and (ligne[3+l]=='E' or ligne[4+l]!='-'):
			disp_x+=str(ligne[4+l])
			l+=1
		if ligne[4+l]=='-':
			disp_y+="-"
		l+=1
		while ligne[4+l]!=' ' and (ligne[3+l]=='E' or ligne[4+l]!='-'):
			disp_y+=str(ligne[4+l])
			l+=1
		while 4+l<len(ligne):
			disp_z+=str(ligne[4+l])
			l+=1
		fichier3.write(str(node)+',1,'+str(disp_x)+"\n"+str(node)+',2,'+str(disp_y)+"\n"+str(node)+',3,'+str(disp_z))
	if j==Line_velo+5 and velo==False:
		fichier3.write('\n*INITIAL CONDITIONS, TYPE=VELOCITY\n')
	if j==Line_velo+5 or velo==True:
		velo=True
		m=0
		node,vel_x,vel_y,vel_z="","","",""
		while ligne[4+m]==' ':
			m+=1
		while ligne[4+m]!=' ' and (ligne[3+m]=='E' or ligne[4+m]!='-'):
			node+=str(ligne[4+m])
			m+=1
		if ligne[4+m]=='-':
			vel_x+="-"
		m+=1
		while ligne[4+m]!=' ' and (ligne[3+m]=='E' or ligne[4+m]!='-'):
			vel_x+=str(ligne[4+m])
			m+=1
		if ligne[4+m]=='-':
			vel_y+="-"
		m+=1
		while ligne[4+m]!=' ' and (ligne[3+m]=='E' or ligne[4+m]!='-'):
			vel_y+=str(ligne[4+m])
			m+=1
		while 4+m<len(ligne):
			vel_z+=str(ligne[4+m])
			m+=1
		fichier3.write(str(node)+',1,'+str(vel_x)+"\n"+str(node)+',2,'+str(vel_y)+"\n"+str(node)+',3,'+str(vel_z))

p=0
for line in lines:
	p+=1
	if p>=n:
		fichier3.write(line)

fichier.close()
fichier2.close()
fichier3.close()
os.remove(str(Nom_inp))
os.rename('Intermed.inp',str(Nom_inp))
