import numpy as np
from math import *
import linecache
import os
import sys


Chemin_GMSH=sys.argv[1] 
Name_surf=sys.argv[2] 
Fix=float(sys.argv[3]) 
Envergure=float(sys.argv[4])
X_r=float(sys.argv[5])
Rotation=float(sys.argv[6])*3.141592/180.0
File=sys.argv[7]
Ratio=float(sys.argv[8])


# Centre of rotation 
Tab=np.fromfile(File,sep=" ")
N=len(Tab)
for i in range(0,N):
    Tab[i]*=Ratio
norm_max=0
norm_min=1e9
for i in range(0,int(N/2.0)):
    norm_i=Tab[2*i]**2+Tab[2*i+1]**2
    if norm_i>norm_max:
        norm_max=norm_i
        X_max=Tab[2*i]
        Y_max=Tab[2*i+1]
    if norm_i<norm_min:
        norm_min=norm_i
        X_min=Tab[2*i]
        Y_min=Tab[2*i+1]
X_c=X_r*(X_max-X_min) 
Y_c=0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Retrieve the number of faces associated with the surface "Name_surf" and the number of the first face         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
f = open('../Fluid/constant/polyMesh/boundary','r')
contenu = f.read()

Nb_surf=contenu.find(Name_surf)
contenu_2=contenu[Nb_surf:]
Nb_faces=contenu_2.find('nFaces')
contenu_3=contenu_2[Nb_faces:]
Nb_point=contenu_3.find(';')
contenu_4=contenu_3[Nb_point+1:]
Nb_point_2=contenu_4.find(';')
ncar=0
nbFaces=""
car=""
while car!=" ":
	ncar+=1
	car=contenu_3[Nb_point-ncar]
	nbFaces+=car
car=""
ncar=0
nstartFace=""
while car!=" ":
	ncar+=1
	car=contenu_4[Nb_point_2-ncar]
	nstartFace+=car    
nFaces=int("".join(reversed(nbFaces))) # Number of faces
startFace=int("".join(reversed(nstartFace))) # Number of the first face

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                              Rewrite the file 'faces' in the "no-frills" file 'faces2'                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Create():
	k=0
	length=0
	Fich='faces2'
	fichier2=open(Fich,'w')
	fichier=open('../Fluid/constant/polyMesh/faces','r')
	lignes=fichier.readlines()
	for ligne in lignes:
		length=length+1
	for ligne in lignes:
		k+=1
		ligne=ligne.replace('(',' ') # We delete the parentheses
		ligne=ligne.replace(')','')
		if k>20 and k<length:        # We delete the 20 first lines which contain no information
			fichier2.write(ligne)
	fichier2.close
	fichier.close()

print("Retrieval of the number of the faces associated with the surface "+str(Name_surf))
Create()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                        First face to select in faces2                                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Table=np.fromfile("faces2",sep=" ")

ligne=0
k=0
i=0

while 1!=2:
	i=i+k
	if Table[i]>=3 and Table[i]<=8:
		k=int(Table[i])+1
	ligne+=1
	if ligne==startFace+1:
		FirstFace=i 
		break

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                    Storage of the number of the faces in ListFace                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ListFace=np.zeros((nFaces,9))
l=FirstFace+1
ListFace[0,0]=Table[l]
ligne=0
k=0
i=0
er=Table[l-1]
while ligne<nFaces:
	if Table[l]>=3 and Table[l]<=9 and er==0:
		  ligne+=1
		  k=0
		  i+=1
		  er=Table[l]
	else:
		  ListFace[i,k]=Table[l]
		  k+=1
		  er-=1
	l+=1
        if l==len(Table):
            break
os.remove("faces2")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Rewrite the file 'points' in the "no-frills" file 'points2'                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Create2():
	k=0
	length=0
	Fich='points2'
	fichier2=open(Fich,'w')
	fichier=open('../Fluid/constant/polyMesh/points','r')
	lignes=fichier.readlines()
	for ligne in lignes:
		length=length+1
	for ligne in lignes:
		k=k+1
		ligne=ligne.replace('(','') 
		ligne=ligne.replace(')','')
		if k>20 and k<length:
			fichier2.write(ligne)
	fichier2.close
	fichier.close()

print("Retrieval of the coordinates of the points of each face belonging to the surface "+str(Name_surf))
Create2()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                           Writing process of the file maillage.txt                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
print("Gathering information in Maillage.txt")
Maillage='Maillage.txt' # Contains the coordinates of all the points of the surface. Each point corresponds to its face
fichier=open(Maillage,'w')
for i in range(0,int(nFaces)):
	b=""
	for k in range(0,10):
		if int(ListFace[i,k])!=0:
			a=linecache.getline("points2",int(ListFace[i,k])+1)		
			b="{0} {1}".format(i,a)
			fichier.write(b)
		else:
			break

os.remove("points2")
fichier.close()
Table3=np.fromfile("Maillage.txt",count=-1,sep=" ")

g_fix=0
g_fin=0
for i in range(0,int(len(Table3)/4.0)):
	if Table3[4*i-1]<=Fix:
		g_fix+=1
	elif Table3[4*i-1]>=Envergure:
		g_fin+=1

os.remove("Maillage.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                            Writing process of the file Maillage.geo                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
print("Storing points, lines et contours. Can take a little time...")
List_Point=[] # Storage of the points coordinates
List_Line=[] # Storage of the lines without duplicates
List_LineLoop=[] # Storage of the faces contours
List_Face_Fixe=[] # Storage of the information about the fix surfaces (0 --> fix faces in the simulation | 1 --> free boundary condition)
List_Line_Fixe=[]
face_n=0
i=0

while i < int(len(Table3)/4):
	Pt=[]
	Line=[]
        fix=0
	List_fixe=[]
	while face_n==int(Table3[4*i]):
		New_point=[Table3[4*i+1],Table3[4*i+2],Table3[4*i+3]]
                if Table3[4*i+3]>Fix:
                        fix=1
		if New_point in List_Point:
			Pt.append(List_Point.index(New_point))
		else:
			List_Point.append(New_point)
			Pt.append(len(List_Point)-1)
		if g_fix<0.9*g_fin and Table3[4*i+3]<=Fix+5e-4:
			List_fixe.append(Pt[-1])
		i+=1
		if i>=int(len(Table3)/4):
			break
        if i<int(len(Table3)/4):
            face_n=Table3[4*i]
	for j in range(0,len(Pt)):
		if j!=len(Pt)-1:
			New_line=[Pt[j],Pt[j+1]]
		else:
			New_line=[Pt[len(Pt)-1],Pt[0]]
		if New_line in List_Line:
			Line.append(List_Line.index(New_line)+1)
		elif [New_line[1],New_line[0]] in List_Line:
			Line.append("-"+str(List_Line.index([New_line[1],New_line[0]])+1))
		else:
			List_Line.append(New_line)
			Line.append(len(List_Line))
		if g_fix<0.9*g_fin:
			if j!=len(Pt)-1: 
				if (Pt[j] in List_fixe) and (Pt[j+1] in List_fixe):
                                        List_Line_Fixe.append(Line[-1])
			else:
				if (Pt[len(Pt)-1] in List_fixe) and (Pt[0] in List_fixe):
                                        List_Line_Fixe.append(Line[-1])
	Loop=[]
	for k in range(len(Line)-len(Pt),len(Line)):
		Loop.append(Line[k])
	List_LineLoop.append(Loop)
        List_Face_Fixe.append(fix)
	if i>len(Table3)/16 and i<len(Table3)/16+4:
		print("Storing 25%")
	elif i>len(Table3)/8 and i<len(Table3)/8+4:
		print("Storing 50%")
	elif i>3*len(Table3)/16 and i<3*len(Table3)/16+4:
		print("Storing 75%")

print("End of the storing process")	
print("Writing process of the input file for GMSH")	
Mesh='Maillage.geo' 
File=open(Mesh,'w')
for i in range(0,len(List_Point)):
	File.write("Point("+str(i+1)+") = {"+str(List_Point[i][0])+", "+str(List_Point[i][1])+", "+str(List_Point[i][2])+", 1.0};\n")

File.write("\n")

for i in range(0,len(List_Line)):
	File.write("Line("+str(i+1)+") = {"+str(List_Line[i][0]+1)+", "+str(List_Line[i][1]+1)+"};\n")

File.write("\n")

for i in range(0,len(List_LineLoop)):
	File.write("Line Loop("+str(i+1)+") = {")
	for k in range(0,len(List_LineLoop[i])-1):
		File.write(str(List_LineLoop[i][k])+", ")
	File.write(str(List_LineLoop[i][len(List_LineLoop[i])-1])+"};\n")
	File.write("Plane Surface("+str(i+1)+") = {"+str(i+1)+"};\n")

if g_fix<0.9*g_fin:
	File.write("Line Loop("+str(len(List_LineLoop)+1)+") = {")
	for i in range(0, len(List_Line_Fixe)):
		File.write(str(List_Line_Fixe[i])+",")
	File.seek(-1, os.SEEK_END)
	File.truncate()
	File.write("};\n")
	File.write("Plane Surface("+str(len(List_LineLoop)+1)+") = {"+str(len(List_LineLoop)+1)+"};\n")

File.write("\nSurface Loop(1) = {1")

for i in range(1,len(List_LineLoop)):
	File.write(", "+str(i+1))

if g_fix<0.9*g_fin:
	File.write(", "+str(len(List_LineLoop)+1))

File.write('};\n\nVolume(1) = {1};\n\nPhysical Volume("Volume") = {1};\n\nPhysical Surface("Aile") = {')
for i in range(0,len(List_LineLoop)):
    if List_Face_Fixe[i]!=0:
	File.write(str(i+1)+",")
File.seek(-1, os.SEEK_END)
File.truncate()
File.write('};\n\nPhysical Surface("Fix") = {')

if g_fix<0.9*g_fin:
	File.write(str(len(List_LineLoop)+1)+",")
else:
	for i in range(0,len(List_LineLoop)):
		if List_Face_Fixe[i]==0:
		    File.write(str(i+1)+",")
File.seek(-1, os.SEEK_END)
File.truncate()

File.write("};\n\nTransfinite Line {1")
for i in range(1,len(List_Line)):
	File.write(", "+str(i+1))
File.write("} = 1 Using Progression 1;\n\n")

for i in range(0,len(List_LineLoop)):
	File.write("Transfinite Surface {"+str(i+1)+"};\n")

File.close()

File=open('Maillage.geo','r')
lignes=File.readlines()
File2=open('Maillage2.geo','w')
for ligne in lignes:
    if 'Physical Surface' in ligne:
        ligne=ligne.replace('Physical Surface','//Physical Surface')
    File2.write(ligne)
File2.close()
File.close()

print("End of the writing process of Maillage.geo")

print("Meshing process with GMSH")

executable=str(Chemin_GMSH)+" Maillage.geo -3 -parametric -order 2 -o Maillage.inp"
os.system(executable) 
executable2=str(Chemin_GMSH)+" Maillage2.geo -3 -parametric -order 2 -o Maillage2.inp"
os.system(executable2)

