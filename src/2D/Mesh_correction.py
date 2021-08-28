import numpy as np
from numpy import *
import sys
import os
import linecache
import matplotlib.pyplot as plt

# Input parameters
Name_surf=sys.argv[1]
H=float(sys.argv[2])
theta=float(sys.argv[3])*3.141592/180.0
Ratio=float(sys.argv[4])
Material=sys.argv[5]
U_inf=float(sys.argv[6])
y_plus=float(sys.argv[7])
R=float(sys.argv[8])
Taille_maille=float(sys.argv[9])
N=int(sys.argv[10])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Retrieve the number of faces associated with the surface "Name_surf" and the number of the first face         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
f = open('constant/polyMesh/boundary','r')
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
startFace=int("".join(reversed(nstartFace))) # Label of the first face

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                               Rewrite the file 'faces' in the "no-frills" file 'faces2'                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Create():
	k=0
	length=0
	Fich='faces2'
	fichier2=open(Fich,'w')
	fichier=open('constant/polyMesh/faces','r')
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
#                                          First face to select in faces2                                               #
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
#                            Rewrite the file 'points' in the "no-frills" file 'points2'                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Create2():
	k=0
	length=0
	Fich='points2'
	fichier2=open(Fich,'w')
	fichier=open('constant/polyMesh/points','r')
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
#                                      Writing process of the file maillage.txt                                         #
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

TablePts=np.fromfile("points2",count=-1,sep=" ")
os.remove("points2")
fichier.close()
Table3=np.fromfile("Maillage.txt",count=-1,sep=" ")
os.remove("Maillage.txt")

# Point corresponding to the trailing edge
x_pts=[]
x=Table3[1]
for i in range(1,int(len(Table3)/4)):
	if Table3[4*i+1]>x:
		x=Table3[4*i+1]
		y=Table3[4*i+2]
		z=Table3[4*i+3]
x_pts.append([x,y,z])

# Points on the wing that will be modified
for j in range(1,N):
	x=Table3[1]
	for i in range(1,int(len(Table3)/4)):
		test=True
		for l in range(0,len(x_pts)):
			if abs(Table3[4*i+1]-x_pts[l][0])<1e-4 or Table3[4*i+2]>(tan(theta)*(x_pts[0][0]-Table3[4*i+1])+x_pts[0][1]):
				test=False
		if Table3[4*i+1]>x and test==True:
			x=Table3[4*i+1]
			y=Table3[4*i+2]
	x_pts.append([x,y,z])
	x=Table3[1]
	for i in range(1,int(len(Table3)/4)):
		test=True
		for l in range(0,len(x_pts)):
			if abs(Table3[4*i+1]-x_pts[l][0])<1e-4 or Table3[4*i+2]<(tan(theta)*(x_pts[0][0]-Table3[4*i+1])+x_pts[0][1]):
				test=False
		if Table3[4*i+1]>x and test==True:
			x=Table3[4*i+1]
			y=Table3[4*i+2]
	x_pts.append([x,y,z])

# 10 points closest to the selected points on the wing 
List_pts=[]
List_suppr=[]
for n in range(0,N-2):
	x,y=x_pts[n][0],x_pts[n][1]
	List_pts_n=[]
	List_suppr_n=[]
	for j in range(0,10):
		min_=1e5
		for i in range(0,int(len(TablePts)/3)):
		    if i not in List_suppr_n:
		        if abs(TablePts[3*i+2]-z)<1e-5:
		            norme=(TablePts[3*i]-x)**2+(TablePts[3*i+1]-y)**2
		            if norme<min_:
		                min_=norme
		                xx=TablePts[3*i]
		                yy=TablePts[3*i+1]
		                ind=i
		List_pts_n.append([xx,yy])
		List_suppr_n.append(ind)
		x,y=xx,yy
	List_pts.append(List_pts_n)
	List_suppr.append(List_suppr_n)

# Linear regression to approximate the equation of the line formed by the previous 10 points
a=[]
b=[]
for n in range(0,N-2):
	x_moy=y_moy=0
	for i in range(0,10):
		x_moy+=0.1*List_pts[n][i][0]
		y_moy+=0.1*List_pts[n][i][1]
	sigma=0
	var=0
	for i in range(0,10):
		sigma+=(List_pts[n][i][0]-x_moy)*(List_pts[n][i][1]-y_moy)  	
		var+=(List_pts[n][i][0]-x_moy)**2
	a.append(sigma/var)
	b.append(y_moy-sigma*x_moy/var)

#~~~~~~~~~~~~~~~~~~~~~ Features of the boundary layer mesh ~~~~~~~~~~~~~~~~~~~~~#
if Material=='Water':
	rho=1000
	mu=1e-3
elif Material=='Air':
	rho=1
	mu=1e-6
else:
	print('Error in the definition of the material')
	sys.exit()

Re=rho*U_inf*Ratio/mu
Cf=0.026*Re**(-1.0/7.0)
t_wall=0.5*Cf*rho*U_inf**2
u_etoile=sqrt(t_wall/rho)
y=2*mu*y_plus/(rho*u_etoile) # Size of the first cell - The '2*' is for the centre to be at y metre from the wall (OpenFOAM uses collocated grids)

Max=int((log(Taille_maille)-log(y))/log(R)) # Number of cells between the wing and the cell with the same size as the parameter Taille_maille
Length=0 # Thickness of this layer
for n in range(0,Max+1):
	Length+=y*R**n

# Equation of the circle
x1,y1=x_pts[N-1][0],x_pts[N-1][1]
x2,y2=x_pts[N-2][0],x_pts[N-2][1]

R=0.5*sqrt((x1-x2)**2+(y1-y2)**2)
xc,yc=0.5*(x1+x2),0.5*(y1+y2)

u=np.zeros(2)
u[0],u[1]=(x1-x2)/(2*R),(y1-y2)/(2*R)
normal=np.zeros(2)
normal[0],normal[1]=(y1-y2)/(2*R),-(x1-x2)/(2*R)

x_obj=[[xc+R*normal[0],yc+R*normal[1]]]
part=(N-3)/2+1
for n in range(1,part):
	angle=-math.pi*n/(2*part)
	x_obj.append([xc+math.cos(angle)*R*normal[0]+math.sin(angle)*R*u[0],yc+math.cos(angle)*R*normal[1]+math.sin(angle)*R*u[1]])
	angle=math.pi*n/(2*part)
	x_obj.append([xc+math.cos(angle)*R*normal[0]+math.sin(angle)*R*u[0],yc+math.cos(angle)*R*normal[1]+math.sin(angle)*R*u[1]])

delta=[]
for n in range(0,N-2):
	delta.append([x_obj[n][0]-x_pts[n][0],x_obj[n][1]-x_pts[n][1]])

ratio=0.5
List_Pts_modif=[]
List_indice_modif=[]
for n in range(0,N-2):
	List_Pts_modif_n=[]
	List_indice_modif_n=[]
	# Coordinates of the points that belong to the lines 
	for i in range(0,int(len(TablePts)/3)):
		if abs(a[n]*TablePts[3*i]+b[n]-TablePts[3*i+1])/sqrt(1+a[n]**2)<1e-4 and ((x_pts[n][0]-TablePts[3*i])**2+(x_pts[n][1]-TablePts[3*i+1])**2)<0.99*Length**2 and TablePts[3*i]>=x_pts[n][0]:
		    List_Pts_modif_n.append([TablePts[3*i],TablePts[3*i+1],TablePts[3*i+2]])
		    List_indice_modif_n.append(i)
	for j in range(0,len(List_Pts_modif_n)):
		if ((x_pts[n][0]-List_Pts_modif_n[j][0])**2+(x_pts[n][1]-List_Pts_modif_n[j][1])**2)<(Ratio*Length)**2:
			List_Pts_modif_n[j][0]+=delta[n][0]
			List_Pts_modif_n[j][1]+=delta[n][1]
		else:
			d=sqrt((x_pts[n][0]-List_Pts_modif_n[j][0])**2+(x_pts[n][1]-List_Pts_modif_n[j][1])**2)
			List_Pts_modif_n[j][0]+=((1-Ratio-d+Ratio*Length)/(1-Ratio))*delta[n][0]
			List_Pts_modif_n[j][1]+=((1-Ratio-d+Ratio*Length)/(1-Ratio))*delta[n][1]
	List_Pts_modif.append(List_Pts_modif_n)
	List_indice_modif.append(List_indice_modif_n)

k=0
indi=[0 for n in range(0,N-2)]
Fich='constant/polyMesh/points2'
fichier2=open(Fich,'w')
fichier=open('constant/polyMesh/points','r')
lignes=fichier.readlines()

for ligne in lignes:
	k+=1
	for n in range(0,N-2):
		if (k-21) in List_indice_modif[n]:
			ligne="("+str(List_Pts_modif[n][indi[n]][0])+" "+str(List_Pts_modif[n][indi[n]][1])+" "+str(List_Pts_modif[n][indi[n]][2])+")\n"
			indi[n]+=1
	fichier2.write(ligne)
fichier2.close
fichier.close()
os.remove('constant/polyMesh/points')
os.rename('constant/polyMesh/points2','constant/polyMesh/points')
