import numpy as np
from numpy import *
import sys
import os
import matplotlib.pyplot as plt

# Input parameters
twoD=sys.argv[1]
Longueur=float(sys.argv[2])
Ratio=float(sys.argv[3])
X_r=float(sys.argv[4])
Rotation=float(sys.argv[5])*3.141592/180.0
File=sys.argv[6]
Material=sys.argv[7]
U_inf=float(sys.argv[8])
mult_chordx=float(sys.argv[9])
mult_chordy=float(sys.argv[10])
mult_span=float(sys.argv[11])
mult_chordH=float(sys.argv[12])
y_plus=float(sys.argv[13])
R=float(sys.argv[14])
Taille_maille=float(sys.argv[15])
Nz=int(sys.argv[16])
xx=int(sys.argv[17])

N_ext=50 # Number of cells in the direction of the normal to the wing in the external layer (far away from the wing)
ratio_loin=2 # Ratio of the successive cell sizes in the external layer

if Material=='Water':
	rho=1000
	mu=1e-3
elif Material=='Air':
	rho=1
	mu=1e-6
else:
	print('Error in the definition of the material')
	sys.exit()

# Writing process of the input file in case of a NACA00xx profile
if xx!=0:
    file_=open(File,'w')
    N0=10
    N1=19
    N2=80
    for k in range(0,N0):
        file_.write(str(k*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*0.0001)-0.126*k*0.0001-0.3516*(k*0.0001)**2+0.2843*(k*0.0001)**3-0.1015*(k*0.0001)**4))+"\n")
    for k in range(0,N1):
        file_.write(str(k*0.001+N0*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*0.001+N0*0.0001)-0.126*(k*0.001+N0*0.0001)-0.3516*(k*0.001+N0*0.0001)**2+0.2843*(k*0.001+N0*0.0001)**3-0.1015*(k*0.001+N0*0.0001)**4))+"\n")
    for k in range(0,N2):
        file_.write(str(k*0.01+N1*0.001+N0*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*0.01+N1*0.001+N0*0.0001)-0.126*(k*0.01+N1*0.001+N0*0.0001)-0.3516*(k*0.01+N1*0.001+N0*0.0001)**2+0.2843*(k*0.01+N1*0.001+N0*0.0001)**3-0.1015*(k*0.01+N1*0.001+N0*0.0001)**4))+"\n")
    k=0
    while k*0.0005+N2*0.01+N1*0.001+N0*0.0001<=1:
        file_.write(str(k*0.0005+N2*0.01+N1*0.001+N0*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*0.0005+N2*0.01+N1*0.001+N0*0.0001)-0.126*(k*0.0005+N2*0.01+N1*0.001+N0*0.0001)-0.3516*(k*0.0005+N2*0.01+N1*0.001+N0*0.0001)**2+0.2843*(k*0.0005+N2*0.01+N1*0.001+N0*0.0001)**3-0.1015*(k*0.0005+N2*0.01+N1*0.001+N0*0.0001)**4))+"\n")
        k+=1
    file_.close()

Table=np.fromfile(File,sep=" ")
N=len(Table)

# Homothety
for i in range(0,N):
    Table[i]*=Ratio

# Centre of rotation
norm_max=0
norm_min=1e9
for i in range(0,int(N/2.0)):
    norm_i=Table[2*i]**2+Table[2*i+1]**2
    if norm_i>norm_max:
        norm_max=norm_i
        X_max=Table[2*i]
        Y_max=Table[2*i+1]
    if norm_i<norm_min:
        norm_min=norm_i
        X_min=Table[2*i]
        Y_min=Table[2*i+1]
X_c=X_min+X_r*(X_max-X_min) 
Y_c=0

# Chord
corde=X_max-X_min

H=mult_chordH*corde # Thickness of the intermediate layer around the wing 
alpha_longueur=mult_chordx*corde
alpha_largeur=mult_chordy*corde

#~~~~~~~~~~~~~~~~~~~~~ Features of the mesh of the boundary layer ~~~~~~~~~~~~~~~~~~~~~#
Re=rho*U_inf*corde/mu
Cf=0.026*Re**(-1.0/7.0)
t_wall=0.5*Cf*rho*U_inf**2
u_etoile=sqrt(t_wall/rho)
y=2*mu*y_plus/(rho*u_etoile) # Size of the first cell - The '2*' is for the centre to be at y metre from the wall (OpenFOAM uses collocated grids)

Max=int((log(Taille_maille)-log(y))/log(R)) # Number of cells between the wing and the cell with the same size as the parameter Taille_maille
Length=0 # Thickness of this layer 
for n in range(0,Max+1):
	Length+=y*R**n 

# Storage of the points of the profile in the vector Profile
Profile=np.zeros((2*int(N/2.0),2))
for i in range(0,int(N/2.0)):
	Profile[i,0],Profile[i,1]=Table[2*i],Table[2*i+1]

# Barycentre
X_g=0.5*(max([Profile[i,0] for i in range(0,int(N/2.0))])-min([Profile[i,0] for i in range(0,int(N/2.0))]))
Y_g=0

# Writing the extended contours
M=9
Contour1=np.zeros((2*(int(N/2.0)+M)-1,2))
Contour2=np.zeros((2*(int(N/2.0)+M)-1,2))
Contour1[0,0],Contour1[0,1]=-Length,0
Contour2[0,0],Contour2[0,1]=-H,0

for i in range(1,int(N/2.0)):
	x_M1=Profile[i-1,0]
	x_M2=Profile[i,0]
	y_M1=Profile[i-1,1]
	y_M2=Profile[i,1]
	norme=sqrt((y_M1-y_M2)**2+(x_M1-x_M2)**2)
	Contour1[i,0],Contour1[i,1]=0.5*(x_M1+x_M2)+Length*(y_M1-y_M2)/norme,0.5*(y_M1+y_M2)+Length*(x_M2-x_M1)/norme
	Contour2[i,0],Contour2[i,1]=0.5*(x_M1+x_M2)+H*(y_M1-y_M2)/norme,0.5*(y_M1+y_M2)+H*(x_M2-x_M1)/norme

# Adding a circle to close the extended contours
R1=abs(Y_max-Contour1[int(N/2.0)-1,1])
R2=abs(Y_max-Contour2[int(N/2.0)-1,1])

for i in range(2,M+1):
	Contour1[int(N/2.0)-2+i,0],Contour1[int(N/2.0)-2+i,1]=X_max+i*R1/(M+1),sqrt(abs(R1**2-(i*R1/(M+1))**2))+Y_max
	Contour2[int(N/2.0)-2+i,0],Contour2[int(N/2.0)-2+i,1]=X_max+i*R2/(M+1),sqrt(abs(R2**2-(i*R2/(M+1))**2))+Y_max

# Symmetry with respect to the axis y=0 (the profile is assumed to be symmetrical)
for i in range(0,int(N/2.0)+M-1):
	Contour1[i+int(N/2.0)+M-1,0],Contour1[i+int(N/2.0)+M-1,1]=Contour1[i,0],-Contour1[i,1]
	Contour2[i+int(N/2.0)+M-1,0],Contour2[i+int(N/2.0)+M-1,1]=Contour2[i,0],-Contour2[i,1]
	if i<int(N/2.0):
		Profile[i+int(N/2.0),0],Profile[i+int(N/2.0),1]=Profile[i,0],-Profile[i,1]

Contour1[2*(int(N/2.0)+M)-2,0],Contour1[2*(int(N/2.0)+M)-2,1]=X_max+Length,0
Contour2[2*(int(N/2.0)+M)-2,0],Contour2[2*(int(N/2.0)+M)-2,1]=X_max+H,0

L=6															
List_profile=np.zeros((L,2))
List_etendu1=np.zeros((L,2))
List_etendu2=np.zeros((L,2))
List_bloc=np.zeros((14,2))												
indice_profile,indice_etendu1,indice_etendu2=[],[],[]

# We look for the point that minimizes the distance with the line X=Xg+-D1_2*chord/2
D1=0.75
D2=0.66
for k in range(0,2):
	norm_profile=1e5
	for i in range(0,int(N/2.0)):
		norm_profile_test=(Profile[i,0]-(X_g+(-1)**k*D1*corde/2))**2
		if norm_profile_test<norm_profile:
			norm_profile=norm_profile_test
			indice=i
	indice_profile.append(indice)
	List_profile[k,0],List_profile[k,1]=Profile[indice_profile[k],0],Profile[indice_profile[k],1]
indice_profile.append(int(N/2)) 											

for k in range(0,2):
	norm_etendu1=1e5
	for i in range(0,int(N/2.0)+M-1):
		norm_etendu1_test=(Contour1[i,0]-(X_g+(-1)**k*D1*corde/2))**2
		if norm_etendu1_test<norm_etendu1:
			norm_etendu1=norm_etendu1_test
			indice=i
	indice_etendu1.append(indice)
	List_etendu1[k,0],List_etendu1[k,1]=Contour1[indice_etendu1[k],0],Contour1[indice_etendu1[k],1]
indice_etendu1.append(int(N/2)+M) 											

for k in range(0,2):
	norm_etendu2=1e5
	for i in range(0,int(N/2.0)+M-1):
		norm_etendu2_test=(Contour2[i,0]-(X_g+(-1)**k*D2*(corde+2*H)/2))**2
		if norm_etendu2_test<norm_etendu2:
			norm_etendu2=norm_etendu2_test
			indice=i
	indice_etendu2.append(indice)
	List_etendu2[k,0],List_etendu2[k,1]=Contour2[indice_etendu2[k],0],Contour2[indice_etendu2[k],1]
indice_etendu2.append(int(N/2)+M) 											

List_profile[2,0],List_profile[2,1]=List_profile[1,0],-List_profile[1,1]
List_profile[3,0],List_profile[3,1]=List_profile[0,0],-List_profile[0,1]
List_profile[4,0],List_profile[4,1]=corde,0 										
List_profile[5,0],List_profile[5,1]=0,0 										
List_etendu1[2,0],List_etendu1[2,1]=List_etendu1[1,0],-List_etendu1[1,1]
List_etendu1[3,0],List_etendu1[3,1]=List_etendu1[0,0],-List_etendu1[0,1]
List_etendu1[4,0],List_etendu1[4,1]=corde+Length,0									
List_etendu1[5,0],List_etendu1[5,1]=-Length,0										
List_etendu2[2,0],List_etendu2[2,1]=List_etendu2[1,0],-List_etendu2[1,1]
List_etendu2[3,0],List_etendu2[3,1]=List_etendu2[0,0],-List_etendu2[0,1]
List_etendu2[4,0],List_etendu2[4,1]=corde+H,0										
List_etendu2[5,0],List_etendu2[5,1]=-H,0										

# Rotation
for i in range(0,2*(int(N/2.0)+M)-1):
	if i<2*int(N/2.0):
		X_profile,Y_profile=Profile[i,0]-X_c,Profile[i,1]-Y_c
		Profile[i,0],Profile[i,1]=X_c+cos(Rotation)*X_profile+sin(Rotation)*Y_profile,Y_c-sin(Rotation)*X_profile+cos(Rotation)*Y_profile
	X_contour1,Y_contour1=Contour1[i,0]-X_c,Contour1[i,1]-Y_c
	Contour1[i,0],Contour1[i,1]=X_c+cos(Rotation)*X_contour1+sin(Rotation)*Y_contour1,Y_c-sin(Rotation)*X_contour1+cos(Rotation)*Y_contour1
	X_contour2,Y_contour2=Contour2[i,0]-X_c,Contour2[i,1]-Y_c
	Contour2[i,0],Contour2[i,1]=X_c+cos(Rotation)*X_contour2+sin(Rotation)*Y_contour2,Y_c-sin(Rotation)*X_contour2+cos(Rotation)*Y_contour2

for j in range(0,L):
	X_profile,Y_profile=List_profile[j,0]-X_c,List_profile[j,1]-Y_c
	List_profile[j,0],List_profile[j,1]=X_c+cos(Rotation)*X_profile+sin(Rotation)*Y_profile,Y_c-sin(Rotation)*X_profile+cos(Rotation)*Y_profile
	X_contour1,Y_contour1=List_etendu1[j,0]-X_c,List_etendu1[j,1]-Y_c
	List_etendu1[j,0],List_etendu1[j,1]=X_c+cos(Rotation)*X_contour1+sin(Rotation)*Y_contour1,Y_c-sin(Rotation)*X_contour1+cos(Rotation)*Y_contour1
	X_contour2,Y_contour2=List_etendu2[j,0]-X_c,List_etendu2[j,1]-Y_c
	List_etendu2[j,0],List_etendu2[j,1]=X_c+cos(Rotation)*X_contour2+sin(Rotation)*Y_contour2,Y_c-sin(Rotation)*X_contour2+cos(Rotation)*Y_contour2

# Adding points on the block
List_bloc[0,0],List_bloc[0,1]=alpha_longueur+X_g,List_etendu2[4,1]
List_bloc[1,0],List_bloc[1,1]=alpha_longueur+X_g,List_etendu2[0,1]
List_bloc[2,0],List_bloc[2,1]=List_etendu2[0,0],alpha_largeur+Y_g
List_bloc[3,0],List_bloc[3,1]=List_etendu2[1,0],alpha_largeur+Y_g
List_bloc[4,0],List_bloc[4,1]=-alpha_longueur+X_g,List_etendu2[1,1]
List_bloc[5,0],List_bloc[5,1]=-alpha_longueur+X_g,List_etendu2[5,1]
List_bloc[6,0],List_bloc[6,1]=-alpha_longueur+X_g,List_etendu2[2,1]
List_bloc[7,0],List_bloc[7,1]=List_etendu2[2,0],-alpha_largeur+Y_g
List_bloc[8,0],List_bloc[8,1]=List_etendu2[3,0],-alpha_largeur+Y_g
List_bloc[9,0],List_bloc[9,1]=alpha_longueur+X_g,List_etendu2[3,1]
List_bloc[10,0],List_bloc[10,1]=alpha_longueur+X_g,alpha_largeur+Y_g
List_bloc[11,0],List_bloc[11,1]=-alpha_longueur+X_g,alpha_largeur+Y_g
List_bloc[12,0],List_bloc[12,1]=-alpha_longueur+X_g,-alpha_largeur+Y_g
List_bloc[13,0],List_bloc[13,1]=alpha_longueur+X_g,-alpha_largeur+Y_g

# List of points that will be in the splines
List_profile_pts_edge0=[]
List_profile_pts_edge1=[]
List_profile_pts_edge2=[]
List_profile_pts_edge3=[]
List_profile_pts_edge4=[]
List_profile_pts_edge5=[]
List_contour1_pts_edge0=[]
List_contour1_pts_edge1=[]
List_contour1_pts_edge2=[]
List_contour1_pts_edge3=[]
List_contour1_pts_edge4=[]
List_contour1_pts_edge5=[]
List_contour2_pts_edge0=[]
List_contour2_pts_edge1=[]
List_contour2_pts_edge2=[]
List_contour2_pts_edge3=[]
List_contour2_pts_edge4=[]
List_contour2_pts_edge5=[]

for i in range(0,int(N/2.0)):
	if i>indice_profile[0]:
		List_profile_pts_edge0.append([Profile[i,0],Profile[i,1]])
	elif i>indice_profile[1]:
		List_profile_pts_edge1.append([Profile[i,0],Profile[i,1]])
	else:
		List_profile_pts_edge2.append([Profile[i,0],Profile[i,1]])
for i in range(2*int(N/2.0)-1,int(N/2.0)-1,-1):
	if i>indice_profile[0]+int(N/2.0):
		List_profile_pts_edge5.append([Profile[i,0],Profile[i,1]])
	elif i>indice_profile[1]+int(N/2.0):
		List_profile_pts_edge4.append([Profile[i,0],Profile[i,1]])
	else: 
		List_profile_pts_edge3.append([Profile[i,0],Profile[i,1]])

for i in range(0,int(N/2)+M-1):
	if i>indice_etendu1[0]:
		List_contour1_pts_edge0.append([Contour1[i,0],Contour1[i,1]])
	elif i>indice_etendu1[1]:
		List_contour1_pts_edge1.append([Contour1[i,0],Contour1[i,1]])
	else:
		List_contour1_pts_edge2.append([Contour1[i,0],Contour1[i,1]])
for i in range(2*(int(N/2.0)+M)-2,int(N/2)+M-1,-1):
	if i>indice_etendu1[0]+int(N/2)+M-1:
		List_contour1_pts_edge5.append([Contour1[i,0],Contour1[i,1]])
	elif i>indice_etendu1[1]+int(N/2)+M-1:
		List_contour1_pts_edge4.append([Contour1[i,0],Contour1[i,1]])
	else:
		List_contour1_pts_edge3.append([Contour1[i,0],Contour1[i,1]])

for i in range(0,int(N/2)+M-1):
	if i>indice_etendu2[0]:
		List_contour2_pts_edge0.append([Contour2[i,0],Contour2[i,1]])
	elif i>indice_etendu2[1]:
		List_contour2_pts_edge1.append([Contour2[i,0],Contour2[i,1]])
	else:
		List_contour2_pts_edge2.append([Contour2[i,0],Contour2[i,1]])
for i in range(2*(int(N/2.0)+M)-2,int(N/2)+M-1,-1):
	if i>indice_etendu2[0]+int(N/2)+M-1:
		List_contour2_pts_edge5.append([Contour2[i,0],Contour2[i,1]])
	elif i>indice_etendu2[1]+int(N/2)+M-1:
		List_contour2_pts_edge4.append([Contour2[i,0],Contour2[i,1]])
	else:
		List_contour2_pts_edge3.append([Contour2[i,0],Contour2[i,1]])

# Calculation of the length of the edges
len0=len1=len2=len3=len4=len5=0

for i in range(1,len(List_contour1_pts_edge0)):
	len0+=sqrt((List_contour1_pts_edge0[i][0]-List_contour1_pts_edge0[i-1][0])**2+(List_contour1_pts_edge0[i][1]-List_contour1_pts_edge0[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge1)):
	len1+=sqrt((List_contour1_pts_edge1[i][0]-List_contour1_pts_edge1[i-1][0])**2+(List_contour1_pts_edge1[i][1]-List_contour1_pts_edge1[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge2)):
	len2+=sqrt((List_contour1_pts_edge2[i][0]-List_contour1_pts_edge2[i-1][0])**2+(List_contour1_pts_edge2[i][1]-List_contour1_pts_edge2[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge3)):
	len3+=sqrt((List_contour1_pts_edge3[i][0]-List_contour1_pts_edge3[i-1][0])**2+(List_contour1_pts_edge3[i][1]-List_contour1_pts_edge3[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge4)):
	len4+=sqrt((List_contour1_pts_edge4[i][0]-List_contour1_pts_edge4[i-1][0])**2+(List_contour1_pts_edge4[i][1]-List_contour1_pts_edge4[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge5)):
	len5+=sqrt((List_contour1_pts_edge5[i][0]-List_contour1_pts_edge5[i-1][0])**2+(List_contour1_pts_edge5[i][1]-List_contour1_pts_edge5[i-1][1])**2)

# Writing the blockMeshDict file
file=open('system/blockMeshDict','w')
file.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      blockMeshDict;\n}\n\n")
file.write("convertToMeters 1;\n\nvertices\n(\n")

file.write("    ("+str(List_profile[4,0])+" "+str(List_profile[4,1])+" 0)   //"+str(0)+"\n")
for i in range(1,3):
	file.write("    ("+str(List_profile[i-1,0])+" "+str(List_profile[i-1,1])+" 0)   //"+str(i)+"\n")
file.write("    ("+str(List_profile[5,0])+" "+str(List_profile[5,1])+" 0)   //"+str(3)+"\n")
for i in range(4,L):
	file.write("    ("+str(List_profile[i-2,0])+" "+str(List_profile[i-2,1])+" 0)   //"+str(i)+"\n")

file.write("    ("+str(List_etendu1[4,0])+" "+str(List_etendu1[4,1])+" 0)   //"+str(L)+"\n")
for i in range(1,3):
	file.write("    ("+str(List_etendu1[i-1,0])+" "+str(List_etendu1[i-1,1])+" 0)   //"+str(i+L)+"\n")
file.write("    ("+str(List_etendu1[5,0])+" "+str(List_etendu1[5,1])+" 0)   //"+str(3+L)+"\n")
for i in range(4,L):
	file.write("    ("+str(List_etendu1[i-2,0])+" "+str(List_etendu1[i-2,1])+" 0)   //"+str(i+L)+"\n")

file.write("    ("+str(List_etendu2[4,0])+" "+str(List_etendu2[4,1])+" 0)   //"+str(2*L)+"\n")
for i in range(1,3):
	file.write("    ("+str(List_etendu2[i-1,0])+" "+str(List_etendu2[i-1,1])+" 0)   //"+str(i+2*L)+"\n")
file.write("    ("+str(List_etendu2[5,0])+" "+str(List_etendu2[5,1])+" 0)   //"+str(3+2*L)+"\n")
for i in range(4,L):
	file.write("    ("+str(List_etendu2[i-2,0])+" "+str(List_etendu2[i-2,1])+" 0)   //"+str(i+2*L)+"\n")

for i in range(0,14):
	file.write("    ("+str(List_bloc[i,0])+" "+str(List_bloc[i,1])+" 0)   //"+str(i+3*L)+"\n")

file.write("    ("+str(List_profile[4,0])+" "+str(List_profile[4,1])+" "+str(Longueur)+")   //"+str(3*L+14)+"\n")
for i in range(1,3):
	file.write("    ("+str(List_profile[i-1,0])+" "+str(List_profile[i-1,1])+" "+str(Longueur)+")   //"+str(i+3*L+14)+"\n")
file.write("    ("+str(List_profile[5,0])+" "+str(List_profile[5,1])+" "+str(Longueur)+")   //"+str(3+3*L+14)+"\n")
for i in range(4,L):
	file.write("    ("+str(List_profile[i-2,0])+" "+str(List_profile[i-2,1])+" "+str(Longueur)+")   //"+str(i+3*L+14)+"\n")

file.write("    ("+str(List_etendu1[4,0])+" "+str(List_etendu1[4,1])+" "+str(Longueur)+")   //"+str(4*L+14)+"\n")
for i in range(1,3):
	file.write("    ("+str(List_etendu1[i-1,0])+" "+str(List_etendu1[i-1,1])+" "+str(Longueur)+")   //"+str(i+4*L+14)+"\n")
file.write("    ("+str(List_etendu1[5,0])+" "+str(List_etendu1[5,1])+" "+str(Longueur)+")   //"+str(3+i+4*L+14)+"\n")
for i in range(4,L):
	file.write("    ("+str(List_etendu1[i-2,0])+" "+str(List_etendu1[i-2,1])+" "+str(Longueur)+")   //"+str(i+4*L+14)+"\n")

file.write("    ("+str(List_etendu2[4,0])+" "+str(List_etendu2[4,1])+" "+str(Longueur)+")   //"+str(5*L+14)+"\n")
for i in range(1,3):
	file.write("    ("+str(List_etendu2[i-1,0])+" "+str(List_etendu2[i-1,1])+" "+str(Longueur)+")   //"+str(i+5*L+14)+"\n")
file.write("    ("+str(List_etendu2[5,0])+" "+str(List_etendu2[5,1])+" "+str(Longueur)+")   //"+str(3+i+5*L+14)+"\n")
for i in range(4,L):
	file.write("    ("+str(List_etendu2[i-2,0])+" "+str(List_etendu2[i-2,1])+" "+str(Longueur)+")   //"+str(i+5*L+14)+"\n")

for i in range(0,14):
	file.write("    ("+str(List_bloc[i,0])+" "+str(List_bloc[i,1])+" "+str(Longueur)+")   //"+str(i+6*L+14)+"\n")

Ny0=int(len0/Taille_maille) # Number of cells in the normal to the wing direction in the zone 0
Ny1=int(len1/Taille_maille) # Number of cells in the normal to the wing direction in the zone 1
Ny2=int(len2/Taille_maille) # Number of cells in the normal to the wing direction in the zone 2
Ny3=int(len3/Taille_maille) # Number of cells in the normal to the wing direction in the zone 3
Ny4=int(len4/Taille_maille) # Number of cells in the normal to the wing direction in the zone 4
Ny5=Ny0 # Number of cells in the normal to the wing direction in the zone 5
N_int=int(H/Taille_maille) # Number of cells in the normal to the wing direction in the intermediate layer

ratio_aile=R**Max 

file.write(");\n\nblocks\n(\n")
file.write("    hex (1 0 6 7 33 32 38 39) ("+str(Ny0)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #0 
file.write("    hex (2 1 7 8 34 33 39 40) ("+str(Ny1)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #1
file.write("    hex (3 2 8 9 35 34 40 41) ("+str(Ny2)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #2
file.write("    hex (4 3 9 10 36 35 41 42) ("+str(Ny3)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #3
file.write("    hex (5 4 10 11 37 36 42 43) ("+str(Ny4)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #4
file.write("    hex (0 5 11 6 32 37 43 38) ("+str(Ny5)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #5
file.write("    hex (7 6 12 13 39 38 44 45) ("+str(Ny0)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #6
file.write("    hex (8 7 13 14 40 39 45 46) ("+str(Ny1)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #7
file.write("    hex (9 8 14 15 41 40 46 47) ("+str(Ny2)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #8
file.write("    hex (10 9 15 16 42 41 47 48) ("+str(Ny3)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #9
file.write("    hex (11 10 16 17 43 42 48 49) ("+str(Ny4)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #10
file.write("    hex (6 11 17 12 38 43 49 44) ("+str(Ny5)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #11
file.write("    hex (13 12 18 19 45 44 50 51) ("+str(Ny0)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_loin)+" 1)\n") #12
file.write("    hex (13 19 28 20 45 51 60 52) ("+str(N_ext)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(ratio_loin)+" "+str(ratio_loin)+" 1)\n") #13
file.write("    hex (14 13 20 21 46 45 52 53) ("+str(Ny1)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_loin)+" 1)\n") #14
file.write("    hex (22 14 21 29 54 46 53 61) ("+str(N_ext)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(ratio_loin)+" "+str(ratio_loin)+" 1)\n") #15
file.write("    hex (15 14 22 23 47 46 54 55) ("+str(Ny2)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(1.0/ratio_loin)+" 1)\n") #16
file.write("    hex (16 15 23 24 48 47 55 56) ("+str(Ny3)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(1.0/ratio_loin)+" 1)\n") #17
file.write("    hex (30 25 16 24 62 57 48 56) ("+str(N_ext)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(ratio_loin)+" "+str(ratio_loin)+" 1)\n") #18
file.write("    hex (17 16 25 26 49 48 57 58) ("+str(Ny4)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(1.0/ratio_loin)+" 1)\n") #19
file.write("    hex (26 31 27 17 58 63 59 49) ("+str(N_ext)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(ratio_loin)+" "+str(ratio_loin)+" 1)\n") #20
file.write("    hex (12 17 27 18 44 49 59 50) ("+str(Ny5)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_loin)+" 1)\n") #21

# Writing edges
file.write(");\n\nedges\n(\n")

# For the profile of the shape
List_=[[List_profile_pts_edge0,List_profile_pts_edge1,List_profile_pts_edge2,List_profile_pts_edge3,List_profile_pts_edge4,List_profile_pts_edge5],[List_contour1_pts_edge0,List_contour1_pts_edge1,List_contour1_pts_edge2,List_contour1_pts_edge3,List_contour1_pts_edge4,List_contour1_pts_edge5],[List_contour2_pts_edge0,List_contour2_pts_edge1,List_contour2_pts_edge2,List_contour2_pts_edge3,List_contour2_pts_edge4,List_contour2_pts_edge5]]
List_pts=[[[1,0],[2,1],[3,2],[4,3],[5,4],[0,5]],[[7,6],[8,7],[9,8],[10,9],[11,10],[6,11]],[[13,12],[14,13],[15,14],[16,15],[17,16],[12,17]]]
for i in range(0,3):
	for k in range(0,6):
		file.write("       spline "+str(List_pts[i][k][0])+" "+str(List_pts[i][k][1])+"\n       (\n             ")
		j=0
		while 1==1:
			try:
				file.write("("+str(List_[i][k][j][0])+" "+str(List_[i][k][j][1])+" 0)\n")
				j+=1
			except IndexError:
				file.write(")\n\n")
				break
		file.write("       spline "+str(32+List_pts[i][k][0])+" "+str(32+List_pts[i][k][1])+"\n       (\n             ")
		j=0
		while 1==1:
			try:
				file.write("("+str(List_[i][k][j][0])+" "+str(List_[i][k][j][1])+" "+str(Longueur)+")\n")
				j+=1
			except IndexError:
				file.write(")\n\n")
				break

Gap=32

file.write(");\nboundary\n(\n")
file.write("Aile\n{\ntype wall;\nfaces\n(\n")
for i in range(0,L-1):
	file.write("("+str(i)+" "+str(i+1)+" "+str(i+1+Gap)+" "+str(i+Gap)+")\n")
file.write("("+str(L-1)+" 0 "+str(Gap)+" "+str(L-1+Gap)+")\n")
file.write("\n);\n}\n")

file.write("Inlet\n{\ntype patch;\nfaces\n(\n(22 29 61 54)\n(23 22 54 55)\n(24 23 55 56)\n(30 24 56 62)\n);\n}\n")
file.write("Outlet\n{\ntype patch;\nfaces\n(\n(19 28 60 51)\n(18 19 51 50)\n(27 18 50 59)\n(31 27 59 63)\n);\n}\n")
file.write("Wall\n{\ntype wall;\nfaces\n(\n(29 21 53 61)\n(21 20 52 53)\n(20 28 60 52)\n(30 25 57 62)\n(25 26 58 57)\n(26 31 63 58)\n);\n}\n")

if twoD=="true":
	file.write("FrontAndBack\n{\ntype empty;\n")
else:
	file.write("Wall2\n{\ntype wall;\n")

file.write("faces\n(\n")
for i in range(0,L-1):
	file.write("("+str(i)+" "+str(i+1)+" "+str(i+1+L)+" "+str(i+L)+")\n")
file.write("("+str(L-1)+" "+str(0)+" "+str(L)+" "+str(2*L-1)+")\n")
for i in range(0,L-1):
	file.write("("+str(L+i)+" "+str(L+i+1)+" "+str(i+1+2*L)+" "+str(i+2*L)+")\n")
file.write("("+str(2*L-1)+" "+str(L)+" "+str(2*L)+" "+str(3*L-1)+")\n")

for i in range(0,L-1):
	file.write("("+str(i+Gap)+" "+str(i+1+Gap)+" "+str(i+1+L+Gap)+" "+str(i+L+Gap)+")\n")
file.write("("+str(L-1+Gap)+" "+str(0+Gap)+" "+str(L+Gap)+" "+str(2*L-1+Gap)+")\n")
for i in range(0,L-1):
	file.write("("+str(L+i+Gap)+" "+str(L+i+1+Gap)+" "+str(i+1+2*L+Gap)+" "+str(i+2*L+Gap)+")\n")
file.write("("+str(2*L-1+Gap)+" "+str(L+Gap)+" "+str(2*L+Gap)+" "+str(3*L-1+Gap)+")\n")

file.write("(12 13 19 18)\n(13 20 28 19)\n(13 14 21 20)\n(14 22 29 21)\n(23 22 14 15)\n(24 23 15 16)\n(30 24 16 25)\n(25 26 17 16)\n(26 31 27 17)\n(17 27 18 12)\n(44 45 51 50)\n(45 52 60 51)\n(45 46 53 52)\n(46 54 61 53)\n(55 54 46 47)\n(56 55 47 48)\n(62 56 48 57)\n(57 58 49 48)\n(58 63 59 49)\n(49 59 50 44) \n);\n}\n);")
file.close()

# To help finding errors
'''
test0=np.zeros((20*int(N/2.0),2))
test1=np.zeros((20*int(N/2.0),2))
test2=np.zeros((20*int(N/2.0),2))
test3=np.zeros((20*int(N/2.0),2))
test4=np.zeros((20*int(N/2.0),2))
test5=np.zeros((20*int(N/2.0),2))
rtest0=np.zeros((20*int(N/2.0)+2*M-2,2))
rtest1=np.zeros((20*int(N/2.0)+2*M-2,2))
rtest2=np.zeros((20*int(N/2.0)+2*M-2,2))
rtest3=np.zeros((20*int(N/2.0)+2*M-2,2))
rtest4=np.zeros((20*int(N/2.0)+2*M-2,2))
rtest5=np.zeros((20*int(N/2.0)+2*M-2,2))

for i in range(0,len(List_profile_pts_edge0)):
        test0[i,0],test0[i,1]=List_profile_pts_edge0[i][0],List_profile_pts_edge0[i][1]
for i in range(0,len(List_profile_pts_edge1)):
        test1[i,0],test1[i,1]=List_profile_pts_edge1[i][0],List_profile_pts_edge1[i][1]
for i in range(0,len(List_profile_pts_edge2)):
        test2[i,0],test2[i,1]=List_profile_pts_edge2[i][0],List_profile_pts_edge2[i][1]
for i in range(0,len(List_profile_pts_edge3)):
        test3[i,0],test3[i,1]=List_profile_pts_edge3[i][0],List_profile_pts_edge3[i][1]
for i in range(0,len(List_profile_pts_edge4)):
        test4[i,0],test4[i,1]=List_profile_pts_edge4[i][0],List_profile_pts_edge4[i][1]
for i in range(0,len(List_profile_pts_edge5)):
        test5[i,0],test5[i,1]=List_profile_pts_edge5[i][0],List_profile_pts_edge5[i][1]
for i in range(0,len(List_contour2_pts_edge0)):
        rtest0[i,0],rtest0[i,1]=List_contour2_pts_edge0[i][0],List_contour2_pts_edge0[i][1]
for i in range(0,len(List_contour2_pts_edge1)):
        rtest1[i,0],rtest1[i,1]=List_contour2_pts_edge1[i][0],List_contour2_pts_edge1[i][1]
for i in range(0,len(List_contour2_pts_edge2)):
        rtest2[i,0],rtest2[i,1]=List_contour2_pts_edge2[i][0],List_contour2_pts_edge2[i][1]
for i in range(0,len(List_contour2_pts_edge3)):
        rtest3[i,0],rtest3[i,1]=List_contour2_pts_edge3[i][0],List_contour2_pts_edge3[i][1]
for i in range(0,len(List_contour2_pts_edge4)):
        rtest4[i,0],rtest4[i,1]=List_contour2_pts_edge4[i][0],List_contour2_pts_edge4[i][1]
for i in range(0,len(List_contour2_pts_edge5)):
        rtest5[i,0],rtest5[i,1]=List_contour2_pts_edge5[i][0],List_contour2_pts_edge5[i][1]

plt.plot(test0[:,0],test0[:,1],'o')
plt.plot(test1[:,0],test1[:,1],'go')
plt.plot(test2[:,0],test2[:,1],'ro')
plt.plot(test3[:,0],test3[:,1],'bo')
plt.plot(test4[:,0],test4[:,1],'yo')
plt.plot(test5[:,0],test5[:,1],'o','indigo')
plt.plot(rtest0[:,0],rtest0[:,1],'o')
plt.plot(rtest1[:,0],rtest1[:,1],'go')
plt.plot(rtest2[:,0],rtest2[:,1],'ro')
plt.plot(rtest3[:,0],rtest3[:,1],'bo')
plt.plot(rtest4[:,0],rtest4[:,1],'yo')
plt.plot(rtest5[:,0],rtest5[:,1],'o','indigo')
plt.show()
plt.plot(Profile[:,0],Profile[:,1],'g-')
plt.plot(List_profile[:,0],List_profile[:,1],'go')
plt.plot(Contour1[:,0],Contour1[:,1],'r-')
plt.plot(List_etendu1[:,0],List_etendu1[:,1],'ro')
plt.plot(Contour2[:,0],Contour2[:,1],'y-')
plt.plot(List_etendu2[:,0],List_etendu2[:,1],'yo')
plt.plot(List_bloc[:,0],List_bloc[:,1],'bo')
plt.show()
'''

