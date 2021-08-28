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
ratio_loinH=2 # Ratio Top | Ratio of the last cell over the first one
ratio_loinG=2 # Ratio Left
ratio_loinB=2 # Ratio Bottom
R_loinD=1.2 # Ratio Right between the size of a cell and the size of the next one

if Material=='Water':
	rho=1000
	mu=1e-3
elif Material=='Air':
	rho=1
	mu=1e-6
else:
	print('Error in the definition of the material')
	sys.exit()

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
delta=corde/2.0 # Thickness of the successive layers that break the propagation of the mesh stemming from the boundary layer mesh

#~~~~~~~~~~~~~~~~~~~~~ Features of the mesh of the boundary layer ~~~~~~~~~~~~~~~~~~~~~#
Re=rho*U_inf*corde/mu
Cf=0.026*Re**(-1.0/7.0)
t_wall=0.5*Cf*rho*U_inf**2
u_etoile=sqrt(t_wall/rho)
y=2*mu*y_plus/(rho*u_etoile) # Size of the first cell - The '2*' is for the centre to be at y metre from the wall (OpenFOAM uses collocated grids)

Max=int((log(Taille_maille)-log(y))/log(R)) # Number of cells between the wing and the cell with the same size as the parameter Taille_maille
Length=0 # Layer thickness
for n in range(0,Max+1):
	Length+=y*R**n 

ratio_aile=R**Max
l_cible=0.00001 # Size of the cell between both areas
d=1.5 # Division by d of the number of cells in each layer on the right of the wing
typ="merge"
if y>l_cible:
	n_casse=0
elif typ=="merge":
	n_casse=int((log(Max)-log(Length/Taille_maille))/log(d))
	str_0=[str(1.0/ratio_aile)]
	str_1=[str(ratio_aile)]
	for i in range(0,n_casse+1):
		str_0.append(str(1.0/ratio_aile))
		str_1.append(str(ratio_aile))
else:
	leng=0
	n=0
	str_0=[str(1.0/ratio_aile)]
	str_1=[str(ratio_aile)]
	Fraction_cell=[]
	while leng<l_cible:
		leng+=y*R**n
		n+=1
	Fraction_len=leng/Length	
	alpha=R**(Max-n)
	beta=R**(n-1.0)
	Fraction_cell.append(0.8*n/(Max+1.0))
	k=0
	while int(n/d**k)>4:
		k+=1
		Fraction_cell.append(0.8*int(n/d**k)/(Max+1.0-n+int(n/d**k)))
	n_casse=k
	for i in range(0,n_casse+1):
		str_i="(("+str(1-Fraction_len)+" "+str(1-Fraction_cell[i])+" "+str(1.0/alpha)+") ("+str(Fraction_len)+" "+str(Fraction_cell[i])+" "+str(1.0/beta)+"))"
		str_0.append(str_i)
		str_i="(("+str(Fraction_len)+" "+str(Fraction_cell[i])+" "+str(beta)+") ("+str(1-Fraction_len)+" "+str(1-Fraction_cell[i])+" "+str(alpha)+"))" 
		str_1.append(str_i)

# Storage of the points of the profile in the vector Profile
Profile=np.zeros((2*int(N/2.0),2))
for i in range(0,int(N/2.0)):
	Profile[i,0],Profile[i,1]=Table[2*i],Table[2*i+1]

# Barycentre
X_g=0.5*(max([Profile[i,0] for i in range(0,int(N/2.0))])-min([Profile[i,0] for i in range(0,int(N/2.0))]))
Y_g=0

# Writing the extended contours
Contour1=np.zeros((N,2))
Contour2=np.zeros((N,2))
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

# Symmetry with respect to the axis y=0 (the profile is assumed to be symmetrical)
for i in range(0,int(N/2.0)):
	Contour1[i+int(N/2.0),0],Contour1[i+int(N/2.0),1]=Contour1[i,0],-Contour1[i,1]
	Contour2[i+int(N/2.0),0],Contour2[i+int(N/2.0),1]=Contour2[i,0],-Contour2[i,1]
	if i<int(N/2.0):
		Profile[i+int(N/2.0),0],Profile[i+int(N/2.0),1]=Profile[i,0],-Profile[i,1]

L=4															
List_profile=np.zeros((L,2))
List_etendu1=np.zeros((L+1,2))
List_etendu2=np.zeros((L+1,2))
List_bloc=np.zeros((16+2*(n_casse+1),2))												
if n_casse>0:
	List_casse=np.zeros(((n_casse+1)*5,2))

indice_profile,indice_etendu1,indice_etendu2=[],[],[]

# We look for the point that minimizes the distance with the line X=Xg+-D1_2*chord/2 
D1=0.75
D2=0.66

List_profile[0,0],List_profile[0,1]=X_max,Y_max
indice_profile.append(int(N/2))	
norm_profile=1e5
for i in range(0,int(N/2.0)):
	norm_profile_test=(Profile[i,0]-(X_g-D1*corde/2))**2
	if norm_profile_test<norm_profile:
		norm_profile=norm_profile_test
		indice=i
indice_profile.append(indice)
List_profile[1,0],List_profile[1,1]=Profile[indice,0],Profile[indice,1]
		
List_etendu1[0,0],List_etendu1[0,1]=Contour1[int(N/2)-1,0],Contour1[int(N/2)-1,1]
indice_etendu1.append(int(N/2))
norm_etendu1=1e5
for i in range(0,int(N/2.0)):
	norm_etendu1_test=(Contour1[i,0]-(X_g-D1*corde/2))**2
	if norm_etendu1_test<norm_etendu1:
		norm_etendu1=norm_etendu1_test
		indice=i
indice_etendu1.append(indice)
List_etendu1[1,0],List_etendu1[1,1]=Contour1[indice,0],Contour1[indice,1]										

for k in range(0,2):
	norm_etendu2=1e5
	for i in range(0,int(N/2.0)):
		norm_etendu2_test=(Contour2[i,0]-(X_g+(-1)**k*D2*(corde+2*H)/2))**2
		if norm_etendu2_test<norm_etendu2:
			norm_etendu2=norm_etendu2_test
			indice=i
	indice_etendu2.append(indice)
	List_etendu2[k,0],List_etendu2[k,1]=Contour2[indice,0],Contour2[indice,1]
											

List_profile[2,0],List_profile[2,1]=0,0
List_profile[3,0],List_profile[3,1]=List_profile[1,0],-List_profile[1,1]					
List_etendu1[2,0],List_etendu1[2,1]=-Length,0
List_etendu1[3,0],List_etendu1[3,1]=List_etendu1[1,0],-List_etendu1[1,1]
List_etendu1[4,0],List_etendu1[4,1]=List_etendu1[0,0],-List_etendu1[0,1]							
List_etendu2[2,0],List_etendu2[2,1]=-H,0
List_etendu2[3,0],List_etendu2[3,1]=List_etendu2[1,0],-List_etendu2[1,1]
List_etendu2[4,0],List_etendu2[4,1]=List_etendu2[0,0],-List_etendu2[0,1]									

# Rotation
for i in range(0,N):
	X_profile,Y_profile=Profile[i,0]-X_c,Profile[i,1]-Y_c
	Profile[i,0],Profile[i,1]=X_c+cos(Rotation)*X_profile+sin(Rotation)*Y_profile,Y_c-sin(Rotation)*X_profile+cos(Rotation)*Y_profile
	X_contour1,Y_contour1=Contour1[i,0]-X_c,Contour1[i,1]-Y_c
	Contour1[i,0],Contour1[i,1]=X_c+cos(Rotation)*X_contour1+sin(Rotation)*Y_contour1,Y_c-sin(Rotation)*X_contour1+cos(Rotation)*Y_contour1
	X_contour2,Y_contour2=Contour2[i,0]-X_c,Contour2[i,1]-Y_c
	Contour2[i,0],Contour2[i,1]=X_c+cos(Rotation)*X_contour2+sin(Rotation)*Y_contour2,Y_c-sin(Rotation)*X_contour2+cos(Rotation)*Y_contour2

for j in range(0,L+1):
	if j<L:
		X_profile,Y_profile=List_profile[j,0]-X_c,List_profile[j,1]-Y_c
		List_profile[j,0],List_profile[j,1]=X_c+cos(Rotation)*X_profile+sin(Rotation)*Y_profile,Y_c-sin(Rotation)*X_profile+cos(Rotation)*Y_profile
	X_contour1,Y_contour1=List_etendu1[j,0]-X_c,List_etendu1[j,1]-Y_c
	List_etendu1[j,0],List_etendu1[j,1]=X_c+cos(Rotation)*X_contour1+sin(Rotation)*Y_contour1,Y_c-sin(Rotation)*X_contour1+cos(Rotation)*Y_contour1
	X_contour2,Y_contour2=List_etendu2[j,0]-X_c,List_etendu2[j,1]-Y_c
	List_etendu2[j,0],List_etendu2[j,1]=X_c+cos(Rotation)*X_contour2+sin(Rotation)*Y_contour2,Y_c-sin(Rotation)*X_contour2+cos(Rotation)*Y_contour2

# Adding points on the block
List_bloc[0,0],List_bloc[0,1]=alpha_longueur+X_g,List_etendu2[0,1] #14
List_bloc[1,0],List_bloc[1,1]=alpha_longueur+X_g,alpha_largeur+Y_g #15
List_bloc[2,0],List_bloc[2,1]=List_etendu2[0,0],alpha_largeur+Y_g #16
List_bloc[3,0],List_bloc[3,1]=List_etendu2[1,0],alpha_largeur+Y_g #17
List_bloc[4,0],List_bloc[4,1]=-alpha_longueur+X_g,alpha_largeur+Y_g #18
List_bloc[5,0],List_bloc[5,1]=-alpha_longueur+X_g,List_etendu2[1,1] #19
List_bloc[6,0],List_bloc[6,1]=-alpha_longueur+X_g,List_etendu2[2,1] #20
List_bloc[7,0],List_bloc[7,1]=-alpha_longueur+X_g,List_etendu2[3,1] #21
List_bloc[8,0],List_bloc[8,1]=-alpha_longueur+X_g,-alpha_largeur+Y_g #22
List_bloc[9,0],List_bloc[9,1]=List_etendu2[3,0],-alpha_largeur+Y_g #23
List_bloc[10,0],List_bloc[10,1]=List_etendu2[4,0],-alpha_largeur+Y_g #24
List_bloc[11,0],List_bloc[11,1]=alpha_longueur+X_g,-alpha_largeur+Y_g #25
List_bloc[12,0],List_bloc[12,1]=alpha_longueur+X_g,List_etendu2[4,1] #26
List_bloc[13,0],List_bloc[13,1]=alpha_longueur+X_g,List_etendu1[4,1] #27
List_bloc[14,0],List_bloc[14,1]=alpha_longueur+X_g,List_profile[0,1] #28
List_bloc[15,0],List_bloc[15,1]=alpha_longueur+X_g,List_etendu1[0,1] #29

if n_casse>0:
	for n in range(0,n_casse+1):
		List_bloc[16+2*n,0],List_bloc[16+2*n,1]=List_bloc[10,0]+delta*(n+1),-alpha_largeur+Y_g
		List_bloc[16+2*n+1,0],List_bloc[16+2*n+1,1]=List_bloc[2,0]+delta*(n+1),alpha_largeur+Y_g
		List_casse[5*n,0],List_casse[5*n,1]=List_etendu2[4,0]+delta*(n+1),List_etendu2[4,1]
		List_casse[5*n+1,0],List_casse[5*n+1,1]=List_etendu1[4,0]+delta*(n+1),List_etendu1[4,1]
		List_casse[5*n+2,0],List_casse[5*n+2,1]=List_profile[0,0]+delta*(n+1),List_profile[0,1]
		List_casse[5*n+3,0],List_casse[5*n+3,1]=List_etendu1[0,0]+delta*(n+1),List_etendu1[0,1]
		List_casse[5*n+4,0],List_casse[5*n+4,1]=List_etendu2[0,0]+delta*(n+1),List_etendu2[0,1]

# List of points that will be in the splines
List_profile_pts_edge0=[]
List_profile_pts_edge1=[]
List_profile_pts_edge2=[]
List_profile_pts_edge3=[]
List_contour1_pts_edge0=[]
List_contour1_pts_edge1=[]
List_contour1_pts_edge2=[]
List_contour1_pts_edge3=[]
List_contour2_pts_edge0=[]
List_contour2_pts_edge1=[]
List_contour2_pts_edge2=[]
List_contour2_pts_edge3=[]

for i in range(0,int(N/2.0)):
	if i>indice_profile[1]:
		List_profile_pts_edge0.append([Profile[i,0],Profile[i,1]])
	else:
		List_profile_pts_edge1.append([Profile[i,0],Profile[i,1]])
for i in range(N-1,int(N/2.0),-1):
	if i>indice_profile[1]+int(N/2.0):
		List_profile_pts_edge3.append([Profile[i,0],Profile[i,1]])
	else: 
		List_profile_pts_edge2.append([Profile[i,0],Profile[i,1]])

for i in range(0,int(N/2)):
	if i>indice_etendu1[1]:
		List_contour1_pts_edge0.append([Contour1[i,0],Contour1[i,1]])
	else:
		List_contour1_pts_edge1.append([Contour1[i,0],Contour1[i,1]])
for i in range(N-1,int(N/2),-1):
	if i>indice_etendu1[1]+int(N/2):
		List_contour1_pts_edge3.append([Contour1[i,0],Contour1[i,1]])
	else:
		List_contour1_pts_edge2.append([Contour1[i,0],Contour1[i,1]])

for i in range(0,int(N/2)):
	if i>indice_etendu2[1]:
		List_contour2_pts_edge0.append([Contour2[i,0],Contour2[i,1]])
	else:
		List_contour2_pts_edge1.append([Contour2[i,0],Contour2[i,1]])
for i in range(N-1,int(N/2),-1):
	if i>indice_etendu2[1]+int(N/2):
		List_contour2_pts_edge3.append([Contour2[i,0],Contour2[i,1]])
	else:
		List_contour2_pts_edge2.append([Contour2[i,0],Contour2[i,1]])

# Calculation of the length of the edges
len0=len1=len2=len3=0

for i in range(1,len(List_contour1_pts_edge0)):
	len0+=sqrt((List_contour1_pts_edge0[i][0]-List_contour1_pts_edge0[i-1][0])**2+(List_contour1_pts_edge0[i][1]-List_contour1_pts_edge0[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge1)):
	len1+=sqrt((List_contour1_pts_edge1[i][0]-List_contour1_pts_edge1[i-1][0])**2+(List_contour1_pts_edge1[i][1]-List_contour1_pts_edge1[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge2)):
	len2+=sqrt((List_contour1_pts_edge2[i][0]-List_contour1_pts_edge2[i-1][0])**2+(List_contour1_pts_edge2[i][1]-List_contour1_pts_edge2[i-1][1])**2)
for i in range(1,len(List_contour1_pts_edge3)):
	len3+=sqrt((List_contour1_pts_edge3[i][0]-List_contour1_pts_edge3[i-1][0])**2+(List_contour1_pts_edge3[i][1]-List_contour1_pts_edge3[i-1][1])**2)

# Writing the file blockMeshDict
file=open('system/blockMeshDict','w')
file.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      blockMeshDict;\n}\n\n")
file.write("convertToMeters 1;\n\nvertices\n(\n")

# Les points
for i in range(0,4):
	file.write("    ("+str(List_profile[i,0])+" "+str(List_profile[i,1])+" 0)   //"+str(i)+"\n")

for i in range(0,5):
	file.write("    ("+str(List_etendu1[i,0])+" "+str(List_etendu1[i,1])+" 0)   //"+str(i+4)+"\n")

for i in range(0,5):
	file.write("    ("+str(List_etendu2[i,0])+" "+str(List_etendu2[i,1])+" 0)   //"+str(i+9)+"\n")

for i in range(0,16):
	file.write("    ("+str(List_bloc[i,0])+" "+str(List_bloc[i,1])+" 0)   //"+str(i+14)+"\n")

if n_casse>0:
	for i in range(0,n_casse+1):
		file.write("    ("+str(List_bloc[16+2*i,0])+" "+str(List_bloc[16+2*i,1])+" 0)   //"+str(7*i+30)+"\n")
		for n in range(0,5):
			file.write("    ("+str(List_casse[n+5*i,0])+" "+str(List_casse[n+5*i,1])+" 0)   //"+str(7*i+31+n)+"\n")
		file.write("    ("+str(List_bloc[17+2*i,0])+" "+str(List_bloc[17+2*i,1])+" 0)   //"+str(7*i+36)+"\n")
	for i in range(0,n_casse+1):
		file.write("    ("+str(List_bloc[16+2*i,0])+" "+str(List_bloc[16+2*i,1])+" 0)   //"+str(7*i+7*n_casse+37)+"\n")
		for n in range(0,5):
			file.write("    ("+str(List_casse[n+5*i,0])+" "+str(List_casse[n+5*i,1])+" 0)   //"+str(7*i+7*n_casse+38+n)+"\n")
		file.write("    ("+str(List_bloc[17+2*i,0])+" "+str(List_bloc[17+2*i,1])+" 0)   //"+str(7*i+7*n_casse+43)+"\n")
	Gap=7*2*n_casse+44
else:
	Gap=30

for i in range(0,4):
	file.write("    ("+str(List_profile[i,0])+" "+str(List_profile[i,1])+" "+str(Longueur)+")   //"+str(i+Gap)+"\n")

for i in range(0,5):
	file.write("    ("+str(List_etendu1[i,0])+" "+str(List_etendu1[i,1])+" "+str(Longueur)+")   //"+str(i+Gap+4)+"\n")

for i in range(0,5):
	file.write("    ("+str(List_etendu2[i,0])+" "+str(List_etendu2[i,1])+" "+str(Longueur)+")   //"+str(i+Gap+9)+"\n")

for i in range(0,16):
	file.write("    ("+str(List_bloc[i,0])+" "+str(List_bloc[i,1])+" "+str(Longueur)+")   //"+str(i+Gap+14)+"\n")

if n_casse>0:
	for i in range(0,n_casse+1):
		file.write("    ("+str(List_bloc[16+2*i,0])+" "+str(List_bloc[16+2*i,1])+" "+str(Longueur)+")   //"+str(7*i+30+Gap)+"\n")
		for n in range(0,5):
			file.write("    ("+str(List_casse[n+5*i,0])+" "+str(List_casse[n+5*i,1])+" "+str(Longueur)+")   //"+str(7*i+31+n+Gap)+"\n")
		file.write("    ("+str(List_bloc[17+2*i,0])+" "+str(List_bloc[17+2*i,1])+" "+str(Longueur)+")   //"+str(7*i+36+Gap)+"\n")
	for i in range(0,n_casse+1):
		file.write("    ("+str(List_bloc[16+2*i,0])+" "+str(List_bloc[16+2*i,1])+" "+str(Longueur)+")   //"+str(7*i+7*n_casse+37+Gap)+"\n")
		for n in range(0,5):
			file.write("    ("+str(List_casse[n+5*i,0])+" "+str(List_casse[n+5*i,1])+" "+str(Longueur)+")   //"+str(7*i+7*n_casse+38+n+Gap)+"\n")
		file.write("    ("+str(List_bloc[17+2*i,0])+" "+str(List_bloc[17+2*i,1])+" "+str(Longueur)+")   //"+str(7*i+7*n_casse+43+Gap)+"\n")

Ny0=int(len0/Taille_maille) # Number of cells in the normal to the wing direction in the zone 0
Ny1=int(len1/Taille_maille) # Number of cells in the normal to the wing direction in the zone 1
Ny2=int(len2/Taille_maille) # Number of cells in the normal to the wing direction in the zone 2
Ny3=int(len3/Taille_maille) # Number of cells in the normal to the wing direction in the zone 3

N_int=int(H/Taille_maille) # Number of cells in the normal to the wing direction in the intermediate layer

N_loinD=int(log(1-(1-R_loinD)*(alpha_largeur-corde)/Taille_maille)/log(R_loinD)-1)+1
ratio_loinD=R_loinD**N_loinD

List1=[[1,0,4,5],[2,1,5,6],[3,2,6,7],[0,3,7,8],[5,4,9,10],[6,5,10,11],[7,6,11,12],[8,7,12,13],[10,9,16,17],[19,10,17,18],[11,10,19,20],[12,11,20,21],[12,21,22,23],[13,12,23,24]]
for i in range(0,len(List1)):
	List1[i][0]+=Gap
	List1[i][1]+=Gap
	List1[i][2]+=Gap
	List1[i][3]+=Gap

file.write(");\n\nblocks\n(\n")
file.write("    hex (1 0 4 5 "+str(List1[0][0])+" "+str(List1[0][1])+" "+str(List1[0][2])+" "+str(List1[0][3])+") ("+str(Ny0)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #0 
file.write("    hex (2 1 5 6 "+str(List1[1][0])+" "+str(List1[1][1])+" "+str(List1[1][2])+" "+str(List1[1][3])+") ("+str(Ny1)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #1
file.write("    hex (3 2 6 7 "+str(List1[2][0])+" "+str(List1[2][1])+" "+str(List1[2][2])+" "+str(List1[2][3])+") ("+str(Ny2)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #2
file.write("    hex (0 3 7 8 "+str(List1[3][0])+" "+str(List1[3][1])+" "+str(List1[3][2])+" "+str(List1[3][3])+") ("+str(Ny3)+" "+str(Max+1)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_aile)+" 1)\n") #3
file.write("    hex (5 4 9 10 "+str(List1[4][0])+" "+str(List1[4][1])+" "+str(List1[4][2])+" "+str(List1[4][3])+") ("+str(Ny0)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #4
file.write("    hex (6 5 10 11 "+str(List1[5][0])+" "+str(List1[5][1])+" "+str(List1[5][2])+" "+str(List1[5][3])+") ("+str(Ny1)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #5
file.write("    hex (7 6 11 12 "+str(List1[6][0])+" "+str(List1[6][1])+" "+str(List1[6][2])+" "+str(List1[6][3])+") ("+str(Ny2)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #6
file.write("    hex (8 7 12 13 "+str(List1[7][0])+" "+str(List1[7][1])+" "+str(List1[7][2])+" "+str(List1[7][3])+") ("+str(Ny3)+" "+str(N_int)+" "+str(Nz)+") simpleGrading (1 1 1)\n") #7
file.write("    hex (10 9 16 17 "+str(List1[8][0])+" "+str(List1[8][1])+" "+str(List1[8][2])+" "+str(List1[8][3])+") ("+str(Ny0)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_loinH)+" 1)\n") #8
file.write("    hex (19 10 17 18 "+str(List1[9][0])+" "+str(List1[9][1])+" "+str(List1[9][2])+" "+str(List1[9][3])+") ("+str(N_ext)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(1.0/ratio_loinG)+" "+str(ratio_loinH)+" 1)\n") #9
file.write("    hex (11 10 19 20 "+str(List1[10][0])+" "+str(List1[10][1])+" "+str(List1[10][2])+" "+str(List1[10][3])+") ("+str(Ny1)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_loinG)+" 1)\n") #10
file.write("    hex (12 11 20 21 "+str(List1[11][0])+" "+str(List1[11][1])+" "+str(List1[11][2])+" "+str(List1[11][3])+") ("+str(Ny2)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(ratio_loinG)+" 1)\n") #11
file.write("    hex (12 21 22 23 "+str(List1[12][0])+" "+str(List1[12][1])+" "+str(List1[12][2])+" "+str(List1[12][3])+") ("+str(N_ext)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(ratio_loinG)+" "+str(1.0/ratio_loinB)+" 1)\n") #12
file.write("    hex (13 12 23 24 "+str(List1[13][0])+" "+str(List1[13][1])+" "+str(List1[13][2])+" "+str(List1[13][3])+") ("+str(Ny3)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading (1 "+str(1.0/ratio_loinB)+" 1)\n") #13

List2=[[24,13,8,0,4,9,16]]
if n_casse>0:
	for i in range(0,2*(n_casse+1)):
		List2.append([30+7*i+n for n in range(0,7)])
	List3=[]
	for i in range(0,n_casse+1):
		List_int=[]
		for j in range(0,6):
			List_int.append([List2[i][j],List2[i][j+1],List2[i+n_casse+2][j+1],List2[i+n_casse+2][j]])
		List3.append(List_int)
	l=n_casse+1
	Ratio_casse=[] 
	N_casse=[]
	R_loinD=1.1
	size=Taille_maille
	for i in range(0,n_casse+1):
		N_loinD=int(log(1-(1-R_loinD)*(delta)/size)/log(R_loinD)-1)+1
		Ratio_casse.append(R_loinD**N_loinD)
		N_casse.append(N_loinD)
		size=size*R_loinD**N_loinD
	N_loinD=int(log(1-(1-R_loinD)*(alpha_largeur-corde-n_casse*delta)/size)/log(R_loinD)-1)+1
	ratio_loinD=R_loinD**N_loinD
	for n in range(0,n_casse+1):
		file.write("    hex ("+str(List3[n][0][3])+" "+str(List3[n][0][2])+" "+str(List3[n][0][1])+" "+str(List3[n][0][0])+" "+str(List3[n][0][3]+Gap)+" "+str(List3[n][0][2]+Gap)+" "+str(List3[n][0][1]+Gap)+" "+str(List3[n][0][0]+Gap)+") ("+str(N_ext)+" "+str(N_casse[n])+" "+str(Nz)+") simpleGrading ("+str(ratio_loinB)+" "+str(1.0/Ratio_casse[n])+" 1)\n")
		file.write("    hex ("+str(List3[n][1][3])+" "+str(List3[n][1][2])+" "+str(List3[n][1][1])+" "+str(List3[n][1][0])+" "+str(List3[n][1][3]+Gap)+" "+str(List3[n][1][2]+Gap)+" "+str(List3[n][1][1]+Gap)+" "+str(List3[n][1][0]+Gap)+") ("+str(N_int)+" "+str(N_casse[n])+" "+str(Nz)+") simpleGrading (1 "+str(1.0/Ratio_casse[n])+" 1)\n")
		file.write("    hex ("+str(List3[n][2][3])+" "+str(List3[n][2][2])+" "+str(List3[n][2][1])+" "+str(List3[n][2][0])+" "+str(List3[n][2][3]+Gap)+" "+str(List3[n][2][2]+Gap)+" "+str(List3[n][2][1]+Gap)+" "+str(List3[n][2][0]+Gap)+") ("+str(int((Max+1)/d**n))+" "+str(N_casse[n])+" "+str(Nz)+") simpleGrading ("+str(str_0[n])+" "+str(1.0/Ratio_casse[n])+" 1)\n")
		file.write("    hex ("+str(List3[n][3][3])+" "+str(List3[n][3][2])+" "+str(List3[n][3][1])+" "+str(List3[n][3][0])+" "+str(List3[n][3][3]+Gap)+" "+str(List3[n][3][2]+Gap)+" "+str(List3[n][3][1]+Gap)+" "+str(List3[n][3][0]+Gap)+") ("+str(int((Max+1)/d**n))+" "+str(N_casse[n])+" "+str(Nz)+") simpleGrading ("+str(str_1[n])+" "+str(1.0/Ratio_casse[n])+" 1)\n")
		file.write("    hex ("+str(List3[n][4][3])+" "+str(List3[n][4][2])+" "+str(List3[n][4][1])+" "+str(List3[n][4][0])+" "+str(List3[n][4][3]+Gap)+" "+str(List3[n][4][2]+Gap)+" "+str(List3[n][4][1]+Gap)+" "+str(List3[n][4][0]+Gap)+") ("+str(N_int)+" "+str(N_casse[n])+" "+str(Nz)+") simpleGrading (1 "+str(1.0/Ratio_casse[n])+" 1)\n")
		file.write("    hex ("+str(List3[n][5][3])+" "+str(List3[n][5][2])+" "+str(List3[n][5][1])+" "+str(List3[n][5][0])+" "+str(List3[n][5][3]+Gap)+" "+str(List3[n][5][2]+Gap)+" "+str(List3[n][5][1]+Gap)+" "+str(List3[n][5][0]+Gap)+") ("+str(N_ext)+" "+str(N_casse[n])+" "+str(Nz)+") simpleGrading ("+str(ratio_loinH)+" "+str(1.0/Ratio_casse[n])+" 1)\n")
else:
	l=0

file.write("    hex (26 "+str(List2[l][1])+" "+str(List2[l][0])+" 25 "+str(26+Gap)+" "+str(List2[l][1]+Gap)+" "+str(List2[l][0]+Gap)+" "+str(25+Gap)+") ("+str(N_loinD+1)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(1.0/ratio_loinD)+" "+str(1.0/ratio_loinB)+" 1)\n") #14
file.write("    hex ("+str(List2[l][1])+" 26 27 "+str(List2[l][2])+" "+str(List2[l][1]+Gap)+" "+str(26+Gap)+" "+str(27+Gap)+" "+str(List2[l][2]+Gap)+") ("+str(N_loinD+1)+" "+str(N_int)+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" 1 1)\n") #15
if n_casse>0:
	file.write("    hex ("+str(List2[l][2])+" 27 28 "+str(List2[l][3])+" "+str(List2[l][2]+Gap)+" "+str(27+Gap)+" "+str(28+Gap)+" "+str(List2[l][3]+Gap)+") ("+str(N_loinD+1)+" "+str(int(Length/Taille_maille))+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" "+str(1.0)+" 1)\n") #16 str(1.0/ratio_aile)
	file.write("    hex ("+str(List2[l][3])+" 28 29 "+str(List2[l][4])+" "+str(List2[l][3]+Gap)+" "+str(28+Gap)+" "+str(29+Gap)+" "+str(List2[l][4]+Gap)+") ("+str(N_loinD+1)+" "+str(int(Length/Taille_maille))+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" "+str(1.0)+" 1)\n") #17 str(ratio_aile)
else:
	file.write("    hex ("+str(List2[l][2])+" 27 28 "+str(List2[l][3])+" "+str(List2[l][2]+Gap)+" "+str(27+Gap)+" "+str(28+Gap)+" "+str(List2[l][3]+Gap)+") ("+str(N_loinD+1)+" "+str(int((Max+1)/d**n_casse))+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" "+str(1.0/ratio_aile)+" 1)\n") #16 
	file.write("    hex ("+str(List2[l][3])+" 28 29 "+str(List2[l][4])+" "+str(List2[l][3]+Gap)+" "+str(28+Gap)+" "+str(29+Gap)+" "+str(List2[l][4]+Gap)+") ("+str(N_loinD+1)+" "+str(int((Max+1)/d**n_casse))+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" "+str(ratio_aile)+" 1)\n") #17 
file.write("    hex ("+str(List2[l][4])+" 29 14 "+str(List2[l][5])+" "+str(List2[l][4]+Gap)+" "+str(29+Gap)+" "+str(14+Gap)+" "+str(List2[l][5]+Gap)+") ("+str(N_loinD+1)+" "+str(N_int)+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" 1 1)\n") #18
file.write("    hex ("+str(List2[l][5])+" 14 15 "+str(List2[l][6])+" "+str(List2[l][5]+Gap)+" "+str(14+Gap)+" "+str(Gap+15)+" "+str(List2[l][6]+Gap)+") ("+str(N_loinD+1)+" "+str(N_ext)+" "+str(Nz)+") simpleGrading ("+str(ratio_loinD)+" "+str(ratio_loinH)+" 1)\n") #19

# Writing the edges
file.write(");\n\nedges\n(\n")

# For the profile of the wing
List_=[[List_profile_pts_edge0,List_profile_pts_edge1,List_profile_pts_edge2,List_profile_pts_edge3],[List_contour1_pts_edge0,List_contour1_pts_edge1,List_contour1_pts_edge2,List_contour1_pts_edge3],[List_contour2_pts_edge0,List_contour2_pts_edge1,List_contour2_pts_edge2,List_contour2_pts_edge3]]
List_pts=[[[1,0],[2,1],[3,2],[0,3]],[[5,4],[6,5],[7,6],[8,7]],[[10,9],[11,10],[12,11],[13,12]]]
for i in range(0,3):
	for k in range(0,4):
		file.write("       spline "+str(List_pts[i][k][0])+" "+str(List_pts[i][k][1])+"\n       (\n             ")
		j=0
		while 1==1:
			try:
				file.write("("+str(List_[i][k][j][0])+" "+str(List_[i][k][j][1])+" 0)\n")
				j+=1
			except IndexError:
				file.write(")\n\n")
				break
		file.write("       spline "+str(Gap+List_pts[i][k][0])+" "+str(Gap+List_pts[i][k][1])+"\n       (\n             ")
		j=0
		while 1==1:
			try:
				file.write("("+str(List_[i][k][j][0])+" "+str(List_[i][k][j][1])+" "+str(Longueur)+")\n")
				j+=1
			except IndexError:
				file.write(")\n\n")
				break

file.write(");\nboundary\n(\n")
file.write("Aile\n{\ntype wall;\nfaces\n(\n")
for i in range(0,L-1):
	file.write("("+str(i)+" "+str(i+1)+" "+str(i+1+Gap)+" "+str(i+Gap)+")\n")
file.write("("+str(L-1)+" 0 "+str(Gap)+" "+str(L-1+Gap)+")\n")
file.write("\n);\n}\n")

file.write("Inlet\n{\ntype patch;\nfaces\n(\n(19 18 "+str(18+Gap)+" "+str(19+Gap)+")\n(20 19 "+str(19+Gap)+" "+str(20+Gap)+")\n(21 20 "+str(20+Gap)+" "+str(21+Gap)+")\n(22 21 "+str(21+Gap)+" "+str(22+Gap)+")\n);\n}\n")
file.write("Outlet\n{\ntype patch;\nfaces\n(\n(25 26 "+str(26+Gap)+" "+str(25+Gap)+")\n(26 27 "+str(27+Gap)+" "+str(26+Gap)+")\n(27 28 "+str(28+Gap)+" "+str(27+Gap)+")\n(28 29 "+str(29+Gap)+" "+str(28+Gap)+")\n(29 14 "+str(14+Gap)+" "+str(29+Gap)+")\n(14 15 "+str(Gap+15)+" "+str(Gap+14)+")\n);\n}\n")

if n_casse>0:
	for n in range(1,n_casse+2):
		file.write("master"+str(n-1)+"\n{\ntype patch;\nfaces\n(\n")
		for j in range(0,6):
			file.write("("+str(List2[n][j])+" "+str(List2[n][j+1])+" "+str(List2[n][j+1]+Gap)+" "+str(List2[n][j]+Gap)+")\n")
		file.write(");\n}\n")
	for n in range(n_casse+2,2*n_casse+3):
		file.write("slave"+str(n-n_casse-2)+"\n{\ntype patch;\nfaces\n(\n")
		for j in range(0,6):		
			file.write("("+str(List2[n][j])+" "+str(List2[n][j+1])+" "+str(List2[n][j+1]+Gap)+" "+str(List2[n][j]+Gap)+")\n")
		file.write(");\n}\n")

file.write("Wall\n{\ntype wall;\nfaces\n(\n(22 23 "+str(23+Gap)+" "+str(22+Gap)+")\n(23 24 "+str(24+Gap)+" "+str(23+Gap)+")\n("+str(List2[l][0])+" 25 "+str(25+Gap)+" "+str(List2[l][0]+Gap)+")\n("+str(List2[l][6])+" 15 "+str(15+Gap)+" "+str(List2[l][6]+Gap)+")\n(17 16 "+str(16+Gap)+" "+str(17+Gap)+")\n(18 17 "+str(17+Gap)+" "+str(18+Gap)+")\n")
if n_casse>0:
	for n in range(0,n_casse+1):
		file.write("("+str(List2[n][0])+" "+str(List2[n+n_casse+2][0])+" "+str(List2[n+n_casse+2][0]+Gap)+" "+str(List2[n][0]+Gap)+")\n")
		file.write("("+str(List2[n][6])+" "+str(List2[n+n_casse+2][6])+" "+str(List2[n+n_casse+2][6]+Gap)+" "+str(List2[n][6]+Gap)+")\n")
	file.write(");\n}\n")
else:
	file.write(");\n}\n")

file.write("\n);\nmergePatchPairs\n(\n")
if n_casse>0:
	for n in range(0,n_casse+1):
		file.write("(master"+str(n)+" slave"+str(n)+")\n")
file.write(")\n;")

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

