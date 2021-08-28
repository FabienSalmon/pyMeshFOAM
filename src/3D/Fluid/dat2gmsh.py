import numpy as np
from numpy import *
import sys
import os

# Input parameters
Chemin_GMSH=sys.argv[1] 
twoD=sys.argv[2]
Longueur=float(sys.argv[3])
Ratio=float(sys.argv[4])
File=sys.argv[5]
X_r=float(sys.argv[6])
Rotation=float(sys.argv[7])*3.141592/180.0
xx=int(sys.argv[8])

nb_discr=100 # Increase this parameter to improve the mesh number in the direction of the span (nb_discr/metre) | useless in 2D
Saumon='Circle' # For the wing tip to be circular | To have a flat wing tip, just write ''

# Writing process of the input file in case of a NACA00xx profile
if xx!=0:
    file_=open(File,'w')
    N0=10
    N1=19
    N2=800
    for k in range(0,N0):
        file_.write(str(k*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*0.0001)-0.126*k*0.0001-0.3516*(k*0.0001)**2+0.2843*(k*0.0001)**3-0.1015*(k*0.0001)**4))+"\n")
    for k in range(0,N1):
        file_.write(str(k*0.001+N0*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*0.001+N0*0.0001)-0.126*(k*0.001+N0*0.0001)-0.3516*(k*0.001+N0*0.0001)**2+0.2843*(k*0.001+N0*0.0001)**3-0.1015*(k*0.001+N0*0.0001)**4))+"\n")
    dk=0.8/N2
    for k in range(0,N2):
        file_.write(str(k*dk+N1*0.001+N0*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*dk+N1*0.001+N0*0.0001)-0.126*(k*dk+N1*0.001+N0*0.0001)-0.3516*(k*dk+N1*0.001+N0*0.0001)**2+0.2843*(k*dk+N1*0.001+N0*0.0001)**3-0.1015*(k*dk+N1*0.001+N0*0.0001)**4))+"\n")
    k=0
    dk2=0.0005
    X_m=0.89 # Modify this parameter to smooth the trailing edge ==> X_m=1 -> pointed shape | X_m<1 -> smoothing
    while k*dk2+N2*dk+N1*0.001+N0*0.0001 <= X_m:
        file_.write(str(k*dk2+N2*dk+N1*0.001+N0*0.0001)+" "+str(xx/20.0*(0.2969*sqrt(k*dk2+N2*dk+N1*0.001+N0*0.0001)-0.126*(k*dk2+N2*dk+N1*0.001+N0*0.0001)-0.3516*(k*dk2+N2*dk+N1*0.001+N0*0.0001)**2+0.2843*(k*dk2+N2*dk+N1*0.001+N0*0.0001)**3-0.1015*(k*dk2+N2*dk+N1*0.001+N0*0.0001)**4))+"\n")
        k+=1
    Y_m=xx/20.0*(0.2969*sqrt(X_m)-0.126*X_m-0.3516*(X_m)**2+0.2843*(X_m)**3-0.1015*(X_m)**4)
    M=20
    for i in range(0,M+1):
        theta=i*math.pi/(2*M+1)
        file_.write(str(X_m+Y_m*sin(theta))+" "+str(Y_m*cos(theta))+"\n")
    file_.close()
    Table0=np.fromfile(File,sep=" ")
    N_len=len(Table0)
    file_=open(File,'a')
    for i in range(0,int(N_len/2.0)):
        file_.write(str(Table0[N_len-2*i-2])+" "+str(-1.0*Table0[N_len-2*i-1])+"\n")
    file_.close()

Table=np.fromfile(File,sep=" ")
N=len(Table)
Table_mem=[]

# Homothety
for i in range(0,N):
    Table[i]*=Ratio
    Table_mem.append(Table[i])

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
X_c=X_r*(X_max-X_min) 
Y_c=0

# Rotation
for i in range(0,int(N/2.0)):
    XX,YY=Table[2*i]-X_c,Table[2*i+1]-Y_c
    Table[2*i],Table[2*i+1]=X_c+cos(Rotation)*XX+sin(Rotation)*YY,Y_c-sin(Rotation)*XX+cos(Rotation)*YY

print("Writing process of NACA.geo")
file=open('NACA.geo','w')

for i in range(0,int(len(Table)/2)):
	file.write("Point("+str(i+1)+") = {"+str(Table[2*i])+", "+str(Table[2*i+1])+", "+str(Longueur)+", 1.0};\n")
	
file.write("\n")

for j in range(0,int(N/2)-1):
	file.write("Line("+str(j+1)+") = {"+str(j+1)+", "+str(j+2)+"};\n")
file.write("Line("+str(int(N/2))+") = {"+str(int(N/2))+", 1};\n")

# Wing tip
def indice(list_,element):
	r=1e5
	for i in range(0,len(list_)):
		d=abs(list_[i]-element)
	        if d<r:
		        r=d
		        ind=i
	return ind

if Saumon=='Circle':
	List_y=[Table_mem[2*i+1] for i in range(0,int(N/2.0))]
	l=0
	Line_loop_saumon=[]
	Line_loop_saumon_first=[]
	Nb_pts=15
	for i in range(0,int(N/4.0)):
		Pts1_x,Pts1_y=Table_mem[2*i],Table_mem[2*i+1]
		ind=indice(List_y,-1*Table_mem[2*i+1])
		Pts2_x,Pts2_y=Table_mem[2*ind],Table_mem[2*ind+1]
		R=0.5*sqrt((Pts1_x-Pts2_x)**2+(Pts1_y-Pts2_y)**2)
		if R>0.000001:
			for k in range(1,Nb_pts):
				x=X_c+cos(Rotation)*(Pts1_x-X_c)+sin(Rotation)*(R*cos(k*np.pi/Nb_pts)-Y_c)
				y=Y_c-sin(Rotation)*(Pts1_x-X_c)+cos(Rotation)*(R*cos(k*np.pi/Nb_pts)-Y_c)
				file.write("Point("+str(int(N/2)+(Nb_pts-1)*l+k)+") = {"+str(x)+", "+str(y)+", "+str(Longueur+R*sin(k*np.pi/Nb_pts))+", 1.0};\n")
			file.write("Line("+str(int(N/2)+Nb_pts*l+1)+") = {"+str(i+1)+", "+str(int(N/2)+(Nb_pts-1)*l+1)+"};\n")
			for k in range(1,Nb_pts-1):
				file.write("Line("+str(int(N/2)+Nb_pts*l+k+1)+") = {"+str(int(N/2)+(Nb_pts-1)*l+k)+", "+str(int(N/2)+(Nb_pts-1)*l+k+1)+"};\n")
			file.write("Line("+str(int(N/2)+Nb_pts*l+Nb_pts)+") = {"+str(int(N/2)+(Nb_pts-1)*l+Nb_pts-1)+", "+str(ind+1)+"};\n")
			l+=1
			if l==1:
				Line_loop_saumon.append(str(-i-1))
				for k in range(0,Nb_pts):
					Line_loop_saumon.append(str(int(N/2)+Nb_pts*(l-1)+k+1))
					Line_loop_saumon_first.append(str(int(N/2)+Nb_pts*(l-1)+k+1))
			else:
				if l>2:
					Line_loop_saumon.append(str(-i))
					for k in range(0,Nb_pts):
						Line_loop_saumon.append(str(int(N/2)+Nb_pts*(l-2)+k+1))
				Line_loop_saumon.append(str(-ind-1))
				for k in range(Nb_pts-1,-1,-1):
					Line_loop_saumon.append(str(-(int(N/2)+Nb_pts*(l-1)+k+1)))
				file.write("Line Loop("+str(l)+") = {")
				for j in range(0,len(Line_loop_saumon)):
					file.write(str(Line_loop_saumon[j])+",")
				file.seek(-1, os.SEEK_END)
				file.write("};\n")
				file.write("Plane Surface("+str(l)+") = {"+str(l)+"};\n")
				Line_loop_saumon=[]
		else:
			pass
	file.write("Line Loop("+str(l+1)+") = {1,")
	for k in range(0,len(Line_loop_saumon_first)):
		file.write(str(Line_loop_saumon_first[k])+",")
	file.write(str(int(N/2)-1)+","+str(int(N/2))+"};\nPlane Surface("+str(l+1)+") = {"+str(l+1)+"};\n")
else:
	pass

file.write("Line Loop(1) = {")
for i in range(1,int(N/2)):
	file.write(str(i)+", ")
file.write(str(int(N/2))+"};\n")
file.write("\nPlane Surface(1) = {1};\n")
file.write("Extrude {0, 0, "+str(-1*Longueur)+"} {\n Surface{1};\nLayers{")
if twoD==False:
	file.write(str(int(Longueur*nb_discr)))
else:
	file.write("1")
file.write("};\nRecombine;\n}\n")

file.close()
print("End of the writing process of NACA.geo")
print("Coarse mesh of the wing with GMSH")
executable=str(Chemin_GMSH)+' NACA.geo -3 -format stl'
os.system(executable) 
print("End of the meshing process")
