import numpy as np
from numpy import *
import sys
import os

File=sys.argv[1]
Ratio=float(sys.argv[2])
Envergure=float(sys.argv[3])
Length=float(sys.argv[4])
Width=float(sys.argv[5])
X_r=float(sys.argv[6])
Rotation=float(sys.argv[7])*3.141592/180.0
xx=int(sys.argv[8])

# Writing process of the NACA.dat file
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

# Homotethie
for i in range(0,N):
    Table[i]*=Ratio

# Acquisition du centre de rotation
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
X_c=X_r*(X_max-X_min) # Car le bord d'attaque et le bord de fuite sont inverses par rapport aux coordonnees selon x
Y_c=X_r*(Y_max-Y_min)

# Rotation
for i in range(0,int(N/2.0)):
    XX,YY=Table[2*i]-X_c,Table[2*i+1]-Y_c
    Table[2*i],Table[2*i+1]=X_c+cos(Rotation)*XX+sin(Rotation)*YY,Y_c-sin(Rotation)*XX+cos(Rotation)*YY

# Ecriture du fichier ftr
fichier=open("NACA.ftr","w")
fichier.write("4\n(\nAile\nwall\n\nInlet\nwall\n\nOutlet\nwall\n\nWall\nwall\n)\n\n"+str(N+8)+"\n(\n")
for i in range(0,int(N/2.0)):
    fichier.write("("+str(Table[2*i])+" "+str(Table[2*i+1])+" 0)\n")
fichier.write("\n")
for i in range(0,int(N/2.0)):
    fichier.write("("+str(Table[2*i])+" "+str(Table[2*i+1])+" "+str(Envergure)+")\n")

CDG_x=0
CDG_y=0
for i in range(0,int(N/2.0)):
    CDG_x+=2*Table[2*i]/N
    CDG_y+=2*Table[2*i+1]/N

fichier.write("\n("+str(CDG_x-Length/2.0)+" "+str(CDG_y-Width/2.0)+" 0)\n") # numero N
fichier.write("("+str(CDG_x+Length/2.0)+" "+str(CDG_y-Width/2.0)+" 0)\n") # N+1
fichier.write("("+str(CDG_x-Length/2.0)+" "+str(CDG_y+Width/2.0)+" 0)\n") # N+2
fichier.write("("+str(CDG_x+Length/2.0)+" "+str(CDG_y+Width/2.0)+" 0)\n") # N+3
fichier.write("("+str(CDG_x-Length/2.0)+" "+str(CDG_y-Width/2.0)+" "+str(Envergure)+")\n") #N+4
fichier.write("("+str(CDG_x+Length/2.0)+" "+str(CDG_y-Width/2.0)+" "+str(Envergure)+")\n") #N+5
fichier.write("("+str(CDG_x-Length/2.0)+" "+str(CDG_y+Width/2.0)+" "+str(Envergure)+")\n") #N+6
fichier.write("("+str(CDG_x+Length/2.0)+" "+str(CDG_y+Width/2.0)+" "+str(Envergure)+")\n") #N+7

fichier.write(")\n\n"+str(len(Table)+8)+"\n(\n")

for i in range(0,int(N/2.0)-1):
    fichier.write('(('+str(i)+' '+str(i+1)+' '+str(i+int(N/2)+1)+') 0)\n')
    fichier.write('(('+str(i)+' '+str(i+int(N/2)+1)+' '+str(i+int(N/2))+') 0)\n')
fichier.write('(('+str(int(N/2)-1)+' '+str(0)+' '+str(int(N/2))+') 0)\n')
fichier.write('(('+str(int(N/2)-1)+' '+str(int(N/2))+' '+str(N-1)+') 0)\n')

# Bloc de fluide autour de l'aile
fichier.write('(('+str(N+2)+' '+str(N+4)+' '+str(N)+') 1)\n')
fichier.write('(('+str(N+2)+' '+str(N+6)+' '+str(N+4)+') 1)\n')
fichier.write('(('+str(N+3)+' '+str(N+5)+' '+str(N+1)+') 2)\n')
fichier.write('(('+str(N+3)+' '+str(N+7)+' '+str(N+5)+') 2)\n')
fichier.write('(('+str(N+3)+' '+str(N+6)+' '+str(N+2)+') 3)\n')
fichier.write('(('+str(N+3)+' '+str(N+7)+' '+str(N+6)+') 3)\n')
fichier.write('(('+str(N+1)+' '+str(N+4)+' '+str(N)+') 3)\n')
fichier.write('(('+str(N+1)+' '+str(N+5)+' '+str(N+4)+') 3)\n')
fichier.write(')')
fichier.close()
