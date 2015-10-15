import os,sys
import subprocess
import numpy as np



springfile=open("springs.txt","r")
coordfile=open("shearcoordinates.txt","r")
energyfile=open("shearenergy.txt","r")


springs=[]
coords=[]
boxdx=[]


for line in springfile:
	springs.append(line.split())

springs=[springs[i][0:4] for i in range(0,len(springs))]
for i in range(0,len(springs)):
	for j in range(0,4):
		springs[i][j]=int(springs[i][j])


for line in energyfile:
		line=line.split()
		boxdx.append(float(line[0]))


lineNR=-1;

for line in coordfile:
	coords=[]
	lineNR=lineNR+1
	line=line.split()
	for j in range(0,len(line)):
		coords.append(float(line[j]))

	numbercoords=len(coords)/2
	numbersprings=len(springs)

	phi=boxdx[lineNR]
	e1x=1
	e1y=0
	e2x=np.tan(phi)
	e2y=1

	plotcoordfile=open("plotcoords.txt","w")
	for k in range(0,len(springs)):
		for gamma in range(-1,2):
			for delta in range(-1,2):
				one=springs[k][0]
				two=springs[k][1]	
				wlr=springs[k][2]
				wud=springs[k][3]
				y1=coords[one+numbercoords]+delta*e2y+gamma*e1y
				y2=coords[two+numbercoords]+wud+delta*e2y+gamma*e1y
				y1b=coords[one+numbercoords]
				y2b=coords[two+numbercoords]+wud	
				x1=coords[one]+gamma*e1x+delta*e2x+np.tan(phi)*y1b
				x2=coords[two]+wlr+gamma*e1x+delta*e2x+np.tan(phi)*y2b
				plotcoordfile.write(str(x1)+"\t"+str(y1)+"\n")
				plotcoordfile.write(str(x2)+"\t"+str(y2)+"\n\n")

	plotcoordfile.close()

	boxfile=open("boxfile.txt","w")
	boxfile.write(str(0)+"\t"+str(0)+"\n")
	boxfile.write(str(1)+"\t"+str(0)+"\n\n")
	boxfile.write(str(1)+"\t"+str(0)+"\n")
	boxfile.write(str(1+np.tan(phi))+"\t"+str(1)+"\n\n")
	boxfile.write(str(1+np.tan(phi))+"\t"+str(1)+"\n")
	boxfile.write(str(np.tan(phi))+"\t"+str(1)+"\n\n")

	boxfile.write(str(np.tan(phi))+"\t"+str(1)+"\n")
	boxfile.write(str(0)+"\t"+str(0)+"\n")
	boxfile.close()


	plotfile=open("plotfile.txt","w")
	plotfile.write("set term png\n")
	plotfile.write("set yrange[-.1 : 1.1]\n")
	plotfile.write("set xrange [-.1:1.1]\n")
	plotfile.write("set size square\n")

	plotfile.write("set output "+"'movie/"+str(lineNR)+".png'\n")

	plotfile.write("plot 'plotcoords.txt' with lines lt black notitle, 'boxfile.txt' with lines lw 3 notitle\n") 

	plotfile.write("set size square\n")
#	plotfile.write("replot\n")
	plotfile.close()

	os.system("gnuplot 'plotfile.txt'")






