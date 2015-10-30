import os,sys
import numpy as np

mikadofile=open("mikado.txt","r")
paramfile=open("params.txt","r")

mikado=[]

for line in paramfile:

	line=line.split('=')
	if line[0]=="<Lstick":
		line[1]=line[1].replace(">","0") #replace the bracket by a zero
		Lstick=float(line[1])
	


for line in mikadofile:
	line=line.split()
	mikado.append(line)

	
plotmikado=open("plotmikado.txt","w")
for i in range(0,len(mikado)):
	#convert to int and floats
	mikado[i][0]=int(mikado[i][0])
	mikado[i][1]=float(mikado[i][1])
	mikado[i][2]=float(mikado[i][2])
	mikado[i][3]=float(mikado[i][3])
	mikado[i][4]=int(mikado[i][4])
	mikado[i][5]=int(mikado[i][5])
	th=mikado[i][3]	
	x0=mikado[i][1]
	y0=mikado[i][2]
	x1=x0+Lstick*np.cos(th);
	y1=y0+Lstick*np.sin(th);
	#write the data to txtfile in all 9 box-images
	plotmikado.write(str(x0)+"\t"+str(y0)+"\n")
	plotmikado.write(str(x1)+"\t"+str(y1)+"\n\n")
	
	plotmikado.write(str(x0+1)+"\t"+str(y0)+"\n")
	plotmikado.write(str(x1+1)+"\t"+str(y1)+"\n\n")
	
	plotmikado.write(str(x0-1)+"\t"+str(y0)+"\n")
	plotmikado.write(str(x1-1)+"\t"+str(y1)+"\n\n")

	plotmikado.write(str(x0)+"\t"+str(y0+1)+"\n")
	plotmikado.write(str(x1)+"\t"+str(y1+1)+"\n\n")

	plotmikado.write(str(x0)+"\t"+str(y0-1)+"\n")
	plotmikado.write(str(x1)+"\t"+str(y1-1)+"\n\n")

	plotmikado.write(str(x0+1)+"\t"+str(y0-1)+"\n")
	plotmikado.write(str(x1+1)+"\t"+str(y1-1)+"\n\n")

	plotmikado.write(str(x0+1)+"\t"+str(y0+1)+"\n")
	plotmikado.write(str(x1+1)+"\t"+str(y1+1)+"\n\n")

	plotmikado.write(str(x0-1)+"\t"+str(y0-1)+"\n")
	plotmikado.write(str(x1-1)+"\t"+str(y1-1)+"\n\n")

	plotmikado.write(str(x0-1)+"\t"+str(y0+1)+"\n")
	plotmikado.write(str(x1-1)+"\t"+str(y1+1)+"\n\n")
	
	

	
		
		

plotmikado.close()

periodicnodes=[]	
periodicnodefile=open("pnodes.txt","w")
nodefile=open("nodes.txt","r")
for line in nodefile:
	line=line.split()
	line[0]=int(line[0])
	line[1]=float(line[1])
	line[2]=float(line[2])
	periodicnodefile.write(str(line[0])+"\t"+str(line[1]+1)+"\t"+str(line[2])+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1]-1)+"\t"+str(line[2])+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1])+"\t"+str(line[2]+1)+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1])+"\t"+str(line[2]-1)+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1]+1)+"\t"+str(line[2]+1)+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1]+1)+"\t"+str(line[2]-1)+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1]-1)+"\t"+str(line[2]+1)+"\n")
	periodicnodefile.write(str(line[0])+"\t"+str(line[1]-1)+"\t"+str(line[2]-1)+"\n")
	

box=open("boxinitial.txt","w")
box.write(str(0)+"\t"+str(0)+"\n")
box.write(str(0)+"\t"+str(1)+"\n\n")
box.write(str(0)+"\t"+str(1)+"\n")
box.write(str(1)+"\t"+str(1)+"\n\n")
box.write(str(1)+"\t"+str(1)+"\n")
box.write(str(1)+"\t"+str(0)+"\n\n")
box.write(str(1)+"\t"+str(0)+"\n")
box.write(str(0)+"\t"+str(0)+"\n\n")
box.close()

#make sure everything is plotted in one file
gnuplotfile=open("plotmikadoperiodicwbox.txt","w")
gnuplotfile.write("set term eps\n")
gnuplotfile.write("set xra[-.1:1.1]\n")
gnuplotfile.write("set yra[-.1:1.1]\n")
gnuplotfile.write("set size square\n")
gnuplotfile.write("set output 'output/pltmikado.eps'\n")

gnuplotfile.write("plot 'plotmikado.txt' with lines lt black notitle , ")
gnuplotfile.write("'nodes.txt' using 2:3 pt 7 ps 0.3 lt rgb 'blue' notitle , ")
gnuplotfile.write("'pnodes.txt' using 2:3 pt 7 ps 0.3 lt rgb 'blue' notitle , \n")
#gnuplotfile.write("'boxinitial.txt' with lines linecolor rgb '#09ad00' lw 3 notitle"+"\n\n")
#gnuplotfile.write("set xrange[-.1 : 1.1] \nset yrange[-.1 : 1.1]"+"\n")
gnuplotfile.write("set size square \n")
#gnuplotfile.write("replot")
gnuplotfile.close()


os.system("gnuplot plotmikadoperiodicwbox.txt")





