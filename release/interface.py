from Tkinter import *
from PIL import Image, ImageTk
import os,sys
#import Tkinter.messagebox

def store():
	params=open("params.txt","w")
	if(rlengshort.get()==""):
		rlengshort.set(".0001")
	params.write("###Network parameters###\n")
	params.write("<rlenshort="+rlengshort.get()+">\n")
	params.write("<rlenlong="+rlenglong.get()+">\n")
	params.write("<stretchf="+stretchf.get()+">\n\n")
		
	params.write("<bendingOn="+str(bendingon.get())+">\n")
	if(bendingon.get()==0):
		kappa.set("0.0")
	params.write("<kappa="+kappa.get()+">\n\n")
	
	params.write("<k1="+K1.get()+">\n")
	params.write("<k2="+K2.get()+">\n")
	
	params.write("<Lstick="+length.get()+">\n")
	params.write("<NumberMikado="+number.get()+">\n")
	
	params.write("###Background Parameters###\n")
	params.write("<backGroundOn="+str(backgroundtype.get())+">\n")
	params.write("TypeOn=0 means off, 1 for a square lattice, 2 for a hexagonal lattice\n\n")
	
	params.write("###Conjugate gradient parameters###\n")
	params.write("<Nit="+Iteration.get()+">\n")
	params.write("<tolGradE="+LenGrad.get()+">\n\n")
	
	params.write("###Shear-parameters###\n")
	params.write("<StepSize="+StepSize.get()+">\n")
	params.write("<NumberStepsRight="+RightSteps.get()+">\n")
	params.write("<NumberStepsLeft="+LeftSteps.get()+">\n\n")
	
	params.write("<END>")

	params.close()
	return

def testrun():	
	params=open("params.txt","w")
	if(rlengshort.get()==""):
		rlengshort.set(".0001")
	params.write("###Network parameters###\n")
	params.write("<rlenshort="+rlengshort.get()+">\n")
	params.write("<rlenlong="+rlenglong.get()+">\n")
	params.write("<stretchf="+stretchf.get()+">\n\n")
		
	params.write("<bendingOn="+str(bendingon.get())+">\n")
	if(bendingon.get()==0):
		kappa.set("0.0")
	params.write("<kappa="+kappa.get()+">\n\n")
	
	params.write("<k1="+K1.get()+">\n")
	params.write("<k2="+K2.get()+">\n")
	
	params.write("<Lstick="+length.get()+">\n")
	params.write("<NumberMikado="+number.get()+">\n")
	
	params.write("###Background Parameters###\n")
	params.write("<backGroundOn="+str(backgroundtype.get())+">\n")
	params.write("TypeOn=0 means off, 1 for a square lattice, 2 for a hexagonal lattice\n\n")
	
	params.write("###Conjugate gradient parameters###\n")
	params.write("<Nit="+Iteration.get()+">\n")
	params.write("<tolGradE="+LenGrad.get()+">\n\n")
	
	params.write("###Shear-parameters###\n")
	params.write("<StepSize="+StepSize.get()+">\n")
	params.write("<NumberStepsRight="+str(1)+">\n")
	params.write("<NumberStepsLeft="+str(0)+">\n\n")
	
	params.write("<END>")

	params.close()
	
	os.system("./mikado")
	os.system("python2 plotmikado.py")
	
	image=Image.open("output/pltmikado.eps")
	
	image=image.resize((400,400),Image.ANTIALIAS)
	photo=ImageTk.PhotoImage(image)
	label=Label(image=photo)
	label.image=photo
	label.place(relx=1,x=2,y=2,anchor=NE)	
	return





app=Tk()


app.title("ImportStuff")
app.geometry('1200x500+200+200')

labelText=StringVar()
labelText.set("Rest Length")
label1=Label(app,textvariable=labelText,height=1)
label1.grid(row=1, column=0)
rlengshort=StringVar()
rlengshort.set("0.1")
RESTLENGHT=Entry(app, textvariable=rlengshort)
RESTLENGHT.grid(row=1,column=1)

labelText2=StringVar()
labelText2.set("Rest Length 2")
label2=Label(app,textvariable=labelText2,height=1)
label2.grid(row=2,column=0)
rlenglong=StringVar()
rlenglong.set("0.01")
RESTLENGTH2=Entry(app,textvariable=rlenglong)
RESTLENGTH2.grid(row=2,column=1)

stretchText=StringVar()
stretchText.set("Stretching Factor")
stretchlabel=Label(app,textvariable=stretchText,height=1)
stretchlabel.grid(row=3,column=0)
stretchf=StringVar()
stretchf.set("1.0")
STRETCHF=Entry(app,textvariable=stretchf)
STRETCHF.grid(row=3,column=1)

Bendingontekst=StringVar()
Bendingontekst.set("Bending on / off")
Bendignonlabel=Label(app,textvariable=Bendingontekst)
Bendignonlabel.grid(row=5,column=0)

bendingon=IntVar()
bendingradio1=Radiobutton(app,text="off",variable=bendingon,value=0)
bendingradio1.grid(row=5,column=1)
bendingradio2=Radiobutton(app,text="on",variable=bendingon,value=1)
bendingradio2.grid(row=5,column=2)

kappatekst=StringVar()
kappatekst.set("kappa")
kappalabel=Label(app,textvariable=kappatekst,height=0)
kappalabel.grid(row=6,column=0)
kappa=StringVar()
kappa.set("1e-9")
KAPPA=Entry(app,textvariable=kappa)
KAPPA.grid(row=6,column=1)

k1text=StringVar()
k1text.set("spring constant 1 ")
k1label=Label(app,textvariable=k1text,height=0)
k1label.grid(row=7,column=0)
k1=StringVar()
k1.set("0.5")
K1=Entry(app,textvariable=k1)
K1.grid(row=7,column=1)

k2text=StringVar()
k2text.set("spring constant 2 ")
k2label=Label(app,textvariable=k2text,height=0)
k2label.grid(row=8,column=0)
k2=StringVar()
k2.set("0.5")
K2=Entry(app,textvariable=k2)
K2.grid(row=8,column=1)

numbertext=StringVar()
numbertext.set("# of sticks")
numberlabel=Label(app,textvariable=numbertext,height=0)
numberlabel.grid(row=10,column=0)
number=StringVar()
number.set("100")
numberentry=Entry(app,textvariable=number)
numberentry.grid(row=10,column=1)

lentext=StringVar()
lentext.set("Length of the sticks")
lenlabel=Label(app,textvariable=lentext,height=0)
lenlabel.grid(row=11,column=0)
length=StringVar()
length.set("0.4")
lenentry=Entry(app,textvariable=length)
lenentry.grid(row=11,column=1)




backgroundontext=StringVar()
backgroundontext.set("Background Network")
backgroundlabel=Label(app,textvariable=backgroundontext,height=0)
backgroundlabel.grid(row=12,column=0)
backgroundtype=IntVar()
backgroundtype.set(0)
backgroundradio1=Radiobutton(app,text="off",variable=backgroundtype,value=0)
backgroundradio1.grid(row=12,column=1)
backgroundradio2=Radiobutton(app,text="Square",variable=backgroundtype,value=1)
backgroundradio2.grid(row=12,column=2)
backgroundradio3=Radiobutton(app,text="Hexagonal",variable=backgroundtype,value=2)
backgroundradio3.grid(row=12,column=3)

Iterationtext=StringVar()
Iterationtext.set("Max Nr of CG Iterations")
Iterationlabel=Label(app,textvariable=Iterationtext,height=0)
Iterationlabel.grid(row=13,column=0)
Iteration=StringVar()
Iteration.set("1e13")
Iterationentry=Entry(app,textvariable=Iteration)
Iterationentry.grid(row=13,column=1)

LenGradtext=StringVar()
LenGradtext.set("Length of the Gradient")
LenGradLabel=Label(app,textvariable=LenGradtext,height=0)
LenGradLabel.grid(row=14,column=0)
LenGrad=StringVar()
LenGrad.set("1e-16")
LenGradentry=Entry(app,textvariable=LenGrad)
LenGradentry.grid(row=14,column=1)

StepSizetext=StringVar()
StepSizetext.set("Size of the shearstep")
StepSizeLabel=Label(app,textvariable=StepSizetext,height=0)
StepSizeLabel.grid(row=15,column=0)
StepSize=StringVar()
StepSize.set("0.0001")
StepSizeEntry=Entry(app,textvariable=StepSize)
StepSizeEntry.grid(row=15,column=1)


RightSteptext=StringVar()
RightSteptext.set("Number of steps to the right")
RightStepLabel=Label(app,textvariable=RightSteptext,height=0)
RightStepLabel.grid(row=16,column=0)
RightSteps=StringVar()
RightSteps.set("100")
RightStepsEntry=Entry(app,textvariable=RightSteps)
RightStepsEntry.grid(row=16,column=1)


LeftSteptext=StringVar()
LeftSteptext.set("Number of steps to the left")
LeftStepLabel=Label(app,textvariable=LeftSteptext,height=0)
LeftStepLabel.grid(row=17,column=0)
LeftSteps=StringVar()
LeftSteps.set("200")
LeftStepsEntry=Entry(app,textvariable=LeftSteps)
LeftStepsEntry.grid(row=17,column=1)



# Here are the buttons 
store=Button(app,text="Save",width=20,command=store)
store.grid(row=19,column=1)

TestRun=Button(app,text="TestRun",width=20,command=testrun)
TestRun.grid(row=19,column=2)

app.mainloop()
