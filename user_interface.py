#!/usr/bin/env python
import Tkinter
import os
from Tkconstants import *

#Copyright (C) 2017					
#Edison Montoya						
#eamonto@gmail.com
#Up to date: Feb 16  2017				

def calculate(*args):
    try:
        f = open("param.txt","w")

        f.write(input_box[0].get()+"\n")
        f.write(input_box[1].get()+"\n")
        f.write(input_box[2].get()+"\n")
        f.write(input_box[3].get()+"\n")
        f.write(input_box[4].get()+"\n")
        f.write(input_box[5].get()+"\n")
#        f.write(input_box[6].get()+"\n")
#        f.write(input_box[7].get()+"\n")
#        f.write(input_box[8].get()+"\n")
        f.write(input_box[9].get()+"\n")
        f.write(input_box[10].get()+"\n")
        f.write(input_box[11].get()+"\n")
        f.write(input_box[12].get()+"\n")

#        f.write(input_box2[0].get()+"\n")
        f.write(input_box2[1].get()+"\n")
        f.write(input_box2[2].get()+"\n")
#        f.write(input_box2[3].get()+"\n")
        f.write(input_box2[4].get()+"\n")
        f.write(input_box2[5].get()+"\n")
#        f.write(input_box2[6].get()+"\n")
        f.write(input_box2[7].get()+"\n")
#        f.write(input_box2[8].get()+"\n")
#        f.write(input_box2[9].get()+"\n")

        f.close()

        print os.system("./exec.out param.txt "+input_box[13].get())

    except ValueError:
        pass


tk = Tkinter.Tk()

#Window Tittle
tk.title("Input Parameters")

#Principal Frame
frame = Tkinter.Frame(tk, relief=RIDGE, borderwidth=2)
frame.grid(column=0, row=0, sticky=(N, W, E, S))

frame2 = Tkinter.Frame(tk, relief=RIDGE, borderwidth=2)
frame2.grid(column=0, row=1, sticky=(N, W, E, S))


#io variables 
l = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
input_box = [Tkinter.StringVar() for n in l]
l2 = [0,1,2,3,4,5,6,7,8,9]
input_box2 = [Tkinter.StringVar() for n in l2]


#Input box
box_entry=[Tkinter.Entry(frame, width=7, textvariable=input_box[n]) for n in l]


#Tkinter.Radiobutton(frame2, text="RK4",   variable=input_box2[0], value="0").grid(column=2, row=0, sticky=W)
#Tkinter.Radiobutton(frame2, text="RK-F",  variable=input_box2[0], value="1").grid(column=3, row=0, sticky=W)
#Tkinter.Radiobutton(frame2, text="RK-CK", variable=input_box2[0], value="2").grid(column=4, row=0, sticky=W)

#Tkinter.Radiobutton(frame2, text="k=1", variable=input_box2[1], value="0").grid(column=2, row=1, sticky=W)
Tkinter.Radiobutton(frame2, text="Bianchi I",   variable=input_box2[1], value="1").grid(column=2, row=1, sticky=W)
Tkinter.Radiobutton(frame2, text="Bianchi IX",  variable=input_box2[1], value="2").grid(column=3, row=1, sticky=W)

Tkinter.Radiobutton(frame2, text="Classical",variable=input_box2[2], value="0").grid(column=2, row=2, sticky=W)
Tkinter.Radiobutton(frame2, text="Effective",variable=input_box2[2], value="1").grid(column=3, row=2, sticky=W)

#Tkinter.Radiobutton(frame2, text="(muici,p_i)",  variable=input_box2[3], value="0").grid(column=2, row=3, sticky=W)
#Tkinter.Radiobutton(frame2, text="(mu1c1,mu2c2,P_phi,p_i)",variable=input_box2[3], value="1").grid(column=3, row=3, sticky=W)

Tkinter.Radiobutton(frame2, text="Edward (N=1)",variable=input_box2[4], value="0").grid(column=2, row=4, sticky=W)
Tkinter.Radiobutton(frame2, text="Asieh (N=V)", variable=input_box2[4], value="1").grid(column=3, row=4, sticky=W)
#Tkinter.Radiobutton(frame2, text="B-S",   variable=input_box2[4], value="2").grid(column=4, row=4, sticky=W)

Tkinter.Radiobutton(frame2, text="not",variable=input_box2[5], value="0").grid(column=2, row=5, sticky=W)
Tkinter.Radiobutton(frame2, text="yes",variable=input_box2[5], value="1").grid(column=3, row=5, sticky=W)

#Tkinter.Radiobutton(frame2, text="Internal (N=V)",variable=input_box2[6], value="0").grid(column=2, row=6, sticky=W)
#Tkinter.Radiobutton(frame2, text="Cosmic (N=1)",  variable=input_box2[6], value="1").grid(column=3, row=6, sticky=W)

Tkinter.Radiobutton(frame2, text="Off",  variable=input_box2[7], value="0").grid(column=2, row=7, sticky=W)
Tkinter.Radiobutton(frame2, text="On",   variable=input_box2[7], value="1").grid(column=3, row=7, sticky=W)

#Tkinter.Radiobutton(frame2, text="Off", variable=input_box2[8], value="0").grid(column=2, row=8, sticky=W)
#Tkinter.Radiobutton(frame2, text="On",  variable=input_box2[8], value="1").grid(column=3, row=8, sticky=W)

#Tkinter.Radiobutton(frame2, text="Inflation", variable=input_box2[9], value="0").grid(column=2, row=9, sticky=W)
#Tkinter.Radiobutton(frame2, text="Moduli (Cyclic)",    variable=input_box2[9], value="1").grid(column=3, row=9, sticky=W)


#Text Labels
Tkinter.Label(frame, text="mu1c1 =").grid( column=1, row=0, sticky=E)
Tkinter.Label(frame, text="mu2c2 =").grid( column=1, row=1, sticky=E)
Tkinter.Label(frame, text="mu3c3 =").grid( column=1, row=2, sticky=E)
Tkinter.Label(frame, text="p1 =").grid( column=1, row=3, sticky=E)
Tkinter.Label(frame, text="p2 =").grid( column=1, row=4, sticky=E)
Tkinter.Label(frame, text="p3 =").grid( column=1, row=5, sticky=E)
#Tkinter.Label(frame, text="phi =").grid(column=1, row=6, sticky=E)
#Tkinter.Label(frame, text=" Tolerance error =").grid(column=1, row=7, sticky=E)
#Tkinter.Label(frame, text="Max iterations =").grid( column=1, row=8, sticky=E)
Tkinter.Label(frame, text="Initial time =").grid(   column=1, row=9, sticky=E)
Tkinter.Label(frame, text="Final time =").grid(     column=1, row=10, sticky=E)
Tkinter.Label(frame, text="dt =").grid(             column=1, row=11, sticky=E)
Tkinter.Label(frame, text="Output every ").grid(    column=1, row=12, sticky=E)
Tkinter.Label(frame, text="Output file ").grid(     column=1, row=13, sticky=E)
Tkinter.Label(frame, text="                                     ").grid(column=2, row=14, sticky=E)


#Tkinter.Label(frame2, text="Integrator").grid(column=1, row=0, sticky=E)
Tkinter.Label(frame2, text="Equations").grid(column=1, row=1, sticky=E)
Tkinter.Label(frame2, text="Level").grid(column=1, row=2, sticky=E)
#Tkinter.Label(frame2, text="Initial Conditions").grid(column=1, row=3, sticky=E)
Tkinter.Label(frame2, text="Effective Equations").grid(column=1, row=4, sticky=E)
Tkinter.Label(frame2, text="Effective to FRW").grid(column=1, row=5, sticky=E)
#Tkinter.Label(frame2, text="Time").grid(column=1, row=6, sticky=E)
Tkinter.Label(frame2, text=" Standard output").grid(column=1, row=7, sticky=E)
#Tkinter.Label(frame2, text="Potential switch").grid(column=1, row=8, sticky=E)
#Tkinter.Label(frame2, text="Select potential").grid(column=1, row=9, sticky=E)


#Button
Tkinter.Button(frame2, text="Calculate", command=calculate).grid(column=2, row=8, sticky=E)


for n in l:
    if n!=6 and n!=7 and n!=8:
        box_entry[n].grid(column=2, row=n, sticky=(W, E))

#Initial Values
input_box[0].set(1.57079633)
input_box[1].set(1.57079633)
input_box[2].set(1.57079633)
input_box[3].set(100.0)
input_box[4].set(200.0)
input_box[5].set(300.0)
#input_box[6].set(0.0)
#input_box[7].set(1.0e-7)
#input_box[8].set(100)
input_box[9].set(0.0)
input_box[10].set(-3.5)
input_box[11].set(-5.0e-4)
input_box[12].set(20)
input_box[13].set("output")

#input_box2[0].set(0)
input_box2[1].set(2)
input_box2[2].set(1)
#input_box2[3].set(0)
input_box2[4].set(0)
input_box2[5].set(0)
#input_box2[6].set(1)
input_box2[7].set(0)
#input_box2[8].set(0)
#input_box2[9].set(0)


#Put the cursor on the input box
box_entry[0].focus()

#If the user press Enter or the button then call the calculate routine
tk.bind('<Return>', calculate)

#Run the window
tk.mainloop()
