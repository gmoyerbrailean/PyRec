####### LICENSE #######
# This code is part of the Recombineering module, written by Gregory
# Moyerbrailean at Michigan State University, Department of Microbiology 
# and Molecular Genetics. 
# Copyright (C) 2010 Gregory Moyerbrailean
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

'''A graphical interface for the PyRec program'''

from Tkinter import *
from recombineering import *
import tkMessageBox
import os

class PyRec:
    def __init__(self, parent=0):
        self.F = Frame(parent)
        self.F.master.title("PyRec")

        ## Rows 0-2: Header Image
        fHeader = Frame(self.F).grid(row=0, rowspan=3, column=0, columnspan=8)

        # bkrd = PhotoImage(file="extras/DNA_backdrop.gif")
        bkrd = PhotoImage(file="Files/extras/DNA_backdrop.gif")
        header = Label(fHeader,image=bkrd)
        header.image = bkrd # keep as reference
        header.grid(row=0, rowspan=3, column=0, columnspan=8)

        ## Rows 3-5: User Input
        fOptions = Frame(self.F).grid(row=3, rowspan=4, column=0, columnspan=8)

        # Row 3: File Entry
        lFile = Label(fOptions, text="File:").grid(row=3,column=0,sticky=E)
        self.eFile = Entry(fOptions,width=75)
        # self.eFile.insert(0,"genes/")
        self.eFile.insert(0, "genes/lactis/mutS.fasta")
        self.eFile.grid(row=3,column=1,columnspan=7,sticky=W)

        # Row 4: Organism, Report All Oligos
        lOrganism = Label(fOptions,text="Organism:").grid(row=4,column=0,
                                                          sticky=W)

        ## Automatically detect genomes in the genome/ directory
        orgs = os.listdir('Files/genome/')
        orgs = [o for o in orgs if not o[0]=='.'] # remove hidden files
        self.vOrganism = StringVar()
        self.mOrganism = apply(OptionMenu, (fOptions, self.vOrganism) +
            tuple(orgs))
        self.mOrganism.grid(row=4,column=1,columnspan=2,sticky=W)
        
        self.outputVar = IntVar()
        cbOutput = Checkbutton(fOptions, text="Report all possible oligos",
                               variable=self.outputVar)
        cbOutput.grid(row=4,column=3,columnspan=3,sticky=W)

        # Row 5: Begin, End
        lBegin = Label(fOptions,text="Begin:").grid(row=5,column=0,sticky=E)
        self.eBegin = Entry(fOptions)
        self.eBegin.insert(0,"40")
        self.eBegin.grid(row=5,column=1,sticky=W)

        lEnd = Label(fOptions,text="End:").grid(row=5
                                                ,column=2)
        self.eEnd = Entry(fOptions)
        self.eEnd.insert(0,"350")
        self.eEnd.grid(row=5,column=3,sticky=W)

        # Row 6: Cutsites
        lCuts = Label(fOptions,text="Cutsites:").grid(row=6,column=0,sticky=E)
        self.eCuts = Entry(fOptions,width=75)
        self.eCuts.insert(0,'GAATTC,GGATCC,AAGCTT')
        self.eCuts.grid(row=6,column=1,columnspan=7,sticky=W)

        ## Row 7: Buttons
        fButtons = Frame(self.F).grid(row=7, column=0, columnspan=8)

        # Run: Passes the inputs to the recombineering module
        bRun = Button(fButtons,text="Run!",command=self.Run)
        bRun.grid(row=7,column=1)#,columnspan=2,sticky=E)

        # Reset: Resets the entry fields
        bReset = Button(fButtons,text="Reset",command=self.Reset)
        bReset.grid(row=7,column=2)#,columnspan=2,sticky=W)

        # Depreciated as of version 2.5
##        # Oligos: Opens the directory of created oligos
##        bOligos = Button(fButtons,text="Oligos",command=self.OpenOligos)
##        bOligos.grid(row=6,column=4)#,columnspan=2,sticky=W)

        # Help: Opens the user manual file
        bHelp = Button(fButtons,text="Help",command=self.Help)
        bHelp.grid(row=7,column=3)#,columnspan=2,sticky=W)

        # Quit: Quits the program
        bQuit = Button(fButtons,text="Quit",command=self.Quit)
        bQuit.grid(row=7,column=4)

        ## Row 8-9: Footer
        fFooter = Frame(self.F).grid(row=9, column=0, columnspan=8)

        # Row 8: Spacer
        lSpace = Label(fFooter,text="").grid(row=8,column=0,columnspan=8)
        
        lFooter = Label(fFooter,text=\
                        "Designed by Gregory Moyerbrailean, 2010   ")
        lFooter.grid(row=9,column=0,columnspan=8,sticky=E)

    def Popup(self, message):
        tkMessageBox.showinfo("Warning!",message)

    def Run(self):

        # Create a control loop to kill the "get" if an error occurs
        while True:

            # Get the required inputs
            fileN = self.eFile.get()
            if not fileN:
                self.Popup("Filename cannot be blank!")
                break
            organism = self.vOrganism.get()
            if not organism:
                self.Popup("Must specify an organism!")
                break
            else:
                organism = organism.replace(' ','')

            cutsites = self.eCuts.get()
            if not cutsites:
                self.Popup("Must specify at least one cutsite!")
                break

            # Beginning and end have default values, use them if the entry
            # fields are blank (should not be, due to .insert method)
            beg = self.eBegin.get()
            if not beg:
                beg = 40
            elif beg.isdigit() and int(beg) > 0:
                beg = int(beg)
            else:
                self.Popup("'Begin' value must be a positive integer!")
                break

            end = self.eEnd.get()
            if not end:
                end = 350
            elif end.isdigit() and  int(end) > 0 and int(end) > int(beg):
                end = int(end)
            else:
                self.Popup("'End' value must be a positive integer larger than"\
                           " 'Begin' value!")
                break

            # Determine whether to restrict the output or not
            select = self.outputVar.get()

            # Run the recombineering module
            iterSeqs(fileN,beg,end,organism,cutsites,select)

            # Finished, exit the loop
            break

            # When finished, reset and wait for another round or exit
            # self.Reset()

    def Reset(self):
        self.eFile.delete(0,END)
        self.eFile.insert(0,'genes/')
        self.eBegin.delete(0,END)
        self.eBegin.insert(0,'40')
        self.eEnd.delete(0,END)
        self.eEnd.insert(0,'350')
        self.eCuts.delete(0,END)
        self.eCuts.insert(0,'GAATTC,GGATCC,AAGCTT')

##    def OpenOligos(self):
##        import subprocess,sys
##        if sys.platform == 'win32':
##            path = sys.path[0] + '..\Oligos'
##            subprocess.call('%SystemRoot%\explorer.exe '+path,shell=True)
##        else:
##            subprocess.call('open ../Oligos/',shell=True)


    def Help(self):
        import subprocess,sys
        
        if sys.platform == 'win32':

            self.Popup("Please refer to the User Manual for help")

            # Opening in Internet Explorer for some odd reason
##            path = sys.path[0] + '\..\PyRec_UserManual.pdf'
##            subprocess.call('%SystemRoot%\explorer.exe '+path,shell=True)
        else:
            subprocess.call('open ../PyRec_UserManual.pdf',shell=True)

    def Quit(self):
        self.F.quit()

app = PyRec()
app.F.mainloop()
