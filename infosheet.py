# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:36:03 2019

@author: Stephen
"""
import random
import matplotlib as m
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
import sys
import os
import argparse
from fpdf import FPDF
from math import sqrt
import re

sys.path.append(os.getcwd())

                                                   
from garnets import generate_stellar_system, random_star, append_starcsv_function
from galaxy import Galaxy, plot_galaxy

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
    

def plot_info_sheet(stellar_system, save = None, show = False):
    w, h = figaspect(1/4)
    fig, ax = plt.subplots(figsize=(w,h)) #Create plotting area
    x=[] #Prep Empty Variables
    y=[]
    r=[]
    t=[]
    cdict = {0: 'Unknown', 1: 'Rock', 2: 'Venusian', 3: 'Terrestrial',
             4: 'Sub-Sub-Gas Giant', 5: 'Sub-Gas Giant', 6: 'Gas Giant',
             7: 'Martian', 8: 'Water', 9: 'Ice', 10: 'Asteroids', 11: 'Unknown'}
    plt.title(stellar_system.name)
    for i in range(0, len(stellar_system.planets)):
        x.append(stellar_system.planets[i].orbit.a) #Orbital distance in AU
        y.append(0)
        r.append((sqrt(stellar_system.planets[i].radius/149597870.7)**2)*400000) #Planet size in AU
        t.append(stellar_system.planets[i].type.value)

    #Add Moons
    for i in range(0, len(stellar_system.planets)):
        if len(stellar_system.planets[i].moons) > 0: #Should be false if there is a value present
            for n in range(0, len(stellar_system.planets[i].moons)): #For each moon in list
                x.append(stellar_system.planets[i].orbit.a) #Set the x to be parent bodies orbit
                # y.append(abs(log10(stellar_system.planets[i].moons[n].orbit.a*1500)) + 0.15) #Set y to the moons orbit from parent a.*1000 changed to a*250
                # y.append((n/len(stellar_system.planets[i].moons))*0.2 + 0.15) # Always plot moons equally between 0.15 and 0.85
                y.append((n+1)*0.05) # Always plot moons equally between 0.15 and 0.85
                r.append((sqrt(stellar_system.planets[i].moons[n].radius/149597870.7)**2)*400000)
                t.append(stellar_system.planets[i].moons[n].type.value)
    
    plt.xscale("log")

    plt.yscale("linear")
    plt.ylim(-.1,max(y)+0.5)
    plt.xlim(0.02,max(x)+(max(r)))

    scatter = ax.scatter(x,y,c=t, s=r)

    ax.axvline(x=1, c = "green", alpha = 0.3) #Earth Oribtal Distance
    plt.ylabel("Relative Planet size (not to scale)")
    plt.yticks([])
    plt.xlabel("Orbital Distance (log$_{10}$(AU))")

    legend1 = ax.legend(*scatter.legend_elements(),
                        loc="lower left", title="Planet Type", prop={'size' : 6})#,
                        #labels = cdict)

    for fob in range(0,len(legend1.texts)):
        ix = [int(s) for s in re.findall(r'\b\d+\b', str(legend1.texts[fob]))]
        ix = ix[2:]
        ix = int(ix[0])
        legend1.texts[fob].set_text(cdict[ix])
    ax.add_artist(legend1)

    if save is not None:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    plt.close()
        
def make_infosheet(planet):
    #Produce the pdf
    pdf = FPDF('P', 'in', 'Letter')
#    pdf.add_font('Kontrapunkt-Light', '', 'Kontrapunkt-Light.ttf', uni=True)
#    pdf.add_font('Kontrapunkt-Bold', '', 'Kontrapunkt-Bold.ttf', uni=True)


    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    
    text = planet + " System Information"
    pdf.cell(8.5, .5, text, 0, 1, 'C')
    pdf.image(fname, y = 1, w = 7.71259843)
    pdf.ln(3)
    
    #Add plantary output (basic)
    if moonlocation == 0:
        text = galaxy._stars[location]['planets'].planets[planetlocation]
    else:
        text = galaxy._stars[location]['planets'].planets[planetlocation].moons[moonlocation]
    text = repr(text)
    text = text.split('\n')
    text = '\n'.join(text)
    
    pdf.set_font('Arial', '', 8)
    pdf.multi_cell(3.856299215, .15, text, border = 'R', align = "L")
    
    #3.856299215
    #Fun details

    text = planet + '.pdf'
    pdf.output(text, 'F')
    
    
#Demo planet ID: +06260V+00274-XXG0XXXX-D520685B
if __name__ == '__main__':
    plt.rcParams['path.simplify']=False
    # Construct the argument parser
    ap = argparse.ArgumentParser(description = "This application generates a system plot, a planetary image, and description. These are then saved into a convenient pdf.")

    # Add the arguments to the parser
    ap.add_argument("-p", "--planet", required=True, dest='pname', help="Exact string name of the system you would like to investigate.")

    options = ap.parse_args()
    sysname = options.pname[:13]
    
    try:
        with open('galactic_data.pkl', 'rb') as input:
            galaxy = pickle.load(input)
    except:
        print("No data found. Please run 'simulate.py' to produce a set of galactic data.")
        raise

    if not os.path.isdir("System Plots"):
        os.makedirs("System Plots")
    print("Now producing system plot.")
    plt.switch_backend('Agg')
    plt.rcParams["image.cmap"] = "tab20"

    location = 0
    for i in range(len(galaxy._stars)):
        if hasattr(galaxy._stars[i]['planets'], 'name'): 
            if galaxy._stars[i]['planets'].name == sysname:
                fname = 'System Plots/' + galaxy._stars[i]['planets'].name + '_infosheet.png'
                plot_info_sheet(galaxy._stars[i]['planets'], save = fname)
                print('\r', "%.2f" % (i/len(galaxy._stars)*100) + "% completed scanning.", end='', flush=True)
                location = i
                
                if len(options.pname) > 25:
                    for j in range(len(galaxy._stars[i]['planets'].planets)):
                        for n in range(len(galaxy._stars[i]['planets'].planets[j].moons)):
                            if (galaxy._stars[i]['planets'].planets[j].moons[n].mapname) == options.pname:
                                moonlocation = n
                                planetlocation = j
                else:
                    if galaxy._stars[i]['planets'].planets[j].mapname == sysname:
                        planetlocation = j
                        moonlocation = 0
                        
    print('\r', "Plotting complete. File saved in System Plots", end='', flush=True)
    
    make_infosheet(options.pname)

    
    