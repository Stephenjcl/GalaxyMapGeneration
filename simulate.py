# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:36:07 2019

@author: Stephen
"""
import random
import matplotlib as m
import matplotlib.pyplot as plt
import sys
import os
import argparse
sys.path.append(os.getcwd())

from garnets import generate_stellar_system, random_star, plot_stellar_system, append_starcsv_function
from galaxy import Galaxy, plot_galaxy

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

m.use('TkAgg')
plt.rcParams['path.simplify']=False

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


'''
galaxy = Galaxy(35000) #Number of particles
galaxy.reset(13000, 4000, 0.0004, 0.85, 0.90, 0.5, 200, 300)
galaxy.update(100000)
#plot_galaxy(galaxy, dparg = 150)

exo = 0
#Loop to produce all planets for stars
for i in range(len(galaxy._stars)):
#    with HiddenPrints():
#        planets.append(generate_stellar_system(random_star()))
    with HiddenPrints():
        galaxy._stars[i]['planets'] = generate_stellar_system(random_star())
    #galaxy._stars[0]['planets'].planets[0] #Accessing a particulat planet
    exo += len(galaxy._stars[i]['planets'].planets)
    #print (exo, "exoplanets found surrounding star located at:", "%.2f" % galaxy._stars[i]['position'][0], "%.2f" % galaxy._stars[i]['position'][1])

    #Add galactic position to Name
    for n in len(galaxy._stars):
        for m in len(galaxy._stars[n]['planets'].planets):
            galaxy._stars[n]['planets'].planets[m].map_attributes.replace('######', ''.join(galaxy._stars[n]['position'][0], galaxy._stars[n]['position'][1]))

    print ('\r', "%.2f" % (i/len(galaxy._stars)*100), "% stars scanned,", exo, "total planets found.", end='', flush=True)

'''


if __name__ == '__main__':
    #m.use('TkAgg')
    plt.rcParams['path.simplify']=False
    # Construct the argument parser
    ap = argparse.ArgumentParser()

    # Add the arguments to the parser
    ap.add_argument("-t", "--test", required=False, dest='testflag', action='store_true', help="Flag to run a small sample of stars compared to the full program. Also generates a csv and sample solar system map.")

    ap.add_argument("-c", "--csv", required=False, dest='filesave', action='store_true', help="Save a csv of the data.")

    ap.add_argument("-p", "--plot", required=False, dest='plotsave', action='store_true', help="Save an image plot of every star system in the map. Takes up a LARGE amount of space!")

    options = ap.parse_args()
    # args = vars(ap.parse_args())

    #Create a generation error log
    errorlog = []
    errorlocation = []
    
    #Check for presence of galactic data file.
    try:
        with open('galactic_data.pkl', 'rb') as input:
            galaxy = pickle.load(input)
    except:
        # Break down input arguments
        if options.testflag is True:
            filename = "Test.csv"
            seed = 0.1234567890
            random.seed(seed)
            galaxy = Galaxy(1000) #Number of particles
            galaxy.reset(13000, 4000, 0.0004, 0.85, 0.90, 0.5, 200, 300)
            galaxy.update(100000)
        else:
            # random.seed('earth')
            filename = "Production Run.csv"
            seed = 0.1234567890
            random.seed(seed)
            if seed:
                random.seed(seed)
            galaxy = Galaxy(35000) #Number of particles
            galaxy.reset(13000, 4000, 0.0004, 0.85, 0.90, 0.5, 200, 300)
            galaxy.update(100000)
    
        exo = 0
        # Loop to produce all planets for stars
        for i in range(len(galaxy._stars)):
            try:
                with HiddenPrints():
                    galaxy._stars[i]['planets'] = generate_stellar_system(random_star())
                #galaxy._stars[0]['planets'].planets[0] #Accessing a particulat planet
                exo += len(galaxy._stars[i]['planets'].planets)
                #print (exo, "exoplanets found surrounding star located at:", "%.2f" % galaxy._stars[i]['position'][0], "%.2f" % galaxy._stars[i]['position'][1])
                
                #Star name:
                galaxy._stars[i]['planets'].name = (str("{:+06.0f}".format(galaxy._stars[i]['position'][0],1)) + "V" + str("{:+06.0f}".format(galaxy._stars[i]['position'][1],1)))    
                
                #Add galactic position to Name
                for z in range(len(galaxy._stars[i]['planets'].planets)):
                    galaxy._stars[i]['planets'].planets[z].mapname = galaxy._stars[i]['planets'].planets[z].mapname.replace('########', galaxy._stars[i]['planets'].name)
                    if len(galaxy._stars[i]['planets'].planets[z].moons) > 0: #Should be false if there is a value present
                        for n in range(0, len(galaxy._stars[i]['planets'].planets[z].moons)): #For each moon in list
                            galaxy._stars[i]['planets'].planets[z].moons[n].mapname = galaxy._stars[i]['planets'].planets[z].moons[n].mapname.replace('########', galaxy._stars[i]['planets'].planets[z].mapname)


    
                if options.filesave is True:
                    append_starcsv_function(galaxy._stars[i]['planets'], filename, i, len(galaxy._stars))
    
    
                print('\r', "%.2f" % (i/len(galaxy._stars)*100) + "% of total stars scanned,", exo, "total planets found.", end='', flush=True)
            except Exception as e:
                #print("Erronous planet generation procedure on star number: ", i, "please investigate further.")
                errorlocation.append(i-1)
                errorlog.append(str(e))
                pass


        print('\r', "100% of total stars scanned,", exo, "total planets found.", end='', flush=True)
        if options.filesave or options.textsave is True:
            print("Data saved.")
        print("Average of", round(exo/len(galaxy._stars),1), "planets per star.")

#        if len(errorlocation) > 0:
#            print("Sample error:", errorlog[0])
#            print("Status of atmosphere string:", galaxy._stars[errorlocation[0]]['planets'].planets[0].atmosphere)
#            print("Errors total:", len(errorlocation))
        save_object(galaxy, 'galactic_data.pkl')
        print("Saved pickle.")

    #plot_stellar_system
    if options.plotsave is True:
        if not os.path.isdir("System Plots"):
            os.makedirs("System Plots")
        print("Now producing stellar plots.")
        plt.switch_backend('Agg')
        plt.rcParams["image.cmap"] = "tab20"
        
        for i in range(len(galaxy._stars)):
            if len(galaxy._stars[i]['planets'].planets) > 0:
                fname = 'System Plots/' + galaxy._stars[i]['planets'].name + '.png'
                plot_stellar_system(galaxy._stars[i]['planets'], save = fname)
            print('\r', "%.2f" % (i/len(galaxy._stars)*100) + "% of total stars plotted.", end='', flush=True)
        print('\r', "Plotting complete.", end='', flush=True)

    #print('\r', "Data saved. Now plotting galactic image.", end ='', flush = True)
    #plot_galaxy(galaxy, save = "Test Image.png", dparg = 150, show = False)
    #plot_stellar_system(galaxy._stars[0]['planets'], save = 'Test Plot.png')
