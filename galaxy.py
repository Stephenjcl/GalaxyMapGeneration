# -----------------------------------------------------------------------------
#  A Galaxy Simulator based on the density wave theory
#  (c) 2012 Ingo Berg
#
#  Simulating a Galaxy with the density wave theory
#  http://beltoforion.de/galaxy/galaxy_en.html
#
#  Python version(c) 2014 Nicolas P.Rougier
# -----------------------------------------------------------------------------
import math
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
#from os import chdir
#chdir('I:\Shared drives\Personal\Stephen\Personal\Creative Writing\World Building Materials\Galaxy Map Generation\Accrete\Galaxy')

class Galaxy(object):
    """ Galaxy simulation using the density wave theory """

    def __init__(self, n=30000):
        """ Initialize galaxy """

        # Excentricity of the innermost ellipse
        self._inner_excentricity = 0.8

        # Excentricity of the outermost ellipse
        self._outer_excentricity = 1.0

        #  Velovity at the innermost core in km/s
        self._center_velocity = 30

        # Velocity at the core edge in km/s
        self._inner_velocity = 200

        # Velocity at the edge of the disk in km/s
        self._outer_velocity = 300

        # Angular offset per parsec
        self._angular_offset = 0.019

        # Inner core radius
        self._core_radius = 6000

        # Galaxy radius
        self._galaxy_radius = 15000

        # The radius after which all density waves must have circular shape
        self._distant_radius = 0

        # Distribution of stars
        self._star_distribution = 0.45

        # Angular velocity of the density waves
        self._angular_velocity = 0.000001

        # Number of stars
        self._stars_count = n

        # Number of dust particles
        self._dust_count = int(self._stars_count * 0.75)

        # Number of H-II regions
        self._h2_count = 200

        # Particles
        dtype = [ ('theta',       float, 1),
                  ('velocity',    float, 1),
                  ('angle',       float, 1),
                  ('m_a',         float, 1),
                  ('m_b',         float, 1),
                  ('size',        float, 1),
                  ('type',        float, 1),
                  ('temperature', float, 1),
                  ('brightness',  float, 1),
                  ('position',    float, 2),
                  ('planets',     list, 1) ]
        n = self._stars_count + self._dust_count + 2*self._h2_count
        self._particles = np.zeros(n, dtype=dtype)

        i0 = 0
        i1 = i0  + self._stars_count
        self._stars = self._particles[i0:i1]
        self._stars['size'] = 4
        self._stars['type'] = 0

        i0 = i1
        i1 = i0  + self._dust_count
        self._dust = self._particles[i0:i1]
        self._dust['size'] = 64
        self._dust['type'] = 1

        i0 = i1
        i1 = i0 + self._h2_count
        self._h2a = self._particles[i0:i1]
        self._h2a['size'] = 64
        self._h2a['type'] = 2

        i0 = i1
        i1 = i0 + self._h2_count
        self._h2b = self._particles[i0:i1]
        self._h2b['size'] = 8
        self._h2b['type'] = 3


    def __len__(self):
        """ Number of particles """

        if self._particles is not None:
            return len(self._particles)
        return 0


    def __getitem__(self, key):
        """ x.__getitem__(y) <==> x[y] """

        if self._particles is not None:
            return self._particles[key]
        return None


    def reset(self, rad, radCore, deltaAng,
                    ex1, ex2, sigma, velInner, velOuter):

        # Initialize parameters
        # ---------------------
        self._inner_excentricity = ex1
        self._outer_excentricity = ex2
        self._inner_velocity = velInner
        self._outer_velocity = velOuter
        self._angular_offset = deltaAng
        self._core_radius = radCore
        self._galaxy_radius = rad
        self._distant_radius = self._galaxy_radius * 2
        self.m_sigma = sigma

        # Initialize stars
        # ----------------
        stars = self._stars
        R = np.random.normal(0, sigma, len(stars)) * self._galaxy_radius
        stars['m_a']        = R
        stars['angle']      = 90 - R * self._angular_offset
        stars['theta']      = np.random.uniform(0, 360, len(stars))
        stars['temperature']= np.random.uniform(3000, 9000, len(stars))
        stars['brightness'] = np.random.uniform(0.1, 0.5, len(stars))
        stars['velocity']   = 0.000005
        for i in range(len(stars)):
            stars['m_b'][i] = R[i]* self.excentricity(R[i])

        # Initialize dust
        # ---------------
        dust = self._dust
        X = np.random.uniform(0, 2*self._galaxy_radius, len(dust))
        Y = np.random.uniform(-self._galaxy_radius, self._galaxy_radius, len(dust))
        R = np.sqrt(X*X+Y*Y)
        dust['m_a']         = R
        dust['angle']       = R * self._angular_offset
        dust['theta']       = np.random.uniform(0, 360, len(dust))
        dust['velocity']    = 0.000005
        dust['temperature'] = 6000 + R/4
        dust['brightness']  = np.random.uniform(0.015,0.025)
        for i in range(len(dust)):
            dust['m_b'][i] = R[i] * self.excentricity(R[i])

        # Initialise H-II
        # ---------------
        h2a, h2b = self._h2a, self._h2b
        X = np.random.uniform(-self._galaxy_radius, self._galaxy_radius, len(h2a))
        Y = np.random.uniform(-self._galaxy_radius, self._galaxy_radius, len(h2a))
        R = np.sqrt(X*X+Y*Y)

        h2a['m_a']        = R
        h2b['m_a']        = R + 1000

        h2a['angle']      = R * self._angular_offset
        h2b['angle']      = h2a['angle']

        h2a['theta']      = np.random.uniform(0, 360, len(h2a))
        h2b['theta']      = h2a['theta']

        h2a['velocity']   = 0.000005
        h2b['velocity']   = 0.000005

        h2a['temperature'] = np.random.uniform(3000,9000,len(h2a))
        h2b['temperature'] = h2a['temperature']

        h2a['brightness']  = np.random.uniform(0.010,0.015, len(h2a))
        h2b['brightness']  = h2a['brightness']

        for i in range(len(h2a)):
            h2a['m_b'][i] = R[i] * self.excentricity(R[i])
        h2b['m_b'] = h2a['m_b']


    def update(self, timestep=100000):
        """ Update simulation """

        self._particles['theta'] += self._particles['velocity'] * timestep

        P = self._particles
        a,b = P['m_a'], P['m_b']
        theta, beta = P['theta'], -P['angle']

        alpha = theta * math.pi / 180.0
        cos_alpha = np.cos(alpha)
        sin_alpha = np.sin(alpha)
        cos_beta  = np.cos(beta)
        sin_beta  = np.sin(beta)
        P['position'][:,0] = a*cos_alpha*cos_beta - b*sin_alpha*sin_beta
        P['position'][:,1] = a*cos_alpha*sin_beta + b*sin_alpha*cos_beta

        D = np.sqrt(((self._h2a['position'] - self._h2b['position'])**2).sum(axis=1))
        S = np.maximum(1,((1000-D)/10) - 50)
        self._h2a['size'] = S
        self._h2b['size'] = np.maximum(S/6,1)


    def excentricity(self, r):

        # Core region of the galaxy. Innermost part is round
        # excentricity increasing linear to the border of the core.
        if  r < self._core_radius:
            return 1 + (r / self._core_radius) * (self._inner_excentricity-1)

        elif r > self._core_radius and r <= self._galaxy_radius:
            a = self._galaxy_radius - self._core_radius
            b = self._outer_excentricity - self._inner_excentricity
            return self._inner_excentricity + (r - self._core_radius) / a * b

        # Excentricity is slowly reduced to 1.
        elif r > self._galaxy_radius and r < self._distant_radius:
            a = self._distant_radius - self._galaxy_radius
            b = 1 - self._outer_excentricity
            return self._outer_excentricity + (r - self._galaxy_radius) / a * b

        else:
            return 1


# for i in range(len(galaxy._particles)):
#     if(galaxy._particles[i]['position'][1]>0):
#         print(galaxy._particles[i]['position'][1])


#How to calculate colours through temperature
def colour_temp(temp):
	# Algorithm for color temp taken from http://www.tannerhelland.com/4435/convert-temperature-rgb-algorithm-code/
	temp = temp / 100
	if temp <= 66:
		r = 255
	else:
		r = temp - 60
		r = 329.698727446 * (r ** -0.1332047592)
		r = min(255, max(0, r))

	if temp < 66:
		g = temp
		g = 99.4708025861 * math.log(g) - 161.1195681661
		g = min(255, max(0, g))
	else:
		g = temp - 60
		g = 288.1221695283 * (g ** -0.0755148492)
		g = min(255, max(0, g))

	if temp >= 65:
		b = 255
	elif temp < 20:
		b = 0
	else:
		b = temp - 10
		b = 138.5177312231 * math.log(b) - 305.0447927307
		b = min(255, max(0, b))

	return (r/255, g/255, b/255)

#Create a simple plot

#m.use('TkAgg')
plt.rcParams['path.simplify']=False

def plot_galaxy(galaxy, save = None, dparg = 900, show = True):
    fig, ax = plt.subplots() #Create plotting area
    x=[] #Prep Empty Variables
    y=[]
    r=[] #Size
    t=[] #Colour
    #a=[] #alpha (brightness?)

    alphabase = galaxy['brightness'].tolist()
    alphabase = [g*2 for g in alphabase]
    #snorm = np.

    for i in range(len(galaxy['position'])):
        x.append(galaxy['position'][i][0])
        y.append(galaxy['position'][i][1])
        # a.append((alphabase[i]))
        #This bit below needed to also store alpha values
        temp = None
        temp = np.array(colour_temp(galaxy['temperature'][i]))
        temp = temp.reshape([1,3])
        temp = np.c_[temp, alphabase[i]]
        t.append(temp)

        r.append(galaxy['size'][i]*0.5)

    #Correct my poor use of numpy
    t = np.vstack(t)
    #size = (galaxy._galaxy_radius*7.69230769e-7)*.035 #Deprecated size calculation
    ax.scatter(x,y,c=t, s = r, marker=".", linewidths = 0)
    plt.xlim(np.percentile(x, 5), np.percentile(x, 95))
    plt.ylim(np.percentile(y, 5), np.percentile(y, 95))
    plt.yticks([])
    plt.xticks([])

    ax.set_facecolor('black')
    #plt.axis('off')

    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.subplots_adjust(bottom = 0)
    fig.subplots_adjust(top = 1)
    fig.subplots_adjust(right = 1)
    fig.subplots_adjust(left = 0)

    if save is not None:
        fig.savefig(save, dpi=dparg, bbox_inches=extent, pad_inches=0)

    if show is True:
        fig.show()
    #plt.close()

####Lets animate the galaxy forming, shall we?####
#Pieces for progress bar

# =============================================================================
galaxy = Galaxy(70000) #Number of particles
galaxy.reset(13000, 4000, 0.0004, 0.90, 0.90, 0.5, 200, 300)
if not os.path.isdir("Frames"):
    os.makedirs("Frames")
    
Frames = 700 #1800 was old number
prog = 0
for i in range(Frames):
    num = str(i)
    num = num.zfill(4)
    fname = "Frames/GD_" + num + ".png"
    #galaxy.update(7222222) #Cover 13 billionyears/60 seconds at 30fps means this many years per frame.
    galaxy.update(250000)
    plot_galaxy(galaxy, save = fname, dparg = 150)
    plt.close()

    #Progress check
    prog += 1
    if prog == int(0.01*Frames):
        print("\r",((i/(Frames-1))*100), "% Completion...")
        prog = 0
print("Frame rendering completed... Compiling to video")

#chdir('I:\Shared drives\Personal\Stephen\Personal\Creative Writing\World Building Materials\Galaxy Map Generation\Accrete\Galaxy\Frames')
os.path.join(os.path, "Frames")
filenames = os.listdir()
with imageio.get_writer('GalaxySimulation.gif', mode='I', fps=30) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
# =============================================================================

#m.use('MacOSX')
#galaxy2 = Galaxy(70000) #Number of particles
#galaxy2.reset(13000, 4000, 0.0004, 0.85, 0.90, 0.5, 200, 300)
#galaxy2.update(100000)
#plot_galaxy(galaxy2, dparg = 150)
