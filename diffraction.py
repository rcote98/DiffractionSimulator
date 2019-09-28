#!/usr/bin/python3

import matplotlib.pyplot as plt
from PIL import Image
import numpy as np

N = 1000 # number of points

global zoom
zoom = 1/4 #zoom in the diffraction graphs

# DEFAULT VALUES
global SIDE
SIDE = 1 # target side, in cm
global RADIUS
RADIUS = 0.5 # circle radius, in mm
global EDGE
EDGE = 1 # square edge, in mm
global LDA
LDA = 600e-9 # wavelength, in m
global Z
Z = 100 # obstacle-sensor distance, in cm
global FNAME
FNAME = "img/atom.png"

#############################################
## GENERACIÓN DE SUPERFICIES DE DIFRACCIÓN ##

def circle(radius):

    """Generates an NxN array of 0s with 1s within the 
    given radius from the center (in cm)"""

    bitmap = np.zeros((N,N))

    for i in range(N):
        for j in range(N):

            deltax, deltay = index_to_distance(i, j) # in cm

            if(np.sqrt(deltax**2 + deltay**2) < radius):
                bitmap[i,j] = 1

    return bitmap

def rectangle(side):

    """Generates an NxN array of 0s with 1s within a centered 
    square of the indicated side (in cm)"""

    bitmap = np.zeros((N,N))

    for i in range(N):
        for j in range(N):

            x, y = index_to_distance(i,j)

            if(np.abs(x) < side/2 and np.abs(y) < side/2):
                bitmap[i,j] = 1

    return bitmap


def index_to_distance(i, j):

    "Returns mapped distance on the array, in cm"
    
    ratio = float(SIDE)/float(N) 
    x = (i - N/2)*ratio
    y = (N/2 - j)*ratio
    
    return [x, y]

#############################################
### DIFRACTION FUNCTIONS ####################


def fraunhoffer(u1):

    """Calculates the intensity in the sensor at z using 
    Fraunhoffer's far-field approximation"""

    u1 = u1.astype(complex)

    fu1 = np.fft.fft2(u1)/(N*N)
    fu1 = np.fft.fftshift(fu1)

    P = np.abs(fu1)

    I = (1/(LDA*(Z/100))**2)*P # z to m

    return I

def fresnel(u1):

    """Calculates the intensity in the sensor at k using 
    Fresnel's approximation"""

    k = 2*np.pi/LDA
    u = u1.astype(complex)

    for i in range(N):
        for j in range(N):
            x, y = index_to_distance(i,j) # in cm
            u[i,j] = u[i,j]*np.exp((1j*k)/(2*Z/100)*((x/100)**2 + (y/100)**2))
    
    fu1 = np.fft.fft2(u)/(N*N)
    fu1 = np.fft.fftshift(fu1)

    P = np.abs(fu1)

    I = (1/(LDA*(Z/100))**2)*P # z to m

    return I
    


#############################################
############### MENU STUFF ##################

def menu():

    print("")
    print("### DIFFRACTION SIMULATOR ### ")
    print("")
    print("Raúl Coterillo")
    print("")
    print("This program calculates the diffraction pattern")
    print("of several shapes using either Fresnel's or")
    print("Fraunhoffer's approximation. For more details,")
    print("refer to the source code.")
    print("")
    print("First, select the shape of the obstacle: (default: Circle)")
    print("[0] - Circle")
    print("[1] - Square") 
    print("[2] - Image")

    choice = input("-> ")
    print("")

    #defaulting to circle
    if(choice == ""):
        choice = "0"

    if choice == "0":

        print("Please, now input the following information:\n")

        global SIDE
        rply = input_def("Obstacle/sensor area side: (default: "+str(SIDE)+"cm)", SIDE)
        SIDE = float(rply)

        global RADIUS
        rply = input_def("Circle radius in mm: (default: "+str(RADIUS)+"mm)", RADIUS)
        RADIUS = float(rply)

        global Z
        rply = input_def("Obstacle-sensor separaton in cm: (default: "+str(Z)+"cm)", Z)
        Z = float(rply)

        global LDA
        rply = input_def("Light wavelength in m: (default: "+str(LDA)+"m)", LDA)
        LDA = float(rply)

        print("Finally, choose the diffraction method: (default: Fresnel's)")
        print("[0] - Fresnel")
        print("[1] - Fraunhoffer")
        method = input("-> ")
        print("")

        if(method == "1"):
            print("Using Fraunhoffer's approximation...")
            obs = circle(RADIUS/10) # converted to cm
            diff = fraunhoffer(obs)
        else:
            print("Using Fresnel's approximation...")
            obs = circle(RADIUS/10) # converted to cm
            diff = fresnel(obs)

        print("DONE!")
        plot_stuff(obs, diff)
        quit()



    if choice == "1":

        print("Please, now input the following information:\n")

        rply = input_def("Obstacle/sensor area side: (default: "+str(SIDE)+"cm)", SIDE)
        SIDE = int(rply)

        global EDGE
        rply = input_def("Square side in mm: (default: "+str(EDGE)+"mm)", EDGE)
        EDGE = int(rply)

        rply = input_def("Obstacle-sensor separaton in cm: (default: "+str(Z)+"cm)", Z)
        Z = float(rply)

        rply = input_def("Light wavelength in m: (default: "+str(LDA)+"m)", LDA)
        LDA = float(rply)

        print("Finally, choose the diffraction method: (default: Fresnel's)")
        print("[0] - Fresnel")
        print("[1] - Fraunhoffer")
        method = input("-> ")
        print("")

        if(method == "1"):
            print("Using Fraunhoffer's approximation...")
            obs = rectangle(EDGE/10) # converted to cm
            diff = fraunhoffer(obs)
        else:
            print("Using Fresnel's approximation...")
            obs = rectangle(EDGE/10) # converted to cm
            diff = fresnel(obs)

        print("DONE!")
        plot_stuff(obs, diff)
        quit()

    if choice == "2":

        print("Please, now input the following information:\n")

        global FNAME
        rply = input_def("Image file name: (default: "+ FNAME +", image must be " 
            + str(N) + "x" + str(N) + " pixels)", FNAME)
        FNAME = rply

        try:
            obs = np.array(Image.open(rply).convert("L"))
        except FileNotFoundError:
            print("\n File does not exist, quitting.")
            quit()
        except OSError:
            print("\n File type not recognized, quitting.")
            quit()

        rply = input_def("Obstacle/sensor area side: (default: "+str(SIDE)+"cm)", SIDE)
        SIDE = int(rply)

        rply = input_def("Obstacle-sensor separaton in cm: (default: "+str(Z)+"cm)", Z)
        Z = float(rply)

        rply = input_def("Light wavelength in m: (default: "+str(LDA)+"m)", LDA)
        LDA = float(rply)

        print("Finally, choose the diffraction method: (default: Fresnel's)")
        print("[0] - Fresnel")
        print("[1] - Fraunhoffer")
        method = input("-> ")
        print("")

        if(method == "1"):
            print("Using Fraunhoffer's approximation...")
            diff = fraunhoffer(obs)
        else:
            print("Using Fresnel's approximation...")
            diff = fresnel(obs)

        print("DONE!")
        plot_stuff(obs, diff)
        quit()

    else:

        print("")
        print("Please, select a valid pattern using the numbers provided.")
        print("")
        quit()


def input_def(ph, def_val):

    print(ph)
    rply = input("-> ")
    print("")

    if(rply != ""):
        return rply
    else:
        return def_val

def plot_stuff(obstacle, difracted):

    plt.figure(figsize=(900, 400))
    factor = 1/(LDA*1000*Z*10)

    plt.suptitle(r"$z$ = " + str(int(Z)) + "cm\n\n" + r"$\lambda$ = " + str(LDA) + "m", fontsize = 20)

    plt.subplot(1,3,1)
    plt.title("Obstacle", fontsize=16)
    plt.imshow(obstacle, cmap="gray", extent=[-SIDE*5,SIDE*5,-SIDE*5, SIDE*5])
    plt.xlabel("y/mm", fontsize=14)
    plt.ylabel("x/mm", fontsize=14)

    plt.subplot(1,3,2)
    plt.title("Diffraction Pattern", fontsize=16)
    plt.imshow(difracted[int(N/2-N*zoom/2):int(N/2+N*zoom/2)-1,int(N/2-N*zoom/2):int(N/2+N*zoom/2)-1],
    cmap="gray", extent=[factor*(-SIDE*5*zoom),factor*(SIDE*5*zoom),factor*(-SIDE*5*zoom),factor*(SIDE*5*zoom)])
    plt.xlabel("y/mm", fontsize=14)
    plt.ylabel("x/mm", fontsize=14)

    plt.subplot(1,3,3)
    plt.title("Diffraction Pattern (log)", fontsize=16)
    plt.imshow(np.log(difracted[int(N/2-N*zoom/2):int(N/2+N*zoom/2)-1,int(N/2-N*zoom/2):int(N/2+N*zoom/2)-1]),
    cmap="gray", extent=[factor*(-SIDE*5*zoom),factor*(SIDE*5*zoom),factor*(-SIDE*5*zoom),factor*(SIDE*5*zoom)])
    plt.xlabel("y/mm", fontsize=14)
    plt.ylabel("x/mm", fontsize=14)

    plt.show()

## MAIN METHOD

if __name__ == "__main__":

    """ Launch this on execute."""

    menu()  

    print("")
