import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from tkinter import ttk

def calculate():
    L = float(L_entry.get())
    r1 = float(r1_entry.get())
    r2 = float(r2_entry.get())
    n = int(n_entry.get())
    m1 = float(m1_entry.get())
    Cp1 = int(Cp1_entry.get())
    rho1 = int(rho1_entry.get())
    m2 = float(m2_entry.get())
    Cp2 = int(Cp2_entry.get())
    rho2 = int(rho2_entry.get())
    T1i = int(T1i_entry.get())
    T2i = int(T2i_entry.get())
    T0 = int(T0_entry.get())
    U = int(U_entry.get())

    # Perform calculations or any other operations with the user inputs
    pi = 3.14159
    Ac1 = pi*r1**2
    Ac2 = pi*(r2**2 - r1**2)
    dx = L/n
    t_final = 1000
    dt = 1
    x = np.linspace(dx/2, L-dx/2, n)
    T1 = np.ones(n)*T0
    T2 = np.ones(n)*T0
    dT1dt = np.zeros(n)
    dT2dt = np.zeros(n)
    t = np.arange(0, t_final, dt)
    for j in range(1, len(t)):
        dT1dt[1:n] = (m1*Cp1*(T1[0:n-1]-T1[1:n])+U*2*pi*r1*dx*(T2[1:n]-T1[1:n]))/(rho1*Cp1*dx*Ac1)
        dT1dt[0] = (m1*Cp1*(T1i-T1[0])+U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho1*Cp1*dx*Ac1)
        dT2dt[1:n] = (m2*Cp2*(T2[0:n-1]-T2[1:n])-U*2*pi*r1*dx*(T2[1:n]-T1[1:n]))/(rho2*Cp2*dx*Ac2)
        dT2dt[0] = (m2*Cp2*(T2i-T2[0])-U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho2*Cp2*dx*Ac2)
        T1 = T1 + dT1dt*dt
        T2 = T2 + dT2dt*dt
        plt.figure(1)
        plt.plot(x, T1, color = "blue", label = "Inside")
        plt.plot(x, T2, color = "red", label = "Outside")
        plt.axis([0, L, 298, 820])
        plt.xlabel("Distance (m)")
        plt.ylabel("Temperature (K)")
        plt.legend(loc = "upper right")
        plt.title("Temperature Distribution in a Pipe")
        plt.pause(0.005)
        plt.clf()
    plt.show()

# Create the main tkinter window
root = tk.Tk()
root.title("Input Form")

# Create and place labels and entry widgets

tk.Label(root, text="Pipe Length (m)").grid(row=0, column=0)
L_entry = tk.Entry(root)
L_entry.grid(row=1, column=0)

tk.Label(root, text="Inner Pipe Radius (m)").grid(row=0, column=1)
r1_entry = tk.Entry(root)
r1_entry.grid(row=1, column=1)

tk.Label(root, text="Outer Pipe Radius (m)").grid(row=0, column=2)
r2_entry = tk.Entry(root)
r2_entry.grid(row=1, column=2)

tk.Label(root, text="Node Number").grid(row=2, column=1)
n_entry = tk.Entry(root)
n_entry.grid(row=3, column=1)

tk.Label(root, text="Fluid1").grid(row=4, column=1)

tk.Label(root, text="Mass Flow Rate (kg/s)").grid(row=5, column=0)
m1_entry = tk.Entry(root)
m1_entry.grid(row=6, column=0)

tk.Label(root, text="Heat Capacity of Fluid").grid(row=5, column=1)
Cp1_entry = tk.Entry(root)
Cp1_entry.grid(row=6, column=1)

tk.Label(root, text="Density of Fluid").grid(row=5, column=2)
rho1_entry = tk.Entry(root)
rho1_entry.grid(row=6, column=2)

tk.Label(root, text="Fluid2").grid(row=7, column=1)

tk.Label(root, text="Mass Flow Rate (kg/s)").grid(row=8, column=0)
m2_entry = tk.Entry(root)
m2_entry.grid(row=9, column=0)

tk.Label(root, text="Heat Capacity of Fluid").grid(row=8, column=1)
Cp2_entry = tk.Entry(root)
Cp2_entry.grid(row=9, column=1)

tk.Label(root, text="Density of Fluid").grid(row=8, column=2)
rho2_entry = tk.Entry(root)
rho2_entry.grid(row=9, column=2)

tk.Label(root, text="Inlet Temperature (K)").grid(row=10, column=0)
T1i_entry = tk.Entry(root)
T1i_entry.grid(row=11, column=0)

tk.Label(root, text="Pipe Inner Surface Temperature (K)").grid(row=10, column=1)
T2i_entry = tk.Entry(root)
T2i_entry.grid(row=11, column=1)

tk.Label(root, text="Innitial Temperature (K)").grid(row=10, column=2)
T0_entry = tk.Entry(root)
T0_entry.grid(row=11, column=2)

tk.Label(root, text="Overall Heat Transfer Coefficient (W/m^2)").grid(row=12, column=1)
U_entry = tk.Entry(root)
U_entry.grid(row=13, column=1)


# Create a button to trigger the calculation
calculate_button = tk.Button(root, text="Calculate", command=calculate)
calculate_button.grid(row=14, column=1)

# Start the tkinter main loop
root.mainloop()
