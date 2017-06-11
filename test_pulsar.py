#!/usr/bin/env python
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import numpy as np


class Pulsar_Test():
    """Test class"""
    def __init__(self, a, b, Porb, nx, ny, P, t0, phi0):
        self.a = a
        self.b = b
        self.e = math.sqrt(1 - self.b**2 / self.a**2)
        self.Porb = Porb
        self.nx = nx
        self.ny = ny
        self.P = P
        self.t0 = t0
        self.phi0 = phi0

    def x(self, t, E):
        return self.a * (math.cos(E) - self.e)

    def y(self, t, E):
        return self.b * math.sin(E)

    def E_to_tobs(self, E):
        t = self.Porb / (2 * math.pi) * (E - self.e * math.sin(E))
        return t - self.nx * self.x(t, E) - self.ny * self.y(t, E) + self.t0

    def E_to_phi(self, E):
        t = self.Porb / (2 * math.pi) * (E - self.e * math.sin(E))
        return t / P + phi0

    def generate(self, start, end, step):
        E_table = np.arange(start, end, step)
        tobs = np.zeros(E_table.shape)
        phi = np.zeros(E_table.shape)
        for i, E in enumerate(E_table):
            tobs[i] = self.E_to_tobs(E)
            phi[i] = self.E_to_phi(E)
        return tobs, phi


class Pulsar():
    """Actual class"""
    def __init__(self, e, Porb, a, b, nx, ny, P, t0, phi0):
        self.e = e
        self.Porb = Porb
        self.a = a
        self.b = b
        self.nx = nx
        self.ny = ny
        self.anx = a * nx
        self.bny = b * ny
        self.P = P
        self.t0 = t0
        self.phi0 = phi0

    def numeric_E(self, E, tobs):
        """Solve for E numerically"""
        t_pulsar = self.Porb / (2 * math.pi) * (E - self.e * math.sin(E))
        return t_pulsar - self.anx * (math.cos(E) - self.e) - self.bny * math.sin(E) + self.t0 - tobs

    def plot(self, tobs):
        """Just for visualization -- plot E against the function that defines it"""
        E_sample = np.arange(0, 2 * math.pi, 0.01)
        f = np.zeros(E_sample.shape)
        for i, E in enumerate(E_sample):
            f[i] = self.f(E, tobs)
        plt.plot(E_sample, f)
        plt.xlabel('E')
        plt.ylabel('function we want to find the root of')
        plt.show()

    def get_phi(self, tobs):
        """Gets the phase model for a particular tobs"""
        # First, generate E from tobs
        E = fsolve(self.numeric_E, math.pi, args=(tobs))
        # Then, solve for t using the pulsar time equation
        t = self.Porb / (2 * math.pi) * (E - self.e * math.sin(E))
        # Finally use t/P + phi0 to get the phase
        return t / self.P + self.phi0

    def test_newtonian(self):
        """Compare the outputs of the pulsar and test pulsar classes"""
        # Create a pulsar test object
        tpulsar = Pulsar_Test(a=self.a, b=self.b, Porb=self.Porb, nx=self.nx, ny=self.ny, P=self.P, t0=self.t0, phi0=self.phi0)
        # Generate a table of tobs and phi values
        tobs, phi = tpulsar.generate(0.5, 4, 0.1)
        # Try generating phi values using the test pulse's tobs
        gen_phi = np.zeros(tobs.shape)
        for i, t in enumerate(tobs):
            gen_phi[i] = self.get_phi(t)
        # Print everything for comparison
        for t, p, gen_p in zip(tobs, phi, gen_phi):
            print t, abs(p - gen_p)

# Prepare all the parameters for both models
# I assume it is not good that the units here aren't consistent, but
# I guess it works for the purposes of this test :)
e = 0.070560
Porb = 2.35769683    # orbital period (days)
a = 2.152813         # semimajor axis (lt s)
b = math.sqrt(a**2 * (1 - e**2))
nx = math.sqrt(2)
ny = math.sqrt(2)
P = 311.49341784442  # pulse frequency (Hz)
t0 = 51602.18629     # epoch of periastron (MJD)
phi0 = 0
anx = a * nx
bny = b * ny

pulsar = Pulsar(e, Porb, a, b, nx, ny, P, t0, phi0)
pulsar.test_newtonian()
