#!/usr/bin/env python
from scipy.optimize import fsolve
from scipy.optimize import minimize
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

    def rel_correction(self, E):
        return -self.a**2 * math.pi / self.Porb * (E + e * math.sin(E))

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

    def get_rel_phi(self, tobs):
        """Gets the phase model for a particular tobs, including relativistic correction"""
        # First, generate E from tobs
        E = fsolve(self.numeric_E, math.pi, args=(tobs))
        # Then, solve for t using the pulsar time equation
        t = self.Porb / (2 * math.pi) * (E - self.e * math.sin(E)) + self.rel_correction(E)
        # Finally use t/P + phi0 to get the phase
        return t / self.P + self.phi0


################################################################################


# Prepare all the parameters for testing the models
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


# First, a unit test for the Newtonian model
def test_newtonian(pulsar):
    """Compare the outputs of the pulsar and test pulsar classes"""
    # Create a pulsar test object
    tpulsar = Pulsar_Test(a=pulsar.a, b=pulsar.b, Porb=pulsar.Porb, nx=pulsar.nx, ny=pulsar.ny, P=pulsar.P, t0=pulsar.t0, phi0=pulsar.phi0)
    # Generate a table of tobs and phi values
    tobs, phi = tpulsar.generate(0.5, 4, 0.1)
    # Try generating phi values using the test pulse's tobs
    gen_phi = np.zeros(tobs.shape)
    for i, t in enumerate(tobs):
        gen_phi[i] = pulsar.get_phi(t)
    # Print everything for comparison
    for t, p, gen_p in zip(tobs, phi, gen_phi):
        print t, abs(p - gen_p)


pulsar = Pulsar(e, Porb, a, b, nx, ny, P, t0, phi0)
test_newtonian(pulsar)


# A degeneracy test for the relativistic model (DOES NOT YET WORK)
def newtonian(tobs, e, Porb, a, b, nx, ny, P, t0, phi0):
    """I have defined a method on its own here to replicate the work of the
    newtonian class, since we will not have values to initialize a pulsar object with."""
    # The lambda function is equivalent to the numeric_E function in the class
    g = lambda E: Porb / (2 * math.pi) * (E - e * math.sin(E)) - a * nx * (math.cos(E) - e) - b * ny * math.sin(E) + t0 - tobs
    # Then from E, we return phi (this next bit is the equivalent of the get_phi method)
    E = fsolve(g, math.pi)
    # Then, solve for t using the pulsar time equation
    t = Porb / (2 * math.pi) * (E - e * math.sin(E))
    # Finally use t/P + phi0 to get the phase
    return t / P + phi0


def delta(args, ni=100, nd=10):
    e, Porb, a, b, nx, ny, P, t0, phi0 = args
    """Defines the function to be minimized for the test_relativistic method:
    the square of the difference between the non-relativistic phi and the
    relativistic phi, calculated over ni subintervals"""
    delta = 0
    tobs = np.linspace(t0, t0 + nd, ni)
    # Here, the 'pulsar' object definied for the newtonian test  will represent our
    # relativistic pulsar (we are optimizing the parameters for the newtonian model --
    # the relativistic pulsar phase model will be calculated as normal so we can still use the class)
    for i in xrange(ni):
        delta += (pulsar.get_phi(tobs[i]) - newtonian(tobs[i], e, Porb, a, b, nx, ny, P, t0, phi0))**2
        # delta += (pulsar.get_rel_phi(tobs[i]) - newtonian(tobs[i], e, Porb, a, b, nx, ny, P, t0, phi0))**2
    print delta
    return delta


def test_rel_degeneracy():
    """Degeneracy test for relativistic model"""
    sol = minimize(delta, [e, Porb, a, b, nx, ny, P, t0, phi0])
    print sol


test_rel_degeneracy()
