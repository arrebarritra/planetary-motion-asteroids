import datetime
from nbody import *
from skyfield.api import load

dims = 3
# Create a timescale and ask the current time.
ts = load.timescale()
tinit = ts.from_datetime(datetime.datetime(
    2022, 1, 1, tzinfo=datetime.timezone.utc))
# tinit = ts.now()

# Load the JPL ephemeris DE421 (covers 1900-2050).
planets = load('de421.bsp')
sfsun, sfearth, sfmoon = planets['sun'], planets['earth'], planets['moon']

sun = Body(M=1.989e30, R=4.2635e-5, r0=sfsun.at(tinit).position.au,
           v0=sfsun.at(tinit).velocity.au_per_d, name="Sun", color='orange')
earth = Body(M=5.972e24, R=4.2635e-5, r0=sfearth.at(tinit).position.au,
             v0=sfearth.at(tinit).velocity.au_per_d, name="Earth", color='blue')
moon = Body(M=7.34767e22, R=1.1613e-5, r0=sfmoon.at(tinit).position.au,
            v0=sfmoon.at(tinit).velocity.au_per_d, name="Moon", color='grey')
sat = Body(M=1e5, R=1e-7, r0=earth.r0 + 40 * np.array([0, earth.R, 0]),
            v0=earth.v0, name="Sat", color='green')


dt = 0.01
solar = System([sun, earth, moon], dims, VerletIntegrator(dt))

tmax = 40
nsteps = int(tmax / dt + 1)

solar.run_animate(nsteps, stepsperframe=1, center=earth, includebodies=[earth, moon], trail=1000)

# while solar.t < tmax:
#     solar.update(printcollisions=True)
#
# solar.plot_trace(center=earth, includebodies=[earth, sat, moon])
# solar.plot_energy()
