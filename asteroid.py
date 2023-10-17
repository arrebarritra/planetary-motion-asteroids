import datetime
from nbody import *
from skyfield.api import load

ts = load.timescale()
tinit = ts.tt(2022, 1, 1, 12, 0)

planets = load('de421.bsp')
sfsun, sfearth, sfmoon = planets['sun'], planets['earth'], planets['moon']

sun = Body(M=1.989e30, R=4.6547e-3, r0=sfsun.at(tinit).position.au,
           v0=sfsun.at(tinit).velocity.au_per_d, name="Sun", color='orange')
earth = Body(M=5.972e24, R=4.2635e-5, r0=sfearth.at(tinit).position.au,
             v0=sfearth.at(tinit).velocity.au_per_d, name="Earth", color='blue')
moon = Body(M=7.348e22, R=1.1613e-5, r0=sfmoon.at(tinit).position.au,
            v0=sfmoon.at(tinit).velocity.au_per_d, name="Moon", color='grey')
astsun = Body(M=1.989e30, R=4.6547e-3, r0=sfsun.at(tinit).position.au,
           v0=sfsun.at(tinit).velocity.au_per_d, name="Sun", color='orange')
astearth = Body(M=5.972e24, R=4.2635e-5, r0=sfearth.at(tinit).position.au,
             v0=sfearth.at(tinit).velocity.au_per_d, name="Earth", color='blue')
astmoon = Body(M=7.348e22, R=1.1613e-5, r0=sfmoon.at(tinit).position.au,
            v0=sfmoon.at(tinit).velocity.au_per_d, name="Moon", color='grey')

astr0 = 1.01 * earth.r0
astv0 = (sfearth.at(tinit + 5).position.au - astr0) / 4.85
ast = Body(M=1e24, R=1.7560e-6, r0=astr0, v0=astv0, name="Asteroid", color='black')


dims = 3
dt = 0.01
solar = System([sun, earth, moon], dims, VerletIntegrator(dt))
astsys = System([astsun, astearth, astmoon, ast], dims, VerletIntegrator(dt))

tmax = 60
nsteps = int(tmax / dt)

# astsys.run_animate(nsteps, stepsperframe=1, center=astearth, includebodies=[astearth, astmoon, ast], trail=1000)
solar.update(nsteps)
astsys.update(nsteps)

astsys.compare_trace(solar, center="Earth", includebodies=["Earth", "Moon"])
# print("Moon pos", np.linalg.norm((moon.r - earth.r) - (astmoon.r - astearth.r)))
# astsys.compare_trace(solar, center="Sun", includebodies=["Sun", "Earth"])
# print("Earth pos", np.linalg.norm((earth.r - sun.r) - (astearth.r - astsun.r)))
# astsys.plot_energy()
