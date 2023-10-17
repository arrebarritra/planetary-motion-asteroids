import datetime
from nbody import *
from skyfield.api import load

def max_error(system, body, sfname):
    max_error = 0
    for i in range(nsteps + 1):
        t = tinit + system.ts[i]
        error = np.linalg.norm(planets[sfname].at(t).position.au
                               - body.rlist[i])
        if error > max_error:
            max_error = error

    return max_error

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

dims = 3
dt = 0.001
solar = System([sun, earth, moon], dims, EulerIntegrator(dt))

tmax = 30 # Set to >366 for Earth around Sun, >28 for Moon around Earth
nsteps = int(tmax / dt)

solar.update(nsteps)
print(solar.orbital_period(earth, moon))
print(max_error(solar, moon, 'moon'))
# print(np.linalg.norm(earth.r0 - moon.r0))
solar.plot_energy()
solar.plot_trace(center=earth, includebodies=[earth, moon])
