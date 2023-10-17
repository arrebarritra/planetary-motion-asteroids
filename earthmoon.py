from planets import *

dims = 2
mtoau = 6.685e-12
mstoaud = 5.775e-7

earth = Body(M=5.972e24, R=4.2635e-5, r0=np.zeros(dims), v0=np.zeros(dims),
             name="Earth", color='b')
moon = Body(M=7.34767e22, R=1.1613e-5, r0=np.array(
    [3.844e8, 0])*mtoau, v0=np.array([0, 1.082e3])*mstoaud, name="Moon", color='grey')

dt = 0.01
eulercromer = EulerCromerIntegrator(dt)
euler = EulerIntegrator(dt)
verlet = VerletIntegrator(dt)

earthmoon = System([earth, moon], dims, verlet)

tmax = 20
nsteps = tmax // dt + 1

while earthmoon.t < tmax:
    earthmoon.update()

earthmoon.plot_trace(center=earth)
earthmoon.plot_energy()

# earthmoon.run_animate(nsteps, center=earth, trail=1000, printcollisions=True, stepsperframe=5)
