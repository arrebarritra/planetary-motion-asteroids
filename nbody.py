import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

G = 1.488e-34  # Gravitational constant G in days and AU
nameindex = 0


class Body:
    """
    Contains properties of gravitational body
    """

    def __init__(self, M, R, r0, v0, name='', color='g'):
        global nameindex

        self.M = M
        self.R = R
        self.r0 = r0
        self.r = np.copy(r0)
        self.rlist = np.array([r0])
        self.v0 = v0
        self.v = np.copy(v0)
        self.vlist = np.array([v0])
        self.F = 0
        self.KE = 0.5 * self.M * np.dot(self.v0, self.v0)

        if name == "":
            self.name = "Planet" + str(nameindex)
            nameindex += 1
        else:
            self.name = name

        self.color = color

    def add_force(self, F):
        """
        Add force F to body
        """
        self.F += F


class Integrator:
    """
    Parent class for integrators
    """

    def __init__(self, dt):
        self.dt = dt

    def integrate(self, body):
        """
        Integrate one time step, updating body position and velocity given the
        total force on it
        """
        pass


class EulerIntegrator(Integrator):
    """
    Applies the Euler integration scheme
    """

    def integrate(self, body):
        body.r += body.v * self.dt
        body.v += body.F / body.M * self.dt

        body.rlist = np.append(body.rlist, [body.r], axis=0)
        body.vlist = np.append(body.vlist, [body.v], axis=0)

        body.KE = 0.5 * body.M * np.dot(body.v, body.v)
        body.F = 0


class EulerCromerIntegrator(Integrator):
    """
    Applies the Euler-Cromer integration scheme
    """

    def integrate(self, body):
        body.v += body.F / body.M * self.dt
        body.r += body.v * self.dt

        body.vlist = np.append(body.vlist, [body.v], axis=0)
        body.rlist = np.append(body.rlist, [body.r], axis=0)

        body.KE = 0.5 * body.M * np.dot(body.v, body.v)
        body.F = 0


class VerletIntegrator(Integrator):
    """
    Applies the basic Verlet integration scheme
    """

    def integrate(self, body):
        if len(body.rlist) == 1:
            body.r += body.v0 * self.dt + 0.5 * body.F / body.M * self.dt**2
            body.v = (body.r - body.r0) / self.dt
        else:
            body.r = 2 * body.rlist[-1] - body.rlist[-2] +\
                     body.F / body.M * self.dt**2
            body.v = (body.r - body.rlist[-1]) / self.dt

        body.rlist = np.append(body.rlist, [body.r], axis=0)
        body.vlist = np.append(body.vlist, [body.v], axis=0)

        body.KE = 0.5 * body.M * np.dot(body.v, body.v)
        body.F = 0


class System:
    """
    System containing n bodies
    """

    def __init__(self, bodies, dims, integrator):
        self.bodies = bodies
        self.n = len(bodies)
        self.dims = dims
        self.t = 0
        self.ts = np.array([0])
        self.integrator = integrator
        self.collisions = []
        self.check_collisions()


        self.KE = 0
        for i in range(self.n):
            self.KE += self.bodies[i].KE
        self.KElist = np.array([self.KE])
        self.PE = 0
        for i in range(self.n):
            for j in range(i + 1, self.n):
                self.PE += System.potential(self.bodies[i], self.bodies[j])
        self.PElist = np.array([self.PE])
        self.E = self.KE + self.PE
        self.Elist = np.array([self.E])

    def update(self, iters=1):
        """Updates body positions and velocities, and energy of the system"""
        for i in range(iters):
            self.t += self.integrator.dt
            self.ts = np.append(self.ts, self.t)
            # Calculate interactions between all bodies
            self.KE = 0
            self.PE = 0
            for i in range(self.n):
                for j in range(i + 1, self.n):
                    F = System.force(self.bodies[i], self.bodies[j])
                    self.PE += System.potential(self.bodies[i],
                                                self.bodies[j])
                    self.bodies[i].add_force(F)
                    self.bodies[j].add_force(-F)

                self.integrator.integrate(self.bodies[i])
                self.KE += self.bodies[i].KE

            self.E = self.KE + self.PE
            self.KElist = np.append(self.KElist, self.KE)
            self.PElist = np.append(self.PElist, self.PE)
            self.Elist = np.append(self.Elist, self.E)
            self.check_collisions()

    def check_collisions(self):
        """Checks for collisions (intersection) between bodies"""
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if (np.linalg.norm(self.bodies[i].r - self.bodies[j].r)
                        < self.bodies[i].R + self.bodies[j].R):
                    self.collisions.append(Collision(self.t,
                                                     self.bodies[i],
                                                     self.bodies[j]))

    def potential(b1, b2):
        """Returns pairwise potential between bodies b1 and b2"""
        return -G*b1.M*b2.M / np.linalg.norm(b1.r - b2.r)

    def force(b1, b2):
        """Returns force on body b1 from body b2"""
        return -G*b1.M*b2.M / np.linalg.norm(b1.r - b2.r)**3 * (b1.r - b2.r)

    def orbital_period(self, c_body, o_body):
        """
        Returns orbital period of central body c_body around orbiting body
        o_body. Only works if the system has been propagated just longer than
        one orbital period.
        """
        rel_start_pos = o_body.r0 - c_body.r0

        min_dist = np.inf
        min_step = 0
        for i in range(len(self.ts)//2, len(self.ts)):
            rel_pos = o_body.rlist[i] - c_body.rlist[i]
            dist = np.linalg.norm(rel_pos - rel_start_pos)
            if dist < min_dist:
                min_dist = dist
                min_step = i

        return self.ts[min_step]

    def get_maxmin_coords(self, center, includeindices):
        maxmin = np.zeros((self.dims, 2))
        maxmin[:,0] = np.inf
        maxmin[:,1] = -np.inf

        for step in range(len(self.ts)):
            for i in includeindices:
                for j in range(self.dims):
                    if self.bodies[i].rlist[step, j] - center[step, j] < maxmin[j, 0]:
                        maxmin[j, 0] = self.bodies[i].rlist[step, j] - center[step, j]
                    elif self.bodies[i].rlist[step, j] - center[step, j] > maxmin[j, 1]:
                        maxmin[j, 1] = self.bodies[i].rlist[step, j] - center[step, j]

        return maxmin

    def center_coords(system, center):
        if center == 0:
            center = np.zeros(system.dims)
        elif isinstance(center, np.ndarray):
            if center.shape == (system.dims,):
                pass
        elif isinstance(center, str):
            for i in range(system.n):
                if system.bodies[i].name == center:
                    center = system.bodies[i].r
        elif isinstance(center, Body):
            center = center.r

        return center

    def center_coords_list(system, center):
        if center == 0:
            center = np.zeros((int(system.t / system.integrator.dt + 1),
                               system.dims))
        elif isinstance(center, np.ndarray):
            if center.shape == (system.dims,):
                center = np.tile(center,
                                 (int(system.t / system.integrator.dt + 1), 1))
        elif isinstance(center, str):
            for i in range(system.n):
                if system.bodies[i].name == center:
                    center = system.bodies[i].rlist
        elif isinstance(center, Body):
            center = center.rlist

        return center

    def body_indices(system, include):
        indices = []
        if include is None:
            return np.arange(system.n)
        elif isinstance(include, np.ndarray) or isinstance(include, list):
            if isinstance(include[0], int):
                return include
            elif isinstance(include[0], str):
                indices = []
                for name in include:
                    for i in range(system.n):
                        if system.bodies[i].name == name:
                            indices.append(i)

                return indices
            elif isinstance(include[0], Body):
                indices = []
                for body in include:
                    for i in range(system.n):
                        if system.bodies[i].name == body.name:
                            indices.append(i)

                return indices

    def plot_trace(self, center=0, includebodies=None):
        center = System.center_coords_list(self, center)

        fig = plt.figure()
        if self.dims == 2:
            ax = fig.add_subplot()
            ax.axis('equal')
        elif self.dims == 3:
            ax = fig.gca(projection='3d')

        includeindices = System.body_indices(self, includebodies)

        for i in includeindices:
            if self.dims == 2:
                ax.plot(self.bodies[i].rlist[:, 0] - center[:, 0],
                        self.bodies[i].rlist[:, 1] - center[:, 1],
                        color=self.bodies[i].color, label=self.bodies[i].name)
                ax.text(0.02, 0.95, 't='+str(np.round(self.t, 2))+' days',
                        transform=ax.transAxes)
            elif self.dims == 3:
                ax.plot(self.bodies[i].rlist[:, 0] - center[:, 0],
                        self.bodies[i].rlist[:, 1] - center[:, 1],
                        self.bodies[i].rlist[:, 2] - center[:, 2],
                        'o-', markevery=[-1],
                        color=self.bodies[i].color, label=self.bodies[i].name)
                ax.text(0.02, 0.95, 0, 't='+str(np.round(self.t, 2))+' days',
                        transform=ax.transAxes)

        plt.legend()
        plt.show()

    def compare_trace(self, compare, center=0, includebodies=None):
        fig = plt.figure()
        if self.dims == 2:
            ax = fig.add_subplot()
            ax.axis('equal')
            ax.text(0.02, 0.95, 't='+str(np.round(self.t, 2))+' days',
                    transform=ax.transAxes)
        elif self.dims == 3:
            ax = fig.gca(projection='3d')
            ax.text(0.02, 0.95, 0, 't='+str(np.round(self.t, 2))+' days',
                    transform=ax.transAxes)

        scenter = System.center_coords_list(self, center)
        ccenter = System.center_coords_list(compare, center)
        includeindices = System.body_indices(self, includebodies)

        for i in includeindices:
            if self.dims == 2:
                ax.plot(self.bodies[i].rlist[:, 0] - scenter[:, 0],
                        self.bodies[i].rlist[:, 1] - scenter[:, 1],
                        color=self.bodies[i].color, label=self.bodies[i].name)
                ax.plot(compare.bodies[i].rlist[:, 0] - ccenter[:, 0],
                        compare.bodies[i].rlist[:, 1] - ccenter[:, 1],
                        color=compare.bodies[i].color, label=compare.bodies[i].name)
            elif self.dims == 3:
                ax.plot(self.bodies[i].rlist[:, 0] - scenter[:, 0],
                        self.bodies[i].rlist[:, 1] - scenter[:, 1],
                        self.bodies[i].rlist[:, 2] - scenter[:, 2],
                        'o-', markevery=[-1],
                        color=self.bodies[i].color, label=self.bodies[i].name)
                ax.plot(compare.bodies[i].rlist[:, 0] - ccenter[:, 0],
                        compare.bodies[i].rlist[:, 1] - ccenter[:, 1],
                        compare.bodies[i].rlist[:, 2] - ccenter[:, 2],
                        'x:', markevery=[-1],
                        color=compare.bodies[i].color, label=compare.bodies[i].name)

        plt.tight_layout()
        plt.legend()
        plt.show()

    def plot_energy(self):
        fig, axs = plt.subplots(3, 1, sharex='col')

        axs[0].plot(self.ts, self.KElist, label="KE")
        axs[0].legend()

        axs[1].plot(self.ts, self.PElist, label="PE")
        axs[1].legend()

        axs[2].plot(self.ts[2:], self.Elist[2:], label="Total energy")
        axs[2].legend()
        plt.xlabel('t')
        plt.show()

    def run_animate(self, nsteps=1000, stepsperframe=1, trail=10, center=0,
                    includebodies=None, savedir=''):
        frames = nsteps // stepsperframe + 1
        self.update(nsteps)
        plt.clf()
        plt.rcParams["figure.figsize"] = (10, 10)

        centerlist = System.center_coords_list(self, center)
        includeindices = System.body_indices(self, includebodies)

        max_coords = self.get_maxmin_coords(centerlist, includeindices)

        if self.dims == 2:
            ax = plt.subplot()
            ax.set_xlim(max_coords[0, 0], max_coords[0, 1])
            ax.set_ylim(max_coords[1, 0],max_coords[1, 1])
            ax.axis('equal')
            lines = ax.plot(np.empty((0, len(includeindices))),
                            np.empty((0, len(includeindices))),
                            'o-', markevery=[-1])
            text = ax.text(0, 0, 't='+str(np.round(self.t, 2))+' days',
                           transform=ax.transAxes)
        elif self.dims == 3:
            ax = plt.subplot(projection='3d')
            ax.set_xlim(max_coords[0, 0], max_coords[0, 1])
            ax.set_ylim(max_coords[1, 0],max_coords[1, 1])
            ax.set_zlim(max_coords[2, 0],max_coords[2, 1])
            lines = [ax.plot([self.bodies[i].rlist[0,0]], [self.bodies[i].rlist[0,1]],
                             [self.bodies[i].rlist[0,2]], 'o-', markevery=[-1],
                             color=self.bodies[i].color,
                             label=self.bodies[i].name)[0]
                     for i in includeindices]
            text = ax.text(0, 0, 0,
                           't='+str(np.round(self.t, 2))+' days',
                           transform=ax.transAxes)

        plt.legend()
        plt.tight_layout()

        anim = animation.FuncAnimation(plt.gcf(), animate,
                                       fargs=[self, lines, text, self.dims,
                                              center, includebodies, trail,
                                              stepsperframe],
                                       frames=frames, interval=1,
                                       blit=True, repeat=False)
        if savedir:
            anim.save(savedir)
        else:
            plt.show()


class Collision:
    """Contains information on collisions between bodies"""

    def __init__(self, t, b1, b2, printmessage=False):
        self.t = t
        self.name1 = b1.name
        self.name2 = b2.name
        self.r1 = b1.rlist[-1]
        self.r2 = b2.rlist[-1]
        self.v1 = b1.vlist[-1]
        self.v2 = b2.vlist[-1]
        self.message = self.name1 + " collided with " + \
            self.name2 + " at t=" + str(self.t)
        print(self)

    def __str__(self):
        return self.message


def animate(framenr, system, lines, text, dims, center, includebodies,
            trail, stepsperframe):
    steps = framenr * stepsperframe
    center = System.center_coords_list(system, center)

    if trail > steps:
        trail = steps

    includeindices = System.body_indices(system, includebodies)
    text.set_text('t='+str(np.round(system.ts[steps], 2))+' days')

    j = 0
    for i in includeindices:
        lines[j].set_data(system.bodies[i].rlist[steps-trail:steps, 0]
                          - center[steps-trail:steps:, 0],
                          system.bodies[i].rlist[steps-trail:steps, 1]
                          - center[steps-trail:steps, 1])
        if dims == 3:
            lines[j].set_3d_properties(
                system.bodies[i].rlist[steps-trail:steps, 2]
                                       - center[steps-trail:steps, 2])

        j += 1

    return np.append(lines, text)
