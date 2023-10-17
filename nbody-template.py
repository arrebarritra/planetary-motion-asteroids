class Body:
    """
    Contains properties of gravitational body
    """

    def __init__(self, M, R, r0, v0):
        self.M = M
        self.R = R
        self.r = r0
        self.v = v0
        self.F = 0

    def add_force(self, F):
        """
        Add force F to body
        """


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


class System:
    """
    System containing n bodies
    """

    def __init__(self, bodies, dims, integrator):
        self.bodies = bodies
        self.n = len(bodies)
        self.dims = dims
        self.t = 0
        self.integrator = integrator

        self.KE = 0
        self.PE = 0
        # Set KE and PE
        self.E = self.KE + self.PE

    def update(self, iters=1):
        """Updates body positions and velocities, and energy of the system"""
        for i in range(self.n):
            for j in range(i + 1, self.n):
                # Calculate interactions between bodies i and j

    def potential(b1, b2):
        """Returns pairwise potential between bodies b1 and b2"""

    def force(b1, b2):
        """Returns force on body b1 from body b2"""
