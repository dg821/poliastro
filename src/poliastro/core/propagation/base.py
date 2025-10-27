from numba import njit as jit
import numpy as np


@jit
def func_twobody(t0, u_, k):
    """Differential equation for the initial value two body problem.

    Parameters
    ----------
    t0 : float
        Time.
    u_ : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Standard gravitational parameter.

    """
    x, y, z, vx, vy, vz = u_
    r3 = (x**2 + y**2 + z**2) ** 1.5

    du = np.array([vx, vy, vz, -k * x / r3, -k * y / r3, -k * z / r3])
    return du


@jit
def func_solar_sail(
    t0,
    u_,
    sail_normal_vec,
    k,
    sail_area_m2=15000,
    sc_mass_kg=500,
    star_flux_N_m2=5.4026e-6,
    ref_dist_km=149597870.691,
):
    def solar_sail_acceleration(
        state,
        sail_normal_vec,
        sail_area_m2,
        sc_mass_kg,
        star_flux_N_m2,
        ref_dist_km,
    ):
        """ """
        x, y, z, vx, vy, vz = state
        r_mag = np.linalg.norm(np.array([x, y, z]))

        u_r = -np.array([x, y, z]) / r_mag  # radial unit vector

        accel_sail = (
            -2 * star_flux_N_m2 * sail_area_m2 / sc_mass_kg * (ref_dist_km / r_mag) ** 2 * np.dot(sail_normal_vec, u_r) ** 2 * sail_normal_vec
        )

        return accel_sail

    du_kep = func_twobody(t0, u_, k)
    ax, ay, az = solar_sail_acceleration(
        u_,
        sail_normal_vec,
        sail_area_m2,
        sc_mass_kg,
        star_flux_N_m2,
        ref_dist_km,
    )
    du_ad = np.array([0, 0, 0, ax, ay, az])
    return du_kep + du_ad
