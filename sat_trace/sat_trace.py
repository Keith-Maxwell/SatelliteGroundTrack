import matplotlib.pyplot as plt
import numpy as np

EARTH_RADIUS = 6378  # km
EARTH_GRAV_CONST = 398600  # km³/s²
ALPHA = 360 / 86164


class SatelliteTrace:
    def __init__(
        self,
        semi_major_axis: float,
        eccentricity: float,
        inclination: float,
        longitude_ascending_node: float,
        argument_of_periapsis: float,
    ):
        self.sma = semi_major_axis  # km
        self.ecc = eccentricity
        self.inc = inclination  # deg
        self.LAN = longitude_ascending_node  # deg
        self.AOP = argument_of_periapsis  # deg
        self._standard_parameters()
        self._initialize()
        _ = self.compute_latitude_longitude(30, 0, 360)

    @property
    def sma(self):
        return self._sma

    @sma.setter
    def sma(self, x):
        if x < EARTH_RADIUS:
            raise ValueError(
                "Semi-major axis cannot be lower than the Earth radius"
            )
        elif x > 1_000_000:
            raise ValueError(
                "Semi-major axis cannot be higher than the Earth "
                "influence sphere"
            )
        self._sma = x

    @property
    def ecc(self):
        return self._ecc

    @ecc.setter
    def ecc(self, e):
        if e >= 1 or e < 0:
            raise ValueError("Eccentricity must be between 0 and 1.")
        self._ecc = max(e, 0.00001)  # to avoid 0

    @property
    def inc(self):
        return self._inc

    @inc.setter
    def inc(self, i):
        if not -180 <= i <= 180:
            raise ValueError("Inclination must be between -180° and 180°")
        self._inc = i

    @property
    def LAN(self):
        return self._LAN

    @LAN.setter
    def LAN(self, la):
        if not 0 <= la <= 360:
            raise ValueError(
                "Longitude of Ascending Node must be between 0° and 360"
            )
        self._LAN = la

    @property
    def AOP(self):
        return self._AOP

    @AOP.setter
    def AOP(self, a):
        if not 0 <= a <= 360:
            raise ValueError(
                "Argument of periapsis must be between 0° and 360"
            )
        self._AOP = a

    def _standard_parameters(self) -> None:
        """This method computes the standard parameters from the orbital
        parameters. The standard parameters include the :
        - radius
        - altitude
        - period
        - velocity
        """
        self.radius_ap = self.sma * (1 + self.ecc)  # km
        self.radius_pe = self.sma * (1 - self.ecc)  # km
        self.altitude_ap = self.radius_ap + EARTH_RADIUS  # km
        self.altitude_pe = self.radius_pe + EARTH_RADIUS  # km
        self.period = (
            2 * np.pi * np.sqrt(np.power(self.sma, 3) / EARTH_GRAV_CONST)
        )  # seconds
        self.velocity_ap = np.sqrt(
            2
            * (
                -(EARTH_GRAV_CONST / (2 * self.sma))
                + (EARTH_GRAV_CONST / self.radius_ap)
            )
        )  # km/s
        self.velocity_pe = np.sqrt(
            2
            * (
                -(EARTH_GRAV_CONST / (2 * self.sma))
                + (EARTH_GRAV_CONST / self.radius_pe)
            )
        )  # km/s

    def _initialize(self) -> None:
        """Compute the critical true anomaly and the time at perigee."""
        # initialize
        self.nu_crit = round(np.degrees(np.arccos(-self.ecc)))
        nu_init = -self.AOP
        self.t_p = self._compute_t(nu_init)

    def compute_latitude_longitude(
        self, step: int = 30, start: int = 0, stop: int = 360
    ) -> tuple:
        """This function computes the latitude and the longitude of the
        satellite in the Earth inertial frame. It updates the values
        stored in the instance as well as returns them.

        Parameters
        ----------
        step : int, optional
            true anomaly step between two points, in degrees, by default 30
        start : int, optional
            start true anomaly in degrees, by default 0.
        stop : int, optional
            stop true anomaly in degrees, by default 360.

        Returns
        -------
        tuple[np.array[float], np.array[float]]
            An array of latitudes and an array of longitudes.
        """
        # compute t in function of nu
        self.nu_list = np.arange(start, stop, step)
        self.t_list = np.array(
            [self._compute_t(nu) - self.t_p for nu in self.nu_list]
        )
        # compute the list of latitudes
        self.latitude = np.degrees(
            np.arcsin(
                np.sin(
                    np.radians(self.inc)
                    * np.sin(np.radians(self.AOP + self.nu_list))
                )
            )
        )
        # compute the list of raw longitudes
        raw_longitude = np.array(
            [
                self._compute_longitude(nu, lat)
                for nu, lat in zip(self.nu_list, self.latitude)
            ]
        )
        # compute the true longitudes, with rotating Earth
        self.longitude = np.array(
            [
                self.LAN + raw_longitude[i] - ALPHA * t
                for i, t in enumerate(self.t_list)
            ]
        )
        return self.latitude, self.longitude

    def _determine_correction_for_t(self, v_c, v):
        """Due to an arcsine function, we must apply corrections"""
        # Thanks to Louis ETIENNE for this code
        if -v_c <= v <= v_c:
            correction = 0
            factor = 1
        elif v < -v_c:
            k = 0
            while v <= -(2 * np.pi) * ((k + 1) // 2) + (-1) ** (k + 1) * v_c:
                k += 1
                correction = -k * np.pi
                factor = (-1) ** k
        else:
            k = 0
            while v >= (2 * np.pi) * ((k + 1) // 2) + (-1) ** k * v_c:
                k += 1
                correction = k * np.pi
                factor = (-1) ** k
        return correction, factor

    def _determine_correction_for_lo(self, aop, v):
        """Due to an arcsine function, we must apply corrections"""
        # Thanks to Louis ETIENNE for this code
        if -aop - np.pi / 2 <= v <= -aop + np.pi / 2:
            correction = 0
            factor = 1
        elif v < -aop - np.pi / 2:
            k = 0
            while v <= -aop - np.pi * (k + 0.5):
                k += 1
                correction = -k * np.pi
                factor = (-1) ** k
        else:
            k = 0
            while v >= -aop + np.pi * (k + 0.5):
                k += 1
                correction = k * np.pi
                factor = (-1) ** k

        if np.radians(self.inc) > np.pi / 2:
            correction *= -1
        return correction, factor

    def _compute_t(self, nu: float) -> float:
        """Function to compute the time associated to each true anomaly"""
        correction, factor = self._determine_correction_for_t(
            np.radians(self.nu_crit), np.radians(nu)
        )
        return np.sqrt(self.sma ** 3 / EARTH_GRAV_CONST) * (
            correction
            + factor
            * np.arcsin(
                np.sqrt((1 - np.power(self.ecc, 2)))
                / (1 + self.ecc * np.cos(np.deg2rad(nu)))
                * np.sin(np.deg2rad(nu))
            )
            - self.ecc
            * (
                np.sqrt((1 - np.power(self.ecc, 2)))
                / (1 + self.ecc * np.cos(np.deg2rad(nu)))
                * np.sin(np.deg2rad(nu))
            )
        )

    def _compute_longitude(self, nu: float, lat: float) -> float:
        """Get the longitude in function of latitude and true anomaly"""
        # input everything in degrees
        correction, factor = self._determine_correction_for_lo(
            np.radians(self.AOP), np.radians(nu)
        )
        return np.degrees(
            correction
            + factor
            * np.arcsin(np.tan(np.radians(lat)) / np.tan(np.radians(self.inc)))
        )

    def get_latitude(self) -> list[float]:
        return self.latitude

    def get_longitude(self) -> list[float]:
        return self.longitude


def plot_trace(lat: list[float], lon: list[float]) -> None:
    # TODO: handle -180 < longitudes < 180
    # TODO: return figure object ?
    plt.plot(lon, lat, marker="o")
    plt.xlabel("Longitudes")
    plt.ylabel("Latitudes")
    plt.show()


def plot_traces(*args: SatelliteTrace):
    if len(args) == 0:
        raise ValueError(
            "Please provide at least 1 instance of SatelliteTrace"
        )
    for sat_trace in args:
        plt.plot(sat_trace.longitude, sat_trace.latitude, marker="o")
    plt.xlabel("Longitudes")
    plt.ylabel("Latitudes")
    plt.show()


if __name__ == "__main__":
    sat1 = SatelliteTrace(40708, 0.8320, 61, 120, 270)
    sat2 = SatelliteTrace(20000, 0.120, 26, 80, 270)
    lat, lon = sat1.get_latitude(), sat1.get_longitude()

    plot_traces(sat1, sat2)
