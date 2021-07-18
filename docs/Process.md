# Process

## Main equation

$$ t - t_p = \sqrt{\frac{a³}{\mu}} * \left[ corr + factor * sin^{-1}\left(\frac{\sqrt{1-e²} * sin(\nu)}{1 + e * cos(\nu)}\right) - e * \left( \frac{\sqrt{1-e²} * sin(\nu)}{1 + e * cos(\nu)}\right) \right] $$

with

- $a$ : semi-major axis.
- $\mu$ : Earth gravitation parameter.
- $e$ : eccentricity.
- $\nu$ : true anomaly.
- $corr$ : a correction value $[\pi]$.
- $factor$ : $-1$ or $+1$.

These values depend on the value of $\nu$, the mean anomaly.

## Step 1 : get the orbital elements

- semi-major axis $a$
- eccentricity $e$
- inclination $i$
- longitude of ascending node $\Omega$
- argument of periapsis $\omega$

## Step 2 : compute some other parameters

- radius at perigee and apogee $r_p$, $r_a$
- altitude at perigee and apogee $z_p$, $z_a$
- Velocity at perigee and apogee $v_p$, $v_a$
- orbital period $T$
- mean motion $\nu$

## Compute $\nu_c$

$$\nu_c = cos^{-1}(-e)$$

## find the corrections

This value of the critical anomaly is needed to determine the correction values $corr$ and $factor$ using the following logic :

...

if $-2\pi -\nu_c < \nu < -2\pi +\nu_c$, then $corr = -2\pi$ and $factor = +1$.

if $-2\pi +\nu_c < \nu < -\nu_c$, then $corr = -\pi$ and $factor = -1$.

if $-\nu_c < \nu < +\nu_c$, then there is no correction.

if $+\nu_c < \nu < 2\pi -\nu_c$, then $corr = +\pi$ and $factor = -1$.

if $2\pi -\nu_c < \nu < 2\pi +\nu_c$, then $corr = +2\pi$ and $factor = +1$.

...

These corrections are due to the presence of the arcsine function, which produces numerical errors in function of the true anomalies used.

## Find $t_p$ the time at perigee.

Suppose $t=0$ and $\nu = \omega$. Apply the main Equation with the appropriate correction to find $t_p$.

## Find $t$ the time at each true anomaly

Now that we know $t_p$, we can map the main formula to a list of desired true anomalies.

## Determine the Latitude and Longitude in Earth Fixed frame

$$ la = sin^{-1} \left(sin(i) * sin(\omega + \nu)\right)$$

$$ L_0 = sin^{-1} \left(\frac{tan(la)}{tan(i)}\right)$$

## Apply corrections to $L_0$

We are in the same case as before, the presence of arcine function induces numerical errors that we need to correct.
...

if $-\omega -360 < \nu < -\omega - 270$, then $corr = -360$ and $factor = +1$.

if $-\omega - 270 < \nu < -\omega - 90$, then $corr = -180$ and $factor = -1$.

if $-\omega - 90 < \nu < -\omega + 90$, then there is no correction.

if $-\omega + 90 < \nu < -\omega + 270$, then $corr = +180$ and $factor = -1$.

if $-\omega + 270< \nu < -\omega + 360$, then $corr = +360$ and $factor = +1$.

...

These corrections are inserted right before the arcsine in the expression of $L_0$.

## Determine the Latitude and Longitude in Earth Inertial frame

The latitude is still correct in the new inertial frame.

The longitude is defined by :
$$ L = \Omega + L_0(t) - \alpha * t$$

with $\alpha$ the rotation rate of the Earth : $\alpha = 360/86164$ °/s

---

And there you have it !
