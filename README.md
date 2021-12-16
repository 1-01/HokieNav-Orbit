# HokieNav-Orbit

Matlab implementation of orbital propagation software for Virginia Tech's ECE senior design project with The Aerospace Corporation.

## What it is and mathematical overview

The motion (6D state vector `x`) of a spacecraft is simulated using the dynamical equations of the perturbed 2-body problem with respect to the ECI (Earth-centered inertial) frame. 
This means that the state `x` is defined

    x = (X, Y, Z, VX, VY, VZ) = (r, v)

(where `r = (X, Y, Z)` is the standard Cartesian position and `v = (VX, VY, VZ)` is the respective velocity)
and is determined in time by the dynamics

     dx
    ---- = f(t, x, u),
     dt

where `t` is time and `u` is the control. In the absence of control, the dynamics of position in Cartesian coordinates are simply related to their respective velocities

     dr
    ---- = v,
     dt
     
but the dynamics of the velocities are given by the accelerations

     dv        - GM
    ---- =  ------------ r + perturbations.
     dt      (r'*r)^1.5
     
Numerically integrating these equations using some Runge-Kutta scheme provides a simulation of the spacecraft's *trajectory* in time.
*Initial conditions* are needed to specify what the spacecraft's state `x` is at the initial time (of integration) `t0`.
These conditions are specified by supplying a *two-line element set* (TLE) as a `.tle` or `.txt` file.

### Perturbations

Perturbations are *small* accelerations that slowly change the orbital behavior and are caused by effects like
* **Aerodynamic drag**
* **Earth's oblateness** and nonuniformity
* Solar radiation pressure (SRP)
* Other-body gravitation (moon, sun, etc.)
* Relativity

The most prominant accelerations of this list are **Earth's oblateness** and **drag**. Earth's oblateness (and even its roughness) are relatively easy to model and be accurate; drag, however, is **not**.
Note that SRP, other-body effects, and relativity require an *ephemeris* (position/velocity of other celestial bodies).

## Some improvements to be made

With spacecraft in LEO orbiting at speeds exceeding *several* km/s, it's somewhat easy to lose position accuracy from modeling errors after even short amounts of time.

* **Slower-varying state**
  * Cartesian coordinates offer an easy way to write the equations of motion that will be singularity-free, but they oscillate quickly and hence aren't the most efficient set of coordinates
  * Some alternatives include: Encke's method, COE (classical/Keplerian orbital elements), EE (equinoctial elements), MEE (modified equinoctial elements)
    * Note that the COE are singular for both circular **and** equitorial orbits
* **Vectorize functions**
  * Doing this generally increases efficiency (at least in Matlab) and improves readability/organization
* **IAU 2000/2006 reduction**
  * Earth's rotation is more complicated than simply spinning about the J2000 +Z axis. See Matlab's `dcmeci2ecef.m` file help (the challenge here is making our own version that's compatible with Coder)
* **Improved SRP**
  * Earth's shadow is a complicated thing to model correctly because of our atmosphere (refraction of light). Getting *this* detailed is likely not necessary, but the current implementation could be generalized a little bit to account for higher orbits
* **TLEs (Use more information and write)**
  * More information than coordinates for initial conditions are given in the TLE that could be useful (e.g. B* for drag)
  * Write own TLE files for downlink
* **Improved events**
  * The time `t` and state `x` at the occurence of a specified event can be marked using the event function `odevents.m`
    * Such an improvement could be when the spacecraft establishes and loses line-of-sight of a location on Earth's surface (e.g. a groundstation)
* **Fully prepare for onboarding**
  * Features like plotting will not be needed when being run autonomously on station
* **Bug squashing**
  * No piece of software is written from scratch without mistakes, and some of them slip through. They're caught in time and through experience

## Things to watch out for

There can be many pitfalls when numerically integrating sets of ODEs - the models need to be correct **and** the integrator options need to be chosen properly. 
* It's possible to have a good model but bad integration options (e.g. tolerances).
  * This could make the solution appear to be bad and point blame towards the model, but the problem is really in the setup (options).
* Some models may not be worthwhile to implement. 
  * It's possible for too much time to be spent evaluating one model (which bottlenecks the program).
  * It's also possible that the model turns out to be *insignificant*. That is, its magnitude is much smaller than other models or the other models' uncertainties (as in, it just gets drowned out).
    * For example, if a perturbation is evaluating to magnitudes beneath machine precision, then the computer itself could be having more of an effect on the solution than the perturbation!

Be sure to take it slow when you need to, carefully evaluate what's happening, and keep a level head!

\- Matt
