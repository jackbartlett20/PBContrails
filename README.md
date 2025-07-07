**PBContrails** is a population balance model for contrails in a 0D perfectly stirred reactor (PSR). The code is adapted from [Stelios Rigopoulos' Chemistry and Particulates MODelling (CPMOD)](https://github.com/CambridgeUniversityPress/Population-Balance-of-Particles-in-Flows). It solves droplet growth, freezing, and ice crystal growth across a volume spectrum by calculating average hygroscopicity and dry volume fraction in each volume interval. It will hopefully also do aggregation if I get round to it.


# Compiling

The code is compiled with `make -f makefile` in the compile directory. This creates the file `cpmod`.


# Executing

Desired initial particle properties, environmental conditions, volume grid properties, and integration settings should be input by amending the files `pbe/pbe.in`, `pbe/species.in`, and `psr/psr.in`.

The program allows for a variable time step to avoid instabilities relating to large droplet growth rates (see below). Numerical errors can quickly cause instability in the program since the droplet growth rate is highly dependent on hygroscopicity and dry volume fraction which subsequently affects the evolution of those properties. If variable time step is enabled, the program will attempt to reduce the time step if the Courant number condition is broken or it is running into unphysical results. To help reduce instability, the Courant number condition can be chosen. Further, a tight Courant number condition is implemented for when the dry volume fraction of any interval exceeds 0.999. Ensure that $Courant_{tight} \leq Courant \leq 1$ if you want the program to do anything useful. I recommend that you choose these parameters by trial and error based on your scenario.

It is likely that the dry volume fraction in an interval will try to go above 1 or below 0 at some point. For this reason, a tolerance on this property can be chosen where $1 < f_{dry} < 1+tolerance$ results in $f_{dry} = 1$ and similarly for below 0. The tolerance should be as small as possible.

The program is executed with `./cpmod`.

Output is in the output directory.


# Notes on model dynamics

Plume dynamics currently follow [Kärcher et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD023491). It can be changed in the `pbe_set_environment` subroutine to suit any needs.

The droplet growth equation follows Pruppacher and Klett / Seinfeld and Pandis. Freezing and ice crystal growth is taken from [Kärcher et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD023491).

Particles species are initialised as lognormally-distributed in volume, but this can be edited without issue.

A uniform volume grid can be chosen, but I haven't tested if it still works and it's not very useful anyway given the enormous range of particle size present in most scenarios.

Runge-Kutta 2nd and 4th order integration is implemented (but needs updating). Only explicit Euler (Runge-Kutta 1st order) has been thoroughly tested.


# Acknowledgements

With thanks to the authors of CPMOD:
Stelios Rigopoulos
Fabian Sewerin
Anxiong Liu
Binxuan Sun
Daniel O’Sullivan