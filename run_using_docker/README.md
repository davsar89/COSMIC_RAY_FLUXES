Using Docker, deploy the run rapidly on any computer (Linux, Mac, Windows).

Requires `docker` and `docker-compose` installed.

Commands to run :
* `docker compose build`
* `docker compose up`

Will grab an Ubuntu 22.04 image and install dependencies (`gfortran`, ...). Will then compile and run executable `electron_fluxes`.

Modify the command and the end of the file `docker-compose.yml` to change input parameters for the executable `electron_fluxes` (`particle ID (0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)`, `min energy threshold (MeV)`, `altitude (km)`, `latitude (deg)` and `longitude`).
