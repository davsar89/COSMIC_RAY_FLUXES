# Docker usage

Run these commands from the project root:

```bash
docker compose -f run_using_docker/docker-compose.yml build
docker compose -f run_using_docker/docker-compose.yml up --abort-on-container-exit
```

The image is built from the reviewed local source tree; it does not clone a mutable remote repository. It is named `cosmic-ray-fluxes:local`, runs as an unprivileged user, and contains only the compiled executable, its runtime library, and the required model input tables.

To run a different calculation without editing the Compose file:

```bash
docker compose -f run_using_docker/docker-compose.yml run --rm cosmic-ray-flux \
  33 0.3 15 20 -80 2019 5 27 0.15 standard
```

Arguments are documented in the project-level `README.md`. Pin `UBUNTU_VERSION` to an immutable image digest in controlled production builds when byte-for-byte reproducibility is required.
