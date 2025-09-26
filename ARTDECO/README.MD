# Tkinter GUI in Docker with X11 Forwarding

This project enables running a Tkinter-based GUI application inside a Docker container while forwarding the graphical interface to the host system using X11.

## Getting Started

### 1. Build the Docker Image

Run the `docker_build.sh` script:

```bash
./docker_build.sh
```

Uses the included Dockerfile.

Detects the UID / GID of the user running the script.

Creates an identical “mirroring” user inside the image so file ownership and X-authority match between host and container.

Ensures the container can safely forward windows to the host’s X11 server.

### 2. Running the Docker Image
```bash
./docker_run.sh
```
The script:

Binds the X11 socket
Mounts /tmp/.X11-unix so the container can reach the X server.

Handles .Xauthority
Copies the host’s .Xauthority, adjusts permissions, and injects it into the container.

Keeps things isolated
Uses a dedicated network namespace and non-privileged user to limit container access while still allowing GUI output.

### 3. Launching the GUI

After docker_run.sh starts, the application’s Tkinter window appears on the host display. Interact with it as you would any local program.