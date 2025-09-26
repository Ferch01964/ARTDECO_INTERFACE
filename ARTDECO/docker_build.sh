#!/bin/bash
sudo usermod -aG docker $USER
docker build --build-arg UID=$(id -u $USER) --build-arg GID=$(id -g $USER) -t app1 .

mkdir -p "$HOME/app1_outputs"

export XSOCK=/tmp/.X11-unix
export XAUTH=/tmp/.docker.xauth
touch $XAUTH
xauth nlist "$DISPLAY" | sed -e 's/^..../ffff/' | xauth -f "$XAUTH" nmerge -

#docker run -it --volume=$XSOCK:$XSOCK:rw --volume=$XAUTH:$XAUTH:rw --env="XAUTHORITY=${XAUTH}" --env="DISPLAY" app1:latest
docker run -it --volume=$XSOCK:$XSOCK:rw --volume=$XAUTH:$XAUTH:rw --volume="$HOME/app1_outputs:/APP/INPUTS/INPUTS"  --env="XAUTHORITY=${XAUTH}" --env="DISPLAY" app1:latest
