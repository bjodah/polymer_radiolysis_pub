#!/bin/bash
set -x
if [[ $IN_DOCKER == "1" ]]; then
    source /etc/bash.bashrc
    source $HOME/.bashrc
    set -e
    addgroup --gid "$HOST_GID" mygroup
    adduser --disabled-password --uid "$HOST_UID" --gid "$HOST_GID" --gecos '' myuser
    mkdir -m 0755 -p ~myuser/.config
    chown myuser:mygroup ~myuser/.config
    adduser myuser sudo
    echo "%sudo ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
    #su -l myuser -c "cd /work; $RUNCMD"
    mkdir /tmp/build; cd /tmp/build
    cp -r /work/external .
    make -C external
    cp /work/src/* .
    cp -r /opt/pykinetgine-0.3.5/external/* .
    cp -r /opt/pykinetgine-0.3.5/kinetgine .
    cp -r /opt/syme-*/{include,lib}/* .
    make
    cd /work/tests
    PATH=/tmp/build:$PATH make -f rad.mk rad-res.png
else
    if groups | grep docker >/dev/null; then
        DOCKERCMD=docker
    else
        DOCKERCMD="sudo docker"
    fi

    ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
    $DOCKERCMD run --rm --cap-add SYS_PTRACE -e IN_DOCKER=1 -e TERM -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) -v "$ABS_REPO_PATH":/work -w /work -it ${DOCKERIMAGE:-polymer_radiolysis_pub} bash --rcfile /work/scripts/$(basename $0)
fi
