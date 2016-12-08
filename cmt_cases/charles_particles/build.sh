#!/bin/bash

if [ ! -f makenek ]; then
    cp ../../nek/makenek makenek
    sed -i '6s/.*/SOURCE_ROOT="$(pwd)\/..\/..\/nek" /' makenek
fi

./makenek box
