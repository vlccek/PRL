#!/bin/bash

# kontrola na pocet argumentu
if [ $# -ne 1 ]; then
  exit 1;
fi;

# pocet procesu == 0, nema cenu pokracovat
if [ $1 -eq 0 ]; then
  exit 2;
fi;

# preklad
mpic++ --prefix /usr/local/share/OpenMPI -o oets oets.cpp

# vygenerovani nahodne posloupnosti cisel, pocet dan prvnim parametrem skriptu
dd if=/dev/random bs=1 count=$1 of=numbers 2>/dev/null

# spusteni aplikace (oversubscribe - vice procesu nez fyzicky k dispozici)
mpirun --oversubscribe --prefix /usr/local/share/OpenMPI -np $1 oets

# uklid
# rm -f oets numbers
