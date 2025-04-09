#!/bin/bash

# ANSI barvy
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'  # No Color

success=0
failure=0

for (( n=2; n<=200; n++ )); do
    # Vytvoří nebo vyprázdní soubor numbers
    > numbers
    echo "Generuji $n náhodných 1-bajtových čísel do souboru 'numbers'..."
    
    for (( i=0; i<n; i++ )); do
        # Generuje číslo v rozmezí 0-255
        echo $(( RANDOM % 256 )) >> numbers
    done

    echo "Spouštím mpiexec s $n procesy..."
    # Spustí program a zachytí jeho výstup
    output=$(mpiexec  --oversubscribe  -n "$n" ./oets)

    echo "Výstup programu:"
    echo "$output"

    # Seřadíme výstup pomocí sort s timeoutem 5 sekund
    sorted=$(echo "$output" | timeout 5s sort -n)
    sort_exit_code=$?

    if [ $sort_exit_code -eq 124 ]; then
        echo -e "${RED}Seřazení vypršelo (timeout)!${NC}"
        ((failure++))
    elif [ "$output" = "$sorted" ]; then
        echo -e "${GREEN}Výstup je správně seřazený.${NC}"
        ((success++))
    else
        echo -e "${RED}Výstup není správně seřazený.${NC}"
        ((failure++))
    fi

    echo "-----------------------------------------"
done

echo "Shrnutí:"
echo -e "Úspěšných testů: ${GREEN}$success${NC}"
echo -e "Neúspěšných testů: ${RED}$failure${NC}"

