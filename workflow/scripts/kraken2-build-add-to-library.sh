#!/bin/bash

for i in $(ls); do kraken2-build --add-to-library $i --db GTDB --threads 8; done


