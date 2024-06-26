#!/bin/bash

for file in *.sh; do
	if [ $file != "all.sh" ]; then
		./$file
	fi
done

