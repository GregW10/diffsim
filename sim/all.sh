#!/bin/bash

for file in *.sh; do
	if [ $file != "all.sh" ]; then
		./$file
		if [ $? -eq 0 ]; then
			echo "Ran $file successfully."
		else
			echo "$file failed."
		fi
	fi
done

