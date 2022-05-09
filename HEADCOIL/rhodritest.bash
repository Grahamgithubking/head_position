#!/bin/bash

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{print $1}')
    SESS=$(echo "$line" | awk -F, '{print $2}')
    echo "Subject is $SUBJ and session is $SESS"
done < "34_preterm.csv"