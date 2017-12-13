#!/bin/bash

if [[ $HOSTNAME = *"gerda-login"* ]]; then

	export GERDA_PHASEII_DATA="/nfs/gerda5/gerda-data/blind/active"
	export GERDA_META_DATA="$GERDA_PHASEII_DATA/meta"

elif [[ $HOSTNAME = *"mpi-hd.mpg.de"* ]]; then

	export GERDA_PHASEII_DATA="/lfs/l3/gerda/Daq/data-phaseII/blind/active"
	export GERDA_META_DATA="$GERDA_PHASEII_DATA/meta"

elif [[ $HOSTNAME = "wildmint" ]]; then

	export GERDA_PHASEII_DATA="/home/sturm/Gerda/Data/v03.00"
	export GERDA_META_DATA="/home/sturm/Programs/gerda-metadata"

else
	echo $HOSTNAME
	echo "Gerda data and meta-data not found on this machine!"

fi

export MU_CAL="$GERDA_META_DATA/config/_aux/geruncfg"
export GERDA_DATA_SETS="$GERDA_META_DATA/data-sets/phy"
