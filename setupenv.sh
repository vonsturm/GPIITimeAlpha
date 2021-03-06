#!/bin/bash

echo "Setting environment for $HOSTNAME ..."

if [[ $HOSTNAME = *"gerda-login"* ]]; then

	export GERDA_PHASEII_DATA="/nfs/gerda5/gerda-data/blind/active"
	export GERDA_META_DATA="$GERDA_PHASEII_DATA/meta"
	export GERDA_DATA_SETS="$GERDA_META_DATA/data-sets/phy"

elif [[ $HOSTNAME = *"mpi-hd.mpg.de"* ]]; then

	export GERDA_PHASEII_DATA="/lfs/l3/gerda/Daq/data-phaseII/blind/active"
	export GERDA_META_DATA="$GERDA_PHASEII_DATA/meta"
	export GERDA_DATA_SETS="/lfs/l2/gerda/Hades/Analysis/Users/sturm/BAT/GPIITimeAlpha/gerda-data-sets"

elif [[ $HOSTNAME = "wildmint" ]]; then

	export GERDA_PHASEII_DATA="/home/sturm/Gerda/Data/v03.00"
	export GERDA_META_DATA="/home/sturm/Programs/gerda-metadata"
	export GERDA_DATA_SETS="$GERDA_META_DATA/data-sets/phy"

else
	echo $HOSTNAME
	echo "Gerda data and meta-data not found on this machine!"

fi

export MU_CAL="$GERDA_META_DATA/config/_aux/geruncfg"

echo "  * GERDA_PHASEII_DATA = $GERDA_PHASEII_DATA"
echo "  * GERDA_META_DATA = $GERDA_META_DATA"
echo "  * GERDA_DATA_SETS = $GERDA_DATA_SETS"
echo "  * MU_CAL = $MU_CAL"
