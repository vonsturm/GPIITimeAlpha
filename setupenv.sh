#!/bin/bash

if [[ $HOSTNAME == *"gerda-login"* ]]; then

	export GERDA_PHASEII_DATA="/nfs/gerda5/gerda-data/blind/active"
	export GERDA_META_DATA="$GERDA_PHASEII_DATA/meta"

else
	export GERDA_PHASEII_DATA="/home/sturm/Gerda/Data/v03.00"
	export GERDA_META_DATA="/home/sturm/Programs/gerda-metadata"
fi

export MU_CAL="$GERDA_META_DATA/config/_aux/geruncfg"
export GERDA_DATA_SETS="$GERDA_META_DATA/data-sets/phy"

