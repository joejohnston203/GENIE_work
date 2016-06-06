#!/bin/bash

. ~/scripts/atlas_commands.sh

scp_export_jpj13 runs/fort_tensor_comp/gen.noRPA . 
scp_export_jpj13 runs/fort_tensor_comp/gen.RPA . 
