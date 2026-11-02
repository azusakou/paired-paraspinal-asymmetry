#!/usr/bin/env bash
set -euo pipefail

data_root="${1:-/Volumes/wang/paper16/ais_sc_plasma/raw}"
px_dir="${data_root}/proteomics/PXD024246"
mkdir -p "${px_dir}/metadata" "${px_dir}/quant"

curl -fL \
  'https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024246&outputMode=XML' \
  -o "${px_dir}/metadata/ProteomeXchange_PXD024246.xml"

curl -fL \
  'https://download.iprox.cn/IPX0002825000/IPX0002825001/Proteins.xml' \
  -o "${px_dir}/quant/Proteins.xml"

shasum -a 256 \
  "${px_dir}/metadata/ProteomeXchange_PXD024246.xml" \
  "${px_dir}/quant/Proteins.xml"

