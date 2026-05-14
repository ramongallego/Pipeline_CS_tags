# scripts/read_params.sh
param_file="$1"

if [[ ! -f "$param_file" ]]; then
  echo "ERROR: Parameter file not found: $param_file"
  exit 1
fi

source "$param_file"
echo "Reading analysis parameters from: $param_file"