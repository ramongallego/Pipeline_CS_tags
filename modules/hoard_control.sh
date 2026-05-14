# modules/hoard_control.sh
# Removes intermediate directories when HOARD=NO.
# Guards each rm with [[ -d ]] so missing dirs don't cause errors.

if [[ "${HOARD}" == "NO" ]]; then
    echo "HOARD=NO: cleaning up intermediate files..."
    for dir in demultiplexed noprimers midfiles swarm_input filtered; do
        target="${OUTPUT_DIR}/${dir}"
        if [[ -d "${target}" ]]; then
            rm -rf "${target}"
            echo "  Removed: ${target}"
        fi
    done
else
    echo "HOARD=${HOARD}: keeping intermediate files."
fi
