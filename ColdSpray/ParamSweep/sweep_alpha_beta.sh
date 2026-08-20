#!/bin/bash
# Sweeps the Robin-Robin contact parameters (alphaPC1, betaPC1, alphaC1P, betaC1P)
# through ColdSprayParameters, running ColdSprayOneSweep once per combination.
#
# Requires in this directory: ColdSprayGeometry, ColdSprayParameters, ColdSprayOne,
# iovtk.so (the compiled iovtk plugin), and ColdSprayOneSweepTemplate (a one-time
# copy of ColdSprayOne that includes "ColdSprayParametersSweep" instead of
# "ColdSprayParameters", has the dt-based alphaPC1/alphaC1P override commented out
# so the swept values actually take effect, and has two literal markers for this
# script to target: SWEEPNAME_PLACEHOLDER (the save-file name prefix) and
# OUTDIR_PLACEHOLDER (the output directory)).
#
# Each combination gets its own directory:
#   results_Params_APC_<alphaPC1>_BPC_<betaPC1>_ACP_<alphaC1P>_BCP_<betaC1P>/
# containing everything for that run together: the VTK frames, the Energy/
# Traction/vonMises CSVs (written directly there by FreeFem via OUTDIR_PLACEHOLDER),
# and the FreeFem++ console log.

# --- Values to sweep. Expand these lists as needed. ---
# All values use scientific notation (e.g. 1e0 for 1.0) - consistent with alpha,
# and avoids embedding decimal points in filenames in the first place.
alphaPC1_list="1e-6"
betaPC1_list="1e0"
alphaC1P_list="1e-6 1e-3"
betaC1P_list="1e0"

for alphaPC1 in $alphaPC1_list; do
    for betaPC1 in $betaPC1_list; do
        for alphaC1P in $alphaC1P_list; do
            for betaC1P in $betaC1P_list; do

                # ParaView's file-series auto-detection gets confused when the
                # "fixed" part of a filename contains an embedded decimal point
                # (e.g. it will misgroup ..._1.0_...0000.vtu, treating "0" and "1"
                # as different series). Values above use scientific notation
                # (1e0, 1e-9, ...) specifically to avoid this at the source, but
                # sanitize "." -> "p" here too as a safety net in case a value
                # with a literal decimal point (e.g. 1.5e0) is used instead.
                # so only the actual frame number varies within a series.
                alphaPC1_safe="${alphaPC1//./p}"
                betaPC1_safe="${betaPC1//./p}"
                alphaC1P_safe="${alphaC1P//./p}"
                betaC1P_safe="${betaC1P//./p}"

                dirtag="Params_APC_${alphaPC1_safe}_BPC_${betaPC1_safe}_ACP_${alphaC1P_safe}_BCP_${betaC1P_safe}"
                resultsDir="results_${dirtag}"
                savename="${dirtag}_ColdSprayOne"
                logfile="${resultsDir}/${savename}.log"

                echo ""
                echo "$resultsDir"
                if [[ ! -e $logfile ]]; then
                    mkdir -p "$resultsDir"

                    sed   "s/^real alphaPC1 = .*/real alphaPC1 = ${alphaPC1};/" ColdSprayParameters \
                    | sed "s/^real betaPC1  = .*/real betaPC1  = ${betaPC1};/" \
                    | sed "s/^real alphaC1P = .*/real alphaC1P = ${alphaC1P};/" \
                    | sed "s/^real betaC1P  = .*/real betaC1P  = ${betaC1P};/" \
                    > ColdSprayParametersSweep.idp

                    sed   "s/SWEEPNAME_PLACEHOLDER/${savename}/" ColdSprayOneSweepTemplate \
                    | sed "s|OUTDIR_PLACEHOLDER|${resultsDir}|" \
                    > ColdSprayOneSweep.edp

                    timeout 1200 FreeFem++ -nw ColdSprayOneSweep.edp > "$logfile" 2>&1
                    echo "  exit code: $?"
                else
                    echo "  already run"
                fi

            done
        done
    done
done

echo ""
echo "Sweep complete. Each combination's VTK, CSVs, and log all live together in:"
echo "  results_Params_APC_<alphaPC1>_BPC_<betaPC1>_ACP_<alphaC1P>_BCP_<betaC1P>/"
echo ""
echo "Summarize with, e.g.:"
echo "  grep -l 'Div by 0\\|Error\\|MATERROR\\|movemesh: .* reverse' results_Params_*/*.log   # crashed combinations"
echo "  grep 'Energy Start' results_Params_*/*.log | tail -1                                   # last logged energy per combination"
echo "  ls results_Params_*/*_Energy.csv                                                        # every combination's energy time series"
echo "  ls results_Params_*/*_Traction.csv                                                       # every combination's traction time series"
echo "  ls results_Params_*/*_vonMises.csv                                                       # every combination's von Mises time series"
