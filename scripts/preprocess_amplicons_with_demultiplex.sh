#!/usr/bin/env bash

# ==============================================================================
# Line-by-line amplicon preprocessing with optional inline demultiplexing
# ==============================================================================
# This script is intentionally explicit and editable.
#
# Supported modes:
#   INPUT_MODE="inline_single"
#       One multiplexed single-end FASTQ containing inline tags and primers.
#       Uses an eDNAFlow/OBITools-style barcode file,