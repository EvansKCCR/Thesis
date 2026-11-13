# Integrin-adhesome workflow 1.0 environment helper.
# Run from this bundle root under WSL2 or Linux (including HPC login/compute nodes).

SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
.ONESHELL:

ENV_NAME ?= adhesome-workflow-v1
ENV_FILE ?= config/environment.yml
REPORT_DIR ?= validation_reports/deployment
USE_ACTIVE_ENV ?= 0
ENV_MANAGER ?= auto

HOST_PROFILE := $(shell \
	if [ -n "$$WSL_DISTRO_NAME" ] || grep -qi microsoft /proc/sys/kernel/osrelease 2>/dev/null || grep -qi microsoft /proc/version 2>/dev/null; then \
		printf WSL; \
	elif [ -n "$$SLURM_JOB_ID$$SLURM_JOBID$$PBS_JOBID$$PBS_ENVIRONMENT$$LSB_JOBID$$LSB_DJOB_NUMPROC$$SGE_JOB_ID$$NSLOTS$$LSF_ENVDIR$$MODULESHOME$$LMOD_CMD" ]; then \
		printf HPC; \
	else printf Linux; fi)
AUTO_THREADS := $(shell \
	for value in "$$SLURM_CPUS_PER_TASK" "$$PBS_NP" "$$PBS_NUM_PPN" "$$NCPUS" "$$NSLOTS" "$$LSB_DJOB_NUMPROC"; do \
		case "$$value" in ''|*[!0-9]*) ;; *) printf '%s' "$$value"; exit 0;; esac; \
	done; nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || printf 4)
THREADS ?= $(AUTO_THREADS)

.PHONY: help info prepare preflight env-create env-update

help:
	@printf '%s\n' \
	  'Integrin-adhesome workflow environment targets:' \
	  '  make info                         Detect WSL, Linux or HPC and show CPU defaults' \
	  '  make prepare                      Create the pinned Conda environment if absent, then preflight' \
	  '  make prepare USE_ACTIVE_ENV=1     Use the already-loaded environment (useful on managed HPC)' \
	  '  make preflight                    Validate tools/configuration without running analysis' \
	  '  make env-create                   Create adhesome-workflow-v1 from config/environment.yml' \
	  '  make env-update                   Explicitly reconcile an existing environment with the YAML'

info:
	@printf 'Detected environment: %s\n' '$(HOST_PROFILE)'
	printf 'Bundle root: %s\n' "$$PWD"
	printf 'Suggested thread budget: %s\n' '$(THREADS)'
	if [ '$(HOST_PROFILE)' = WSL ] && [[ "$$PWD" == /mnt/* ]]; then \
	  printf '%s\n' 'WSL note: this bundle is on a Windows-mounted path. For large FASTA/OrthoFinder runs, use a Linux-filesystem copy.'; \
	fi
	if [ '$(HOST_PROFILE)' = HPC ]; then \
	  printf '%s\n' 'HPC note: use one compute node and the scheduler CPU allocation; site module names are intentionally not guessed.'; \
	fi
	if [ '$(USE_ACTIVE_ENV)' = 1 ]; then \
	  printf 'Environment mode: current shell (%s)\n' "$${CONDA_DEFAULT_ENV:-no Conda env marker}"; \
	else \
	  printf 'Managed environment: %s (%s)\n' '$(ENV_NAME)' '$(ENV_FILE)'; \
	fi

prepare: info
	@set -euo pipefail
	if [ '$(USE_ACTIVE_ENV)' = 1 ]; then
	  echo 'Using the current environment; package installation is skipped.'
	else
	  manager='$(ENV_MANAGER)'
	  if [ "$$manager" = auto ]; then
	    if command -v mamba >/dev/null 2>&1; then manager=mamba;
	    elif command -v conda >/dev/null 2>&1; then manager=conda;
	    elif command -v micromamba >/dev/null 2>&1; then manager=micromamba;
	    else
	      echo 'No Conda-compatible environment manager found.' >&2
	      echo 'Install Miniforge/Mamba, or load your site-managed environment and use: make prepare USE_ACTIVE_ENV=1' >&2
	      exit 2
	    fi
	  fi
	  if ! command -v "$$manager" >/dev/null 2>&1; then
	    echo "Requested environment manager not found: $$manager" >&2; exit 2
	  fi
	  if "$$manager" env list | awk -v name='$(ENV_NAME)' '$$1 == name { found=1 } END { exit !found }'; then
	    echo "Environment '$(ENV_NAME)' already exists; preserving it. Use 'make env-update' to reconcile it explicitly."
	  else
	    echo "Creating '$(ENV_NAME)' from $(ENV_FILE) with $$manager."
	    "$$manager" env create -f '$(ENV_FILE)' -y
	  fi
	fi
	$(MAKE) preflight USE_ACTIVE_ENV='$(USE_ACTIVE_ENV)' ENV_NAME='$(ENV_NAME)' ENV_MANAGER='$(ENV_MANAGER)' REPORT_DIR='$(REPORT_DIR)'

env-create:
	@set -euo pipefail
	manager='$(ENV_MANAGER)'
	if [ "$$manager" = auto ]; then
	  if command -v mamba >/dev/null 2>&1; then manager=mamba;
	  elif command -v conda >/dev/null 2>&1; then manager=conda;
	  elif command -v micromamba >/dev/null 2>&1; then manager=micromamba;
	  else echo 'Install/load Conda, Mamba or Micromamba first.' >&2; exit 2; fi
	fi
	"$$manager" env create -f '$(ENV_FILE)' -y

env-update:
	@set -euo pipefail
	manager='$(ENV_MANAGER)'
	if [ "$$manager" = auto ]; then
	  if command -v mamba >/dev/null 2>&1; then manager=mamba;
	  elif command -v conda >/dev/null 2>&1; then manager=conda;
	  elif command -v micromamba >/dev/null 2>&1; then manager=micromamba;
	  else echo 'Install/load Conda, Mamba or Micromamba first.' >&2; exit 2; fi
	fi
	"$$manager" env update -n '$(ENV_NAME)' -f '$(ENV_FILE)' --prune

preflight:
	@set -euo pipefail
	if [ '$(USE_ACTIVE_ENV)' = 1 ]; then
	  python -B maintenance/validate_deployment.py --dry-run --clean-state --phylogeny --report-dir '$(REPORT_DIR)'
	else
	  manager='$(ENV_MANAGER)'
	  if [ "$$manager" = auto ]; then
	    if command -v mamba >/dev/null 2>&1; then manager=mamba;
	    elif command -v conda >/dev/null 2>&1; then manager=conda;
	    elif command -v micromamba >/dev/null 2>&1; then manager=micromamba;
	    else echo 'No environment manager found; use USE_ACTIVE_ENV=1 for a loaded site environment.' >&2; exit 2; fi
	  fi
	  "$$manager" run -n '$(ENV_NAME)' python -B maintenance/validate_deployment.py --dry-run --clean-state --phylogeny --report-dir '$(REPORT_DIR)'
	fi
	printf 'Preflight report: %s/%s and %s/%s\n' "$$PWD" '$(REPORT_DIR)' '$(REPORT_DIR)' 'deployment_audit.json'
