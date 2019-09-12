default: cosmomc

Debug: cosmomc_debug
Release: cosmomc
cleanDebug: clean delete
cleanRelease: clean delete

rebuild: clean delete isitgr

cosmomc: BUILD ?= gcc
cosmomc_debug: BUILD ?= gcc
isitgr: BUILD ?= gcc

getdist: ./source/*.*90
	cd ./source && make getdist BUILD=$(BUILD)

cosmomc: ./source/*.*90 ./camb/*.*90
	cd ./source && make cosmomc BUILD=$(BUILD)

cosmomc_debug: ./source/*.*90 ./camb/*.*90
	cd ./source && make cosmomc_debug OUTPUT_DIR=Debug BUILD=$(BUILD)

camspec: ./source/*.*90 ./camb/*.*90
	cd ./source && make highL=../highL PLANCKLIKE=cliklike_CamSpec

clean:
	cd ./source && make clean

isitgr: ./source/*.*90 ./camb/*.*90 ./bdndz_code/*.*c
	cd ./source && make isitgr BUILD=$(BUILD)

all: isitgr getdist

delete:
	rm -f cosmomc
	rm -f cosmomc_debug
	rm -f getdist
