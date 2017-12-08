# Sets up all the tools necessary for the analyses
# with the exception of USEARCH:
# It has to be downloaded manually (http://www.drive5.com/usearch/)
# and the path to the binary has to be set accordingly in
# scripts/analysis_quality.sh.



all: sampler swarm-v1 swarm-v2 gefast
	

### Compile FastaSampler
sampler:
	g++ -std=c++11 -O3 -o tools/FastaSampler tools/FastaSampler.cpp


### Download & build Swarm (v1)
swarm-v1:
	git clone https://github.com/torognes/swarm.git tools/swarm-v1
	git -C tools/swarm-v1 checkout e3eea08ffed9dbf3a93afdec41605e81c4602a7b
	cd tools/swarm-v1/; make
	cp tools/swarm-v1/swarm tools/Swarm-1.2.3


### Download & build Swarm (v2)
swarm-v2:
	git clone https://github.com/torognes/swarm.git tools/swarm-v2
	git -C tools/swarm-v2 checkout f131cc02dd3f7e2cbac58cbd16d52ddd53e777eb
	cd tools/swarm-v2/src/; make
	cp tools/swarm-v2/bin/swarm tools/Swarm-2.1.13


### Download & build GeFaST (branch used in analyses of BMC paper)
gefast:
	git clone https://github.com/romueller/gefast tools/gefast
	git -C tools/gefast checkout -b gefast-paper-analysis origin/gefast-paper-analysis
	cd tools/gefast/; make
	cp tools/gefast/build/GeFaST tools/GeFaST


### Remove all created files
.PHONY: clean
clean:
	rm -f tools/FastaSampler
	rm -rf tools/swarm-*
	rm -f tools/Swarm-*
	rm -rf tools/gefast
	rm -f tools/GeFaST