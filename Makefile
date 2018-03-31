# Sets up all the tools necessary for the analyses
# with the exception of USEARCH:
# It has to be downloaded manually (http://www.drive5.com/usearch/)
# and the path to the binary has to be set accordingly in scripts/others.sh,
# scripts/analysis_quality.sh and scripts/analysis_subsamples_fixed.sh.



all: sampler swarm-v1 swarm-v2 gefast vsearch cdhit dnaclust sumaclust
	

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


### Download & build VSEARCH (v2.7.1)
vsearch:
	git clone https://github.com/torognes/vsearch.git tools/vsearch-repo
	git -C tools/vsearch-repo checkout f9960d249b5637d0704e2cc3b9e3fa71573c8ee3
	cd tools/vsearch-repo/; ./autogen.sh; ./configure; make
	cp tools/vsearch-repo/bin/vsearch tools/vsearch

### Download & build CD-HIT (v4.6.8)
cdhit:
	wget -P tools/ https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz
	cd tools; tar -zxvf cd-hit-v4.6.8-2017-1208-source.tar.gz
	cd tools/cd-hit-v4.6.8-2017-1208/; make
	cp tools/cd-hit-v4.6.8-2017-1208/cd-hit tools/cd-hit
	rm tools/cd-hit-v4.6.8-2017-1208-source.tar.gz

### Download & build DNACLUST (release 3)
dnaclust:
	wget -P tools/ https://downloads.sourceforge.net/project/dnaclust/parallel_release_3/dnaclust_linux_release3.zip
	unzip tools/dnaclust_linux_release3.zip -d tools/
	cp tools/dnaclust_linux_release3/dnaclust tools/dnaclust
	rm tools/dnaclust_linux_release3.zip

### Download & build Sumaclust (v1.0.31)
sumaclust:
	wget -P tools/ https://git.metabarcoding.org/obitools/sumaclust/uploads/59ff189079b9e318f07b9ff9d5fee54b/sumaclust_v1.0.31.tar.gz
	cd tools; tar -zxvf sumaclust_v1.0.31.tar.gz
	cd tools/sumaclust_v1.0.31/; make
	cp tools/sumaclust_v1.0.31/sumaclust tools/sumaclust
	rm tools/sumaclust_v1.0.31.tar.gz



### Remove all created files
.PHONY: clean
clean:
	rm -f tools/FastaSampler
	rm -rf tools/swarm-*
	rm -f tools/Swarm-*
	rm -rf tools/gefast
	rm -f tools/GeFaST
	rm -rf tools/cd-hit*
	rm -rf tools/dnaclust*
	rm -rf tools/sumaclust*
	rm -rf tools/vsearch*
