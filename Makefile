.PHONY: all rldockscore rldockscore-debug rldockscore-profile rldockscore-test util clear

all:
	$(MAKE) rldockscore
	#$(MAKE) rldockscore-debug
	#$(MAKE) rldockscore-profile
	#$(MAKE) rldockscore-test
	$(MAKE) util
	chmod +x ./bin/*

rldockscore:
	mkdir -p bin
	sed 's/#define\ BUILD_MODE/#define\ MODE_NORMAL/g' ./src/config.template > ./src/config.h
	g++ ./src/main.cpp ./src/mol.cpp ./src/vec3d.cpp -o ./bin/rldockscore -I ./src/ -O3 -std=c++11 -static -static-libgcc -static-libstdc++ -lstdc++fs -lpthread
	rm ./src/config.h

rldockscore-debug:
	mkdir -p bin
	sed 's/#define\ BUILD_MODE/#define\ MODE_DEBUG/g' ./src/config.template > ./src/config.h
	g++ ./src/main.cpp ./src/mol.cpp ./src/vec3d.cpp -o ./bin/rldockscore_debug -I ./src/ -g -std=c++11 -static -static-libgcc -static-libstdc++ -lstdc++fs -lpthread
	rm ./src/config.h

rldockscore-profile:
	mkdir -p bin
	sed 's/#define\ BUILD_MODE/#define\ MODE_PROFILE/g' ./src/config.template > ./src/config.h
	g++ ./src/main.cpp ./src/mol.cpp ./src/vec3d.cpp -o ./bin/rldockscore_profile -I ./src/ -O3 -std=c++11 -lstdc++fs -lpthread -lprofiler
	rm ./src/config.h

rldockscore-test:
	mkdir -p bin
	sed 's/#define\ BUILD_MODE/#define\ MODE_TEST/g' ./src/config.template > ./src/config.h
	g++ ./src/main.cpp ./src/mol.cpp ./src/vec3d.cpp -o ./bin/rldockscore_test -I ./src/ -O3 -std=c++11 -static -static-libgcc -static-libstdc++ -lstdc++fs -lpthread
	rm ./src/config.h

util:
	cp ./util/openeye_get_info ./bin/
	cp ./util/check_atom_order ./bin/
	cp ./util/convert_rdock_pose ./bin/
	cp ./util/convert_vina_pose ./bin/

clear:
	rm -f bin/rldockscore
	rm -f bin/rldockscore_debug
	rm -f bin/rldockscore_profile
	rm -f bin/rldockscore_test
	rm -f bin/check_atom_order
	rm -f bin/openeye_get_info
	rm -f bin/convert_rdock_pose
	rm -r bin/convert_vina_pose
