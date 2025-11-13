.SILENT:

MAKEFLAGS += --no-print-directory

all: 
	@echo running all tests
	mkdir -p out/
	$(MAKE) stdout
	$(MAKE) vector
	$(MAKE) quaternion
	$(MAKE) matrix
	$(MAKE) transform
	$(MAKE) intersection

vector:
	g++ vector.cpp -o out/vector.out -I../bml -O2
	./out/vector.out

quaternion:
	g++ quaternion.cpp -o out/quaternion.out -I../bml -O2
	./out/quaternion.out

matrix:
	g++ matrix.cpp -o out/matrix.out -I../bml -O2
	./out/matrix.out

preformance:
	g++ preformance.cpp -o out/preformance.out -I../bml -O2
	./out/preformance.out

stdout:
	g++ stdout.cpp -o out/stdout.out -I../bml -O2
	./out/stdout.out

transform:
	g++ transform.cpp -o out/transform.out -I../bml -O2
	./out/transform.out

intersection:
	g++ intersection.cpp -o out/intersection.out -I../bml -O2
	./out/intersection.out

playground:
	mkdir -p out/
	g++ playground.cpp -o out/playground.out -I../bml -O2
	./out/playground.out

clean:
	rm out/matrix.out
	rm out/quaternion.out
	rm out/vector.out
	rm out/preformance.out
	rm out/stdout.out
	rm out/transform.out