# update *.h files
for file in include/*.h
do
 cp ../C++/$file include
done

# update *.cpp files. note ushuffle.c will not be updated
for file in src/*.cpp
do
 cp ../C++/$file src
done

cd src
make clean
make

