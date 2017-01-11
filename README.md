# countkmers
## External libraries
https://github.com/sparsehash/sparsehash

https://github.com/ArashPartow/bloom


## Compile
g++ -std=c++11 countkmers.cpp -o countkmers

## Usage
./countkmers filename kmersize topcount

## TODO
Encoding Nucleobase names as 2 bits since there are only 4 of them
