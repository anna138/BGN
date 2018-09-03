# BGN
BGN Cryptosystem for Research

Compile Mult Test
g++ -std=c++0x -g -o multRun mainmulttest.cpp ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp -lgmp -lsodium -lm -I ~/.local/include/pbc -L ~/.local/lib -Wl,-rpath ~/.local/lib -l pbc

Compile Add Test
g++ -std=c++0x -g -o addRun mainaddtest.cpp ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp -lgmp -lsodium -lm -I ~/.local/include/pbc -L ~/.local/lib -Wl,-rpath ~/.local/lib -l pbc

Compile Main:
g++ -std=c++0x -g -o runMain main.cpp ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp -lgmp -lsodium -lm -I ~/.local/include/pbc -L ~/.local/lib -Wl,-rpath ~/.local/lib -l pbc
