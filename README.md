This project covers 4 files： 1.halo2 with C-FFI 2. C version NTT 3. C version MSM 4. C version setup. Each file can be run independently.

The file halo2 was cloned from [Zcash](https://github.com/zcash/halo2). We integrate it with our operators in C language to demonstrate how to call external code in Rust, for more details please check ```/halo2_proofs/src/arithmetics.rs``` and ```/halo2_proofs/src/build.rs```.  

To run a demo, use ```cargo test --package halo2_proofs --test plonk_api -- plonk_api --exact --show-output```  

The other files are the source code of our C language operators. To run our operator, use:

```
#for NTT  
g++ main.cpp fr.cpp arithmetics.cpp -o main.out
#for MSM
g++ main.cpp arithmetics.cpp fr.cpp curve.cpp fields.cpp -o main.out
#for setup
g++ main.cpp setup.cpp fq2-bn256.cpp curve.cpp fields.cpp -o main.out
#then
./main.out
```
The NTT code uses random numbers as input and requires the user to enter the data size they want to test. We provide a small set of input($$2^5$$) as test data for MSM and setup, they support larger scale but require you to generate input data.
