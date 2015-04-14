STEP 1. download the RQMC generator from the URL: http://www.uni-kl.de/AG-Heinrich/SamplePack.html, 
   and put all the files in the same folder.


STEP 2. compile it: g++ -O3 MainTest.cpp DigitalNetsBase2.C -o main 

NOTE: Some errors might take place in  #include <fstream.h> of DigitalNetsBase2.C 
and #include <iostream.h> of DigitalNetsBase2.h. 
If so, please remove the ".h" involved in the two include commands and add 
"using namespace std;" into the header file DigitalNetsBase2.h.

STEP 3. run main.exe

STEP 4. get the *.txt data 

NOTE THAT THE CODE HAS ONLY BEEN TESTED ON WINDOWS 64 BITS.

Please feel free to report bugs, problems, and suggestions for improvements to hezhijian87@gmail.com.