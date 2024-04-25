#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <cctype>



typedef std::size_t HASH_INDEX_T;

struct MyStringHash {
    HASH_INDEX_T rValues[5] { 983132572, 1468777056, 552714139, 984953261, 261934300 };
    MyStringHash(bool debug = true)
    {
        if(false == debug){
            generateRValues();
        }
    }
    // hash function entry point (i.e. this is h(k))
    HASH_INDEX_T operator()(const std::string& k) const
    {
        // Add your code here
        int power = 0;
        int groupcount = 0;

        std::vector<HASH_INDEX_T> w(5,0); // storing 5 0~4 groups
        HASH_INDEX_T wvalue = 0;



        for (int i= k.size()-1; i>=0; i--) { // start from end of string -> front
            int idx = letterDigitToNumber(k[i]);
            // wvalue += 36^(6-power)*idx;
            wvalue += idx * std::pow(36,power);
            power++; // upping power value by 1
            std::cout << "character: " << k[i] << " index value " << idx << " current wvalue "<< wvalue << std::endl;


            if (power%6 == 0 || i==0) { // if we've reached the end (beginning) of the string, 
            // or we reached end of the group of 6
                int groupidx = 4-groupcount; // need to start from back 4 -> 3-> 2 -> 1
                w[groupidx] = wvalue;
                power = 0; // resetting power
                groupcount += 1;
                wvalue = 0; //resetting wvalue sum 

             }

        }


        //hashing
        HASH_INDEX_T hashvalue = 0;
        for (int i=0; i<w.size(); i++) { // summing hash value of each of w's elements
            // hashvalue += rand[i] * w[i]; // rand generator 
            hashvalue += rValues[i] * w[i]; // rand generator 
        }

        return hashvalue;



    }

    // A likely helper function is to convert a-z,0-9 to an integral value 0-35
    HASH_INDEX_T letterDigitToNumber(char letter) const
    {
        // Add code here or delete this helper function if you do not want it
        // first change all upper to lower case
        letter = tolower(letter); // correct library ?

        if (letter >= 'a' && letter <= 'z') {
            return letter - 'a';
        }
        else if (letter >= '0' && letter <= '9') {
            return (letter - '0') + 26; // 0~9 -> 26~35
        }

        return 0;
    }

    // Code to generate the random R values
    void generateRValues()
    {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 generator (seed);  // mt19937 is a standard random number generator

        // Simply call generator() [it has an operator()] to get another random number
        for(int i{ 0 }; i < 5; ++i)
        {
            rValues[i] = generator();
        }
    }
};

#endif


//1hr