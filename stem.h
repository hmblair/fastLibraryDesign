// stem.h

#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <unordered_set>
#include <fstream>
#include <filesystem>






std::string generateRandomSequence(int length, std::vector<std::string> bases = RNA_BASES) {
    // create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // create a uniform distribution for the nucleic bases
    std::uniform_int_distribution<int> dist(0, 3);

    // create a vector to store the sequence
    std::string sequence = "";

    // loop over the length of the sequence and sample a random nucleic base
    for (int i = 0; i < length; i++) {
        sequence += bases[dist(rd)];
    }

    return sequence;
}


std::vector<int> sampleBitVector(int length) {
    // create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // create a uniform distribution for the nucleic bases
    std::uniform_int_distribution<int> dist(0, 1);

    // create a vector to store the sequence
    std::vector<int> sequence;

    // loop over the length of the sequence and sample a random nucleic base
    for (int i = 0; i < length; i++) {
        sequence.push_back(dist(rd));
    }

    return sequence;
}




std::vector<int> sampleVectorOfIntegersWithOccurenceConstraints(
    int length, 
    int maxValue, 
    std::vector<int> maxOccurences
    ) {

    // create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // create a vector containing maxOccurences[i] copies of i
    std::vector<int> values;
    for (int i = 0; i < maxOccurences.size(); i++) {
        for (int j = 0; j < maxOccurences[i]; j++) {
            values.push_back(i);
        }
    }

    // shuffle the values vector
    std::shuffle(values.begin(), values.end(), rd);

    // get the first length elements of the values vector
    std::vector<int> sampledValues(values.begin(), values.begin() + length);

    return sampledValues;
}



class Barcode {
    public:
        std::vector<int> basePairs;
        std::string stemLoop;

        Barcode(std::vector<int> basePairs, std::string stemLoop) {

            // initialise a barcode with the given base pairs and stem loop
            this->basePairs = basePairs;
            this->stemLoop = stemLoop;
        }

        Barcode(int length, std::vector<int> maxOccurences, std::string stemLoop) {

            // initialise a barcode with a random sequence of base pairs and
            // the given stem loop

            // generate a random sequence of 0s, 1, and 2s to represent the three
            // possible base pairs modulo orientation
            this -> basePairs = sampleVectorOfIntegersWithOccurenceConstraints(length, NUM_PAIRS_MODULO_ORIENTATION, maxOccurences);

            // create a random bit vector to represent the orientation of the base pairs, 
            // multiply the base pairs by two and add the bit vector to get the final base pairs
            std::vector<int> bitVector = sampleBitVector(length);
            for (int i = 0; i < length; i++) {
                this->basePairs[i] = 2 * this->basePairs[i] + bitVector[i];
            }

            // set the stem loop
            this->stemLoop = stemLoop;
        }


        std::string toString() {
            std::string stem = stemLoop;
            for (int i = 0; i < basePairs.size(); i++) {
                stem = RNA_PAIRS[basePairs[i]][0] + stem + RNA_PAIRS[basePairs[i]][1];
            }
            return stem;
        }


        std::vector<Barcode> hammingOneBall() {

            // generate the unit hamming ball of the stem barcode; that is, all
            // stem barcodes that are hamming distance one away from the current
            // stem barcode. To do this, loop over the elements of basePairs. 
            // For each element, we individually make the replacements

            // 0 -> 4
            // 1 -> 5
            // 2 -> 5
            // 3 -> 4
            // 4 -> 0
            // 4 -> 3
            // 5 -> 1
            // 5 -> 2

            // the final four require creating two new vectors rather than just one. 

            // create a vector to store the new stem barcodes
            std::vector<Barcode> stems;

            // loop over the elements of basePairs
            for (int i = 0; i < basePairs.size(); i++) {

                // make two copies of basePairs
                std::vector<int> newBasePairs = basePairs;
                std::vector<int> newBasePairs2 = basePairs;

                // make the replacements depending on the value of basePairs[i]
                if (basePairs[i] == 0) {
                    newBasePairs[i] = 4;
                    stems.push_back(Barcode(newBasePairs, this->stemLoop));
                } else if (basePairs[i] == 1) {
                    newBasePairs[i] = 5;
                    stems.push_back(Barcode(newBasePairs, this->stemLoop));
                } else if (basePairs[i] == 2) {
                    newBasePairs[i] = 5;
                    stems.push_back(Barcode(newBasePairs, this->stemLoop));
                } else if (basePairs[i] == 3) {
                    newBasePairs[i] = 4;
                    stems.push_back(Barcode(newBasePairs, this->stemLoop));
                } else if (basePairs[i] == 4) {
                    newBasePairs[i] = 0;
                    stems.push_back(Barcode(newBasePairs, this->stemLoop));
                    newBasePairs2[i] = 3;
                    stems.push_back(Barcode(newBasePairs2, this->stemLoop));
                } else if (basePairs[i] == 5) {
                    newBasePairs[i] = 1;
                    stems.push_back(Barcode(newBasePairs, this->stemLoop));
                    newBasePairs2[i] = 2;
                    stems.push_back(Barcode(newBasePairs2, this->stemLoop));
                } else {
                    std::cout << "Error: invalid base pair: " << basePairs[i] << std::endl;
                }

            }

            // add the original stem to the list of stems
            stems.push_back(*this);

            return stems;

        }

    
    bool verifyHammingDistance(std::unordered_set<std::string>& barcodesToAvoid) {

        // check if the current stem barcode has a hamming distance of at least 1
        // from all the stem barcodes in the set

        // get the hamming one ball
        std::vector<Barcode> hammingOneBall = this->hammingOneBall();

        // check if any of the elements of the hamming one ball are in the set
        for (int i = 0; i < hammingOneBall.size(); i++) {
            if (barcodesToAvoid.find(hammingOneBall[i].toString()) != barcodesToAvoid.end()) {
                return false;
            }
        }

        return true;
    }
};


std::string getPadding(
    int paddingRequired,
    int minStemLength,
    int maxStemLength,
    std::vector<int> maxOccurences,
    std::string stemLoop
    ) {
        // initialise a string to store the padding
        std::string padding;

        // while the padding required is greater than zero, add either random bases
        // or a stem barcode to the padding, depending on the value of paddingRequired
        while (paddingRequired > 0) {
            // if the padding required is less than the minimum stem length, add random bases
            if (paddingRequired < minStemLength) {
                padding += generateRandomSequence(paddingRequired);

                // set the padding required to zero
                paddingRequired = 0;
            } else {
                // otherwise, add a stem barcode of length min(maxStemLength, paddingRequired)
                int largestStemRequried = (paddingRequired - stemLoop.size()) / 2;
                int stemLength = std::min(maxStemLength, largestStemRequried);

                Barcode stemBarcode = Barcode(stemLength, maxOccurences, stemLoop);

                std::string stem = stemBarcode.toString();

                padding += stem;

                // update the padding required
                paddingRequired -= stem.size();

            }
        }

        return padding;

    }




std::string padSequence(
    std::string sequence, 
    int padToLength,
    int minStemLength,
    int maxStemLength,
    std::vector<int> maxOccurences,
    std::string stemLoop
    ) {
        // calculate the amount of padding needed
        int paddingRequired = padToLength - sequence.size();

        // if the padding is negative, give an error
        if (paddingRequired < 0) {
            std::cout << "Error: sequence is longer than the padding length" << std::endl;
            exit(EXIT_FAILURE);
        }

        // otherwise, loop until the padding is zero
        while (paddingRequired > 0) {
            // if the padding required is less than the minimum stem length, add random bases
            if (paddingRequired < minStemLength) {
                sequence += generateRandomSequence(paddingRequired);

                // set the padding required to zero
                paddingRequired = 0;
            } else {
                // otherwise, add a stem barcode of length min(maxStemLength, paddingRequired)
                int largestStemRequried = (paddingRequired - stemLoop.size()) / 2;
                int stemLength = std::min(maxStemLength, largestStemRequried);

                Barcode stemBarcode = Barcode(stemLength, maxOccurences, stemLoop);

                std::string stem = stemBarcode.toString();

                sequence += stem;

                // update the padding required
                paddingRequired -= stem.size();

            }
        }

        return sequence;
    }