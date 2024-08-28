// library.h

#include "fasta.h"
#include "stem.h"
#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <unordered_set>
#include <fstream>
#include <filesystem>


std::vector<std::string> splitByDelimiter(std::string s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token = "";
    bool insideQuotes = false;
    for (char c : s) {
        if (c == delimiter && !insideQuotes) {
            tokens.push_back(token);
            token = "";
        } else if (c == '"') {
            insideQuotes = !insideQuotes;
        } else {
            token += c;
        }
    }
    tokens.push_back(token);
    return tokens;
}

bool isValidNucleicAcid(std::string sequence) {
    std::unordered_set<char> nucleicAcidBases = {'A', 'T', 'C', 'G', 'U', 'N'};
    for (char base : sequence) {
        if (nucleicAcidBases.find(base) == nucleicAcidBases.end()) {
            return false;
        }
    }
    return true;
}


class LibrarySequence {
    public:

        // the various library constructs
        std::string fivePrimeConstantRegion;
        std::string fivePrimePadding;
        std::string designRegion;
        std::string threePrimePadding;
        std::string barcode;
        std::string threePrimeConstantRegion;

        // the name of the sequence in the library, and the name of the 
        // sublibrary it belongs to
        std::string name;

        LibrarySequence(
            const std::string& fivePrimeConstantRegion = "",
            std::string fivePrimePadding = "",
            std::string designRegion = "",
            std::string threePrimePadding = "",
            std::string barcode = "",
            const std::string& threePrimeConstantRegion = "",
            std::string name = ""
        ) {
            this->fivePrimeConstantRegion = fivePrimeConstantRegion;
            this->fivePrimePadding = fivePrimePadding;
            this->designRegion = designRegion;
            this->threePrimePadding = threePrimePadding;
            this->barcode = barcode;
            this->threePrimeConstantRegion = threePrimeConstantRegion;
            this->name = name;
        }

        void trimDesignReigionOnFivePrimeEnd(int length) {
            this->designRegion = this->designRegion.substr(length);
        }

        void trimDesignReigionOnThreePrimeEnd(int length) {
            this->designRegion = this->designRegion.substr(0, this->designRegion.size() - length);
        }

        void addThreePrimePadding(
            int paddingLength,
            int minStemLength,
            int maxStemLength,
            std::vector<int> maxOccurences,
            std::string barcodeStemLoop
            ) {
            this->threePrimePadding = getPadding(paddingLength, minStemLength, maxStemLength, maxOccurences, barcodeStemLoop);
        }


        void addFivePrimePadding(
            int paddingLength,
            int minStemLength,
            int maxStemLength,
            std::vector<int> maxOccurences,
            std::string barcodeStemLoop
            ) {
            this->fivePrimePadding = getPadding(paddingLength, minStemLength, maxStemLength, maxOccurences, barcodeStemLoop);
        }


        void addBarcode(
            int barcodeLength, 
            std::vector<int> maxOccurences, 
            std::unordered_set<std::string>& barcodes, 
            std::string barcodeStemLoop
            ) {

            // create a barcode object
            Barcode barcode = Barcode(barcodeLength, maxOccurences, barcodeStemLoop);

            // while the barcode has a hamming distance less than two from all
            // other barcodes, generate a new barcode
            while (!barcode.verifyHammingDistance(barcodes)) {
                barcode = Barcode(barcodeLength, maxOccurences, barcodeStemLoop);
            }

            // set the barcode of the library sequence
            this->barcode = barcode.toString();

            // add the barcode to the set of barcodes
            barcodes.insert(this->barcode);
        }




        std::string toString() {
            return this->fivePrimeConstantRegion + this->fivePrimePadding + this->designRegion + this->threePrimePadding + this->barcode + this->threePrimeConstantRegion;
        }

        std::string toSeparatedString() {
            std::string separator = " / ";
            return this->fivePrimeConstantRegion + separator + this->fivePrimePadding + separator + this->designRegion + separator + this->threePrimePadding + separator + this->barcode + separator + this->threePrimeConstantRegion;
        }

        int length() {
            return this->fivePrimeConstantRegion.size() + this->fivePrimePadding.size() + this->designRegion.size() + this->threePrimePadding.size() + this->barcode.size() + this->threePrimeConstantRegion.size();
        }

        int designRegionLength() {
            return this->designRegion.size();
        }

        int paddedDesignRegionLength() {
            return this->fivePrimePadding.size() + this->designRegion.size() + this->threePrimePadding.size();
        }

        void toRNA() {
            this->fivePrimeConstantRegion = ::toRNA(this->fivePrimeConstantRegion);
            this->fivePrimePadding = ::toRNA(this->fivePrimePadding);
            this->designRegion = ::toRNA(this->designRegion);
            this->threePrimePadding = ::toRNA(this->threePrimePadding);
            this->barcode = ::toRNA(this->barcode);
            this->threePrimeConstantRegion = ::toRNA(this->threePrimeConstantRegion);
        
        }

        void toDNA() {
            this->fivePrimeConstantRegion = ::toDNA(this->fivePrimeConstantRegion);
            this->fivePrimePadding = ::toDNA(this->fivePrimePadding);
            this->designRegion = ::toDNA(this->designRegion);
            this->threePrimePadding = ::toDNA(this->threePrimePadding);
            this->barcode = ::toDNA(this->barcode);
            this->threePrimeConstantRegion = ::toDNA(this->threePrimeConstantRegion);
        }

        bool verifyIsValidNucleicAcid() {
            return isValidNucleicAcid(this->fivePrimeConstantRegion) && isValidNucleicAcid(this->fivePrimePadding) && isValidNucleicAcid(this->designRegion) && isValidNucleicAcid(this->threePrimePadding) && isValidNucleicAcid(this->barcode) && isValidNucleicAcid(this->threePrimeConstantRegion);
        }


        void removeBarcode() {
            this->barcode = "";
        }

};



class Library {
    public:
        std::vector<LibrarySequence> librarySequnces;
        std::unordered_set<std::string> barcodes;
        std::string barcodeStemLoop;

        Library(
            std::vector<LibrarySequence> librarySequnces = {},
            std::unordered_set<std::string> barcodes = {},
            std::string barcodeStemLoop = ""
        ) {
            this->librarySequnces = librarySequnces;
            this->barcodes = barcodes;
            this->barcodeStemLoop = barcodeStemLoop;
        }


        Library(
            std::string pathToCSV,
            std::string barcodeStemLoop = ""
        ) {
            this->librarySequnces = this->readFromCSV(pathToCSV);

            // verify that every sequence has a design region
            for (LibrarySequence librarySequence : this->librarySequnces) {
                if (librarySequence.designRegion.size() == 0) {
                    std::cout << "Error: " << librarySequence.toSeparatedString() << " does not have a design region.\n";
                }
            }

            this->barcodes = {};
            int nonUniqueBarcodes = 0;
            int nullBarcodes = 0;

            // create an iterator to loop over the vector
            auto it = this->librarySequnces.begin();

            // loop over the vector using the iterator
            while (it != this->librarySequnces.end()) {
                LibrarySequence& librarySequence = *it;                
                // if the barcode is not empty, add it to the set of barcodes
                if (librarySequence.barcode.size() > 0) {

                    // if the barcode is not the null barcode, check if it is unique.
                    // if it is not unique, remove it from the sequence, and increment
                    // the number of non-unique barcodes
                    if (librarySequence.barcode != "N") {

                        if (librarySequence.barcode.size() != 30) {
                            std::cout << "Error: " << librarySequence.toSeparatedString() << " does not have a 30 nt barcode.\n";
                        }

                        if (this->barcodes.find(librarySequence.barcode) != this->barcodes.end()) {
                            librarySequence.removeBarcode();
                            nonUniqueBarcodes++;
                        } else {
                            this->barcodes.insert(librarySequence.barcode);
                        }
                    } else {
                        this->barcodes.insert(librarySequence.barcode);
                        nullBarcodes++;
                    }

                    
                }

                // increment the iterator
                it++;

            }

            // print the number of barcodes
            std::cout << "There were " << this->barcodes.size() - 1 + nonUniqueBarcodes << " existing non-null (N) barcodes. Of these, " << nonUniqueBarcodes << " were not unique and so were removed. Moreover, there were " << nullBarcodes << " null barcodes." << std::endl;

            // store the barcode stem loop
            this->barcodeStemLoop = barcodeStemLoop;

        }




        void replaceFivePrimeConstantRegion(std::string sequence) {
            for (LibrarySequence& librarySequence : *this) {
                librarySequence.fivePrimeConstantRegion = sequence;
            }
        }

        void replaceThreePrimeConstantRegion(std::string sequence) {
            for (LibrarySequence& librarySequence : *this) {
                librarySequence.threePrimeConstantRegion = sequence;
            }
        }


        void padAllToLengthOnFivePrimeEnd(
            int length, 
            int minStemLength, 
            int maxStemLength, 
            std::vector<int> maxOccurences
            ) {
            for (LibrarySequence& librarySequence : *this) {
                librarySequence.addFivePrimePadding(
                    length - librarySequence.paddedDesignRegionLength(),
                    minStemLength, 
                    maxStemLength, 
                    maxOccurences, 
                    this->barcodeStemLoop
                    );
            }

        }


        int removeBarcode(std::string barcode) {

            // a counter to keep track of the number of barcodes removed
            int numBarcodesRemoved = 0;

            // create an iterator to loop over the vector
            auto it = this->librarySequnces.begin();

            // loop over the vector using the iterator
            while (it != this->librarySequnces.end()) {
                LibrarySequence& librarySequence = *it;
                if (librarySequence.barcode == barcode) {
                    librarySequence.removeBarcode();
                    numBarcodesRemoved++;
                }
                it++;
            }

            // remove the barcode from the set of barcodes
            this->barcodes.erase(barcode);

            return numBarcodesRemoved;
        }


        void padAllToLengthOnThreePrimeEnd(
            int length, 
            int minStemLength, 
            int maxStemLength, 
            std::vector<int> maxOccurences
            ) {
            for (LibrarySequence& librarySequence : *this) {
                librarySequence.addThreePrimePadding(
                    length - librarySequence.paddedDesignRegionLength(),
                    minStemLength, 
                    maxStemLength, 
                    maxOccurences, 
                    this->barcodeStemLoop
                    );
            }
        }


        void barcode(
            int barcodeLength, 
            std::vector<int> maxOccurences
            ) {
            int n = 0;
            for (LibrarySequence& librarySequence : *this) {
                // if the sequence already has a barcode, skip it, else add one
                if (librarySequence.barcode.size() > 0) {
                    continue;
                } else {
                    // create a barcode object
                    Barcode barcode = Barcode(barcodeLength, maxOccurences, barcodeStemLoop);

                    // while the barcode has a hamming distance less than two from all
                    // other barcodes, generate a new barcode
                    while (!barcode.verifyHammingDistance(barcodes)) {
                        barcode = Barcode(barcodeLength, maxOccurences, barcodeStemLoop);
                    }

                    // set the barcode of the library sequence
                    librarySequence.barcode = barcode.toString();

                    // add the barcode to the set of barcodes
                    this->barcodes.insert(librarySequence.barcode);
                }

                n++;
                if (n % 100000 == 0) {
                    std::cout << "Added barcode to " << n << " sequences." << std::endl;
                }

            }
        }


        int barcodeDiscrepancy() {
            return this->librarySequnces.size() - this->barcodes.size();
        }

        int lengthDiscrepancy(int length) {
            int n = 0;
            for (LibrarySequence librarySequence : *this) {
                if (librarySequence.length() != length) {
                    n++;
                }
            }
            return n;
        }


        void writeToCSV(std::string filename) {
            std::ofstream file(filename);
            file << "Name,5' Constant Region,5' Padding,Design Region,3' Padding,Barcode,3' Constant Region\n";
            for (LibrarySequence librarySequence : *this) {
                file << librarySequence.name << "," << librarySequence.fivePrimeConstantRegion << "," << librarySequence.fivePrimePadding << "," << librarySequence.designRegion << "," << librarySequence.threePrimePadding << "," << librarySequence.barcode << "," << librarySequence.threePrimeConstantRegion << "\n";
            }
            file.close();
        }


        std::vector<LibrarySequence> readFromCSV(std::string filename) {
            std::ifstream file(filename);
            std::string line;
            std::getline(file, line);
            std::vector<LibrarySequence> librarySequnces;
            while (std::getline(file, line)) {
                std::vector<std::string> tokens = splitByDelimiter(line, ',');
                LibrarySequence librarySequence = LibrarySequence(
                    tokens[1],
                    tokens[2],
                    tokens[3],
                    tokens[4],
                    tokens[5],
                    tokens[6],
                    tokens[0]
                );

                librarySequnces.push_back(librarySequence);
            }

            return librarySequnces;

        }




        void writeToFasta(std::string filename) {
            std::ofstream file(filename);
            for (LibrarySequence librarySequence : *this) {

                // if the name does not start with a >, add it, otherwise, write
                // the name as is, since it is already in fasta format
                if (librarySequence.name[0] != '>') {
                    file << ">";
                }
                file << librarySequence.name << "\n";
                file << librarySequence.toString() << "\n";
            }
            file.close();
        }


        void toDNA() {
            for (LibrarySequence& librarySequence : *this) {
                librarySequence.toDNA();
            }
        }


        void toRNA() {
            for (LibrarySequence& librarySequence : *this) {
                librarySequence.toRNA();
            }
        }


        void verifyIsValidNucleicAcid() {
            for (LibrarySequence librarySequence : *this) {
                if (!librarySequence.verifyIsValidNucleicAcid()) {
                    std::cout << "Error: " << librarySequence.toSeparatedString() << " is not a DNA sequence.\n";
                }
            }
        }


        std::vector<LibrarySequence>::iterator begin() {
            return this->librarySequnces.begin();
        }

        std::vector<LibrarySequence>::iterator end() {
            return this->librarySequnces.end();
        }

        int size() {
            return this->librarySequnces.size();
        }
};