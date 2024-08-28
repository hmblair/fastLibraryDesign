// main.cpp

#include "library.h"

int main(int argc, char** argv) {
    // ensure there are exactly 2 arguments
    if (argc != 2) {
        std::cout << "Usage: ./main <filename>" << std::endl;
        return 1;
    }

    // get the filename from the command line arguments
    std::string filename = argv[1];

    // set the barcode stem
    std::string barcodeStemLoop = "UUCG";

    // set the five-prime and three-prime constant regions
    std::string fivePrimeConstantRegion = "ACTCGAGTAGAGTCGAAAA";
    std::string threePrimeConstantRegion = "AAAAGAAACAACAACAACAAC";

    // set the min and max stem length
    int minStemLength = 4;
    int maxStemLength = 16;

    // set the padding length
    int padToLength = 100;

    // set the barcode length
    int barcodeLength = 13;

    // set the final desired length of the sequences
    int finalLength = 170;

    // set the maximum number of each base pair in the barcode
    std::vector<int> maxBasePairCounts = {barcodeLength, 5, 1};

    // create a library object
    Library library = Library(
        filename,
        barcodeStemLoop
        );
    
    // print the length of the library
    std::cout << "Number of records: " << library.size() << std::endl;
    std::cout << "----------------------" << std::endl;

    // verify that all sequences are valid nucleic acid sequences
    library.verifyIsValidNucleicAcid();

    // add padding to the five prime end of the barcode
    library.padAllToLengthOnFivePrimeEnd(
        padToLength,
        minStemLength,
        maxStemLength,
        maxBasePairCounts
        );

    std::cout << "All sequences are padded to a length of " << padToLength << " nt." << std::endl;
    std::cout << "----------------------" << std::endl;

    // add barcodes to the library
    library.barcode(
        barcodeLength,
        maxBasePairCounts
        );

    std::cout << "----------------------" << std::endl;

    // remove the null barcode
    std::string nullBarcode = "N";
    int numBarcodesRemoved = library.removeBarcode(nullBarcode);

    std::cout << "Removed " << numBarcodesRemoved << " null barcodes from sequences." << std::endl;
    std::cout << "----------------------" << std::endl;

    // add constant regions to the library
    library.replaceFivePrimeConstantRegion(fivePrimeConstantRegion);
    library.replaceThreePrimeConstantRegion(threePrimeConstantRegion);

    // print the length discrepancy
    int lengthDiscrepancy = library.lengthDiscrepancy(finalLength);
    std::cout << "There are " << lengthDiscrepancy << " sequences that are not of the correct length, which is " << finalLength << std::endl;
    std::cout << "----------------------" << std::endl;

    // print the barcode discrepancy
    int barcodeDiscrepancy = library.barcodeDiscrepancy();
    std::cout << "There are " << barcodeDiscrepancy << " sequences without a unique barcode." << std::endl;
    std::cout << "----------------------" << std::endl;

    // convert to DNA
    library.toDNA();

    // write the library to a csv and a fasta file
    library.writeToCSV("output.csv");
    library.writeToFasta("output.fasta");

    return 0;
    
}