// main.cpp

#include <argparse/argparse.hpp>
#include "library.h"

using std::string;

int main(int argc, char** argv) {

	argparse::ArgumentParser program("fastLibraryDesign");
	
	program.add_argument("filename");

	program.add_argument("--barcodeLength")
		.default_value(13)
		.scan<'d', int>();

	program.add_argument("--barcodeStemLoop")
		.default_value("UUCG");

	program.add_argument("--padToLength")
		.default_value(100)
		.scan<'d', int>();

	program.add_argument("--seq5")
		.default_value("ACTCGAGTAGAGTCGAAAA");	
	program.add_argument("--seq3")
		.default_value("AAAAGAAACAACAACAACAAC");		

	program.add_argument("--minStemLength")
		.default_value(4)		
		.scan<'d', int>();
	program.add_argument("--maxStemLength")
		.default_value(16)		
		.scan<'d', int>();

	try {
	  program.parse_args(argc, argv);
	}
	catch (const std::exception& err) {
	  std::cerr << err.what() << std::endl;
	  std::cerr << program;
	  std::exit(1);
	}

	string filename = program.get<string>("filename");
	
	int barcodeLength = program.get<int>("--barcodeLength");
	string barcodeStemLoop = program.get<string>("--barcodeStemLoop");

	int padToLength = program.get<int>("--padToLength");

	string fivePrimeConstantRegion = program.get<string>("--seq5");
	string threePrimeConstantRegion = program.get<string>("--seq3");

	int minStemLength = program.get<int>("--minStemLength");
	int maxStemLength = program.get<int>("--maxStemLength");


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
