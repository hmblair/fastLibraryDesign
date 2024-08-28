// fasta.h

#include "constants.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>

// remove spaces at the beginning and end of a string
std::string strip(std::string s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) {
        return !std::isspace(c);
    }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {
        return !std::isspace(c);
    }).base(), s.end());
    return s;
}

bool validateNucleicSequence(std::string sequence) {
    for (char base : sequence) {
        if (NUCLEIC_BASES.find(std::string(1, base)) == NUCLEIC_BASES.end()) {
            return false;
        }
    }
    return true;
}

std::string toRNA(std::string sequence) {
    std::string rna = "";
    for (char base : sequence) {
        if (base == 'T') {
            rna += 'U';
        } else {
            rna += base;
        }
    }
    return rna;
}


std::string toDNA(std::string sequence) {
    std::string dna = "";
    for (char base : sequence) {
        if (base == 'U') {
            dna += 'T';
        } else {
            dna += base;
        }
    }
    return dna;
}


typedef struct {
    std::string header;
    std::string sequence;
} FastaRecord;


class FastaFile {
    public:
        std::vector<FastaRecord> records;

        // create an empty FastaFile
        FastaFile() {
            this->records = {};
        }

        // create a FastaFile from a file
        FastaFile(std::string filename) {
            this->records = this->read(filename);
        }

        // read a FastaFile from a file
        std::vector<FastaRecord> read(
            std::string filename,
            bool validate=true
            ) {

            // open the file, and validate that it opened correctly
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Unable to open file: " << filename << std::endl;
                exit(EXIT_FAILURE);
            }

            // get the length of the file
            file.seekg(0, std::ios::end);
            int fileLength = file.tellg();

            // reset the file to the beginning
            file.seekg(0, std::ios::beg);

            // read the file line by line
            std::string line;
            std::vector<FastaRecord> records;
            int numInvalidSequences = 0;

            while (std::getline(file, line)) {
                
                if (line[0] == '>') {

                    // create a new record
                    FastaRecord record;

                    // add the header to the record, and initialize the sequence
                    record.header = line;
                    record.sequence = "";

                    // add the record to the list of records
                    records.push_back(record);
                } else {

                    // strip whitespace from the line
                    line = strip(line);

                    // validate the sequence
                    if (validate && !validateNucleicSequence(line)) {
                        // count the number of invalid sequences
                        numInvalidSequences++;

                        // remove the last record and continue to the next line
                        records.pop_back();
                        continue;
                    }

                    // add the line to the last record
                    records.back().sequence += line;
                }
            }

            // print the number of invalid sequences
            std::cout << "----------------------------------------" << std::endl;
            std::cout << numInvalidSequences << " invalid sequences were ignored from " << filename << "." << std::endl;
            std::cout << "----------------------------------------" << std::endl;

            return records;
        }

        void write(std::string filename) {
            std::ofstream file(filename);
            for (FastaRecord record : this->records) {
                file << record.header << std::endl;
                file << record.sequence << std::endl;
            }
        }

        void push_back(FastaRecord record) {
            this->records.push_back(record);
        }

        void concatenate(FastaFile fastaFile) {
            for (FastaRecord record : fastaFile) {
                this->push_back(record);
            }
        }

        // attach a sequence to the 5' end of each record
        void attachToFivePrimeRegion(std::string sequence) {
            for (FastaRecord &record : this->records) {
                record.sequence = sequence + record.sequence;
            }
        }

        // attach a sequence to the 3' end of each record
        void attachToThreePrimeRegion(std::string sequence) {
            for (FastaRecord &record : this->records) {
                record.sequence = record.sequence + sequence;
            }
        }

        // convert all sequences to RNA
        void toRNA() {
            for (FastaRecord &record : this->records) {
                record.sequence = ::toRNA(record.sequence);
            }
        }

        // convert all sequences to DNA
        void toDNA() {
            for (FastaRecord &record : this->records) {
                record.sequence = ::toDNA(record.sequence);
            }
        }

        // get the unique lengths of the sequences in the file
        std::vector<int> getUniqueLengths() {
            std::vector<int> lengths;
            for (FastaRecord record : this->records) {
                int length = record.sequence.size();
                if (std::find(lengths.begin(), lengths.end(), length) == lengths.end()) {
                    lengths.push_back(length);
                }
            }

            // sort the lengths
            std::sort(lengths.begin(), lengths.end());

            return lengths;
        }

        // remove duplicate sequences from the file
        int removeDuplicates() {
            std::unordered_set<std::string> uniqueSequences;
            std::vector<FastaRecord> uniqueRecords;
            for (FastaRecord record : this->records) {
                if (uniqueSequences.find(record.sequence) == uniqueSequences.end()) {
                    uniqueSequences.insert(record.sequence);
                    uniqueRecords.push_back(record);
                }
            }

            // store the number of duplicates that were removed
            int numDuplicates = this->records.size() - uniqueRecords.size();

            // update the records to only include unique records
            this->records = uniqueRecords;

            // return the number of duplicates that were removed
            return numDuplicates;
        }

        std::vector<FastaRecord>::iterator begin() {
            return this->records.begin();
        }

        std::vector<FastaRecord>::iterator end() {
            return this->records.end();
        }

        int size() {
            return this->records.size();
        }

        FastaRecord operator[](int i) {
            return this->records[i];
        }

};