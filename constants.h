// constants.h

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

std::vector<std::string> N_BASES = {"A", "C", "G", "U"};
std::vector<std::string> R_BASES = {"A", "G"};
std::vector<std::string> Y_BASES = {"C", "U"};
std::vector<std::string> K_BASES = {"G", "U"};
std::vector<std::string> M_BASES = {"A", "C"};
std::vector<std::string> S_BASES = {"C", "G"};
std::vector<std::string> W_BASES = {"A", "U"};
std::vector<std::string> V_BASES = {"A", "C", "G"};
std::vector<std::string> D_BASES = {"A", "G", "U"};
std::vector<std::string> H_BASES = {"A", "C", "U"};
std::vector<std::string> B_BASES = {"C", "G", "U"};

std::unordered_map<std::string, std::vector<std::string> > BASES_DICT = {
    {"A", {"A"}},
    {"C", {"C"}},
    {"G", {"G"}},
    {"U", {"U"}},
    {"T", {"T"}},
    {"N", N_BASES},
    {"R", R_BASES},
    {"Y", Y_BASES},
    {"K", K_BASES},
    {"M", M_BASES},
    {"S", S_BASES},
    {"W", W_BASES},
    {"V", V_BASES},
    {"D", D_BASES},
    {"H", H_BASES},
    {"B", B_BASES},
    {"_", N_BASES}
};


std::string replacePolybaseWithRandomBase(std::string base) {
    if (BASES_DICT.find(base) == BASES_DICT.end()) {
        return base;
    }
    return BASES_DICT[base][rand() % BASES_DICT[base].size()];
}

std::string replaceAllPolybasesWithRandomBases(std::string sequence) {
    std::string newSequence = "";
    for (int i = 0; i < sequence.size(); i++) {
        newSequence += replacePolybaseWithRandomBase(sequence.substr(i, 1));
    }
    return newSequence;
}






std::unordered_set<std::string> NUCLEIC_BASES = {"A", "C", "G", "U", "T"};




std::unordered_set<std::string> DNA_BASES = {"A", "C", "G", "T"};

std::vector<std::string> RNA_BASES = {"A", "C", "G", "U"};
std::vector<std::vector<std::string> > RNA_PAIRS = {{"A", "U"}, {"U", "A"}, {"C", "G"}, {"G", "C"}, {"G", "U"}, {"U", "G"}};
std::vector<int> RNA_PAIRS_INDEX = {0, 1, 2, 3, 4, 5};

int NUM_BASES = RNA_BASES.size();
int NUM_PAIRS = RNA_PAIRS.size();
int NUM_PAIRS_MODULO_ORIENTATION = NUM_PAIRS / 2;

std::vector<std::vector<int> > EDIT_DISTANCES = {
    {0, 2, 2, 2, 1, 2},
    {2, 0, 2, 2, 2, 1},
    {2, 2, 0, 2, 2, 1},
    {2, 2, 2, 0, 1, 2},
    {1, 2, 2, 1, 0, 2},
    {2, 1, 1, 2, 2, 0}
};

std::unordered_map<char, char> COMPLEMENTS = {
    {'A', 'U'},
    {'U', 'A'},
    {'C', 'G'},
    {'G', 'C'},
};