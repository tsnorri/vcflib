#ifndef __VARIANT_H
#define __VARIANT_H

#include <vector>
#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include <stack>
#include <queue>
#include <set>
#include <functional>
#include <cstdio>
#include <regex>
#include "split.h"
#include "join.h"
#include "tabix.hpp"
#include "SmithWatermanGotoh.h"
#include "disorder.h"
#include "ssw_cpp.h"
#include "convert.h"
#include "multichoose.h"
#include "Fasta.h"
extern "C" {
    #include "filevercmp.h"
}

namespace vcflib {

class Variant;

enum VariantFieldType { FIELD_FLOAT = 0
                      , FIELD_INTEGER
                      , FIELD_BOOL
                      , FIELD_STRING
                      , FIELD_UNKNOWN
                      };

enum VariantFieldNumber { ALLELE_NUMBER = -2
                        , GENOTYPE_NUMBER = -1
                        };

const int INDEX_NONE = -1;
const int NULL_ALLELE = -1;

VariantFieldType typeStrToFieldType(std::string& typeStr);
std::ostream& operator<<(std::ostream& out, VariantFieldType type);

typedef std::map<std::string, std::map<std::string, std::vector<std::string> > > Samples;
typedef std::vector<std::pair<int, std::string> > Cigar;

class VariantCallFile {

public:

    std::istream* file;
    tabixpp::Tabix* tabixFile;

    bool usingTabix;
    std::string vcf_header;


    std::string header;
    std::string line; // the current line
    std::string fileformat;
    std::string fileDate;
    std::string source;
    std::string reference;
    std::string phasing;
    std::map<std::string, VariantFieldType> infoTypes;
    std::map<std::string, int> infoCounts;
    std::map<std::string, VariantFieldType> formatTypes;
    std::map<std::string, int> formatCounts;
    std::vector<std::string> sampleNames;
    bool parseSamples;
    bool _done;

    void updateSamples(std::vector<std::string>& newSampleNames);
    std::string headerWithSampleNames(std::vector<std::string>& newSamples); // non-destructive, for output
    void addHeaderLine(std::string line);
    void removeInfoHeaderLine(std::string line);
    void removeGenoHeaderLine(std::string line);
    std::vector<std::string> infoIds(void);
    std::vector<std::string> formatIds(void);
    void reset(void);

    bool open(std::string& filename) {
        std::vector<std::string> filenameParts = split(filename, ".");
        if (filenameParts.back() == "gz" || filenameParts.back() == "bgz") {
            return openTabix(filename);
        } else {
            return openFile(filename);
        }
    }

    bool openFile(std::string& filename) {
        file = &_file;
        _file.open(filename.c_str(), std::ifstream::in);
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool openTabix(std::string& filename) {
        usingTabix = true;
        tabixFile = new tabixpp::Tabix(filename);
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool open(std::istream& stream) {
        file = &stream;
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool open(std::ifstream& stream) {
        file = &stream;
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool openForOutput(std::string& headerStr) {
        parsedHeader = parseHeader(headerStr);
        return parsedHeader;
    }

VariantCallFile(void) : usingTabix(false), parseSamples(true), justSetRegion(false), parsedHeader(false) { }
    ~VariantCallFile(void) {
        if (usingTabix) {
            delete tabixFile;
        }
    }

    bool is_open(void) { return parsedHeader; }

    bool eof(void) { return _file.eof(); }

    bool done(void) { return _done; }

    bool parseHeader(std::string& headerStr);

    bool parseHeader(void);

    bool getNextVariant(Variant& var);

    bool setRegion(std::string region);
    bool setRegion(std::string seq, long int start, long int end = 0);
    std::vector<std::string> getHeaderLinesFromFile();

private:
    bool firstRecord;
    bool justSetRegion;
    bool usingFile;
    std::ifstream _file;
    bool parsedHeader;

};

class VariantAllele {
    friend std::ostream& operator<<(std::ostream& out, VariantAllele& var);
    friend bool operator<(const VariantAllele& a, const VariantAllele& b);
    friend VariantAllele operator+(const VariantAllele& a, const VariantAllele& b);
public:
    std::string ref;
    std::string alt;
    std::string repr;
    long position;
    /* // TODO
    bool isSNP(void);
    bool isMNP(void);
    bool isInsertion(void);
    bool isDeletion(void);
    bool isIndel(void);
    */
    VariantAllele(std::string r, std::string a, long p)
        : ref(r), alt(a), position(p)
    {
        std::stringstream s;
        s << position << ":" << ref << "/" << alt;
        repr = s.str();
    }
};

class Variant {

    friend std::ostream& operator<<(std::ostream& out, Variant& var);
    
public:

    std::string sequenceName;
    long position;
    long zeroBasedPosition(void);
    std::string id;
    std::string ref;
    std::vector<std::string> alt;      // a list of all the alternate alleles present at this locus
    std::vector<std::string> alleles;  // a list all alleles (ref + alt) at this locus
                                  // the indicies are organized such that the genotype codes (0,1,2,.etc.)
                                  // correspond to the correct offest into the allelese vector.
                                  // that is, alleles[0] = ref, alleles[1] = first alternate allele, etc.
    std::string vrepr(void);  // a comparable record of the variantion described by the record
    std::set<std::string> altSet(void);  // set of alleles, rather than vector of them
    std::map<std::string, int> altAlleleIndexes;  // reverse lookup for alleles
    std::map<std::string, std::vector<VariantAllele> > parsedAlternates(bool includePreviousBaseForIndels = false,
                                                         bool useMNPs = false,
                                                         bool useEntropy = false,
                                                         float matchScore = 10.0f,
                                                         float mismatchScore = -9.0f,
                                                         float gapOpenPenalty = 15.0f,
                                                         float gapExtendPenalty = 6.66f,
                                                         float repeatGapExtendPenalty = 0.0f,
                                                         std::string flankingRefLeft = "",
                                                         std::string flankingRefRight = "");
    // the same output format as parsedAlternates, without parsing
    std::map<std::string, std::vector<VariantAllele> > flatAlternates(void);

    std::map<std::string, std::string> extendedAlternates(long int newPosition, long int length);

    // Convert a structural variant the canonical VCF4.2 format using a reference.
    // returns true if the variant is canonicalized, false otherwise.
    bool canonicalize_sv(FastaReference& ref, std::vector<FastaReference*> insertions, int interval_sz = -1);

    std::string originalLine; // the literal of the record, as read
    // TODO
    // the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j
    // std::vector<std::pair<int, int> > genotypes;  // indexes into the alleles, ordered as per the spec
    std::string filter;
    double quality;
    VariantFieldType infoType(std::string& key);
    std::map<std::string, std::vector<std::string> > info;  // vector<std::string> allows for lists by Genotypes or Alternates
    std::map<std::string, bool> infoFlags;
    VariantFieldType formatType(std::string& key);
    std::vector<std::string> format;
    std::map<std::string, std::map<std::string, std::vector<std::string> > > samples;  // std::vector<std::string> allows for lists by Genotypes or Alternates
    std::vector<std::string> sampleNames;
    std::vector<std::string> outputSampleNames;
    VariantCallFile* vcf;

    //void addInfoInt(std::string& tag, int value);
    //void addInfoFloat(std::string& tag, double value);
    //void addInfoString(std::string& tag, std::string& value);

    void removeAlt(std::string& altallele);

public:

    Variant() { }

    Variant(VariantCallFile& v)
        : sampleNames(v.sampleNames)
        , outputSampleNames(v.sampleNames)
        , vcf(&v)
    { }

    void setVariantCallFile(VariantCallFile& v);
    void setVariantCallFile(VariantCallFile* v);

    void parse(std::string& line, bool parseSamples = true);
    void addFilter(std::string& tag);
    bool getValueBool(std::string& key, std::string& sample, int index = INDEX_NONE);
    double getValueFloat(std::string& key, std::string& sample, int index = INDEX_NONE);
    std::string getValueString(std::string& key, std::string& sample, int index = INDEX_NONE);
    bool getSampleValueBool(std::string& key, std::string& sample, int index = INDEX_NONE);
    double getSampleValueFloat(std::string& key, std::string& sample, int index = INDEX_NONE);
    std::string getSampleValueString(std::string& key, std::string& sample, int index = INDEX_NONE);
    bool getInfoValueBool(std::string& key, int index = INDEX_NONE);
    double getInfoValueFloat(std::string& key, int index = INDEX_NONE);
    std::string getInfoValueString(std::string& key, int index = INDEX_NONE);
    void printAlt(std::ostream& out);      // print a comma-sep list of alternate alleles to an ostream
    void printAlleles(std::ostream& out);  // print a comma-sep list of *all* alleles to an ostream
    int getAltAlleleIndex(std::string& allele);
    void updateAlleleIndexes(void);
    void addFormatField(std::string& key);
    void setOutputSampleNames(std::vector<std::string>& outputSamples);
    std::map<std::pair<int, int>, int> getGenotypeIndexesDiploid(void);
    int getNumSamples(void);
    int getNumValidGenotypes(void);
    std::string getGenotype(std::string const& sample);
    bool isPhased(void);
    // TODO
    //void setInfoField(std::string& key, std::string& val);

private:

    std::string lastFormat;

};


// from BamTools
// RuleToken implementation

class RuleToken {

public:

    // enums
    enum RuleTokenType { OPERAND = 0
                       , NUMBER
                       , BOOLEAN_VARIABLE
                       , NUMERIC_VARIABLE
                       , STRING_VARIABLE
                       , AND_OPERATOR
                       , OR_OPERATOR
                       , ADD_OPERATOR
                       , SUBTRACT_OPERATOR
                       , MULTIPLY_OPERATOR
                       , DIVIDE_OPERATOR
                       , NOT_OPERATOR
                       , EQUAL_OPERATOR
                       , GREATER_THAN_OPERATOR
                       , LESS_THAN_OPERATOR
                       , LEFT_PARENTHESIS
                       , RIGHT_PARENTHESIS
                       };

    // constructor
    RuleToken(std::string token, std::map<std::string, VariantFieldType>& variables);
    RuleToken(void) 
        : type(BOOLEAN_VARIABLE)
        , state(false)
    { }

    // data members
    RuleTokenType type;
    std::string value;

    double number;
    std::string str;
    bool state;

    bool isVariable; // if this is a variable
    //bool isEvaluated; // when we evaluate variables

    RuleToken apply(RuleToken& other);

};

inline int priority(const RuleToken& token) {
    switch ( token.type ) {
        case ( RuleToken::MULTIPLY_OPERATOR )     : return 8;
        case ( RuleToken::DIVIDE_OPERATOR )       : return 8;
        case ( RuleToken::ADD_OPERATOR )          : return 7;
        case ( RuleToken::SUBTRACT_OPERATOR )     : return 7;
        case ( RuleToken::NOT_OPERATOR )          : return 6;
        case ( RuleToken::EQUAL_OPERATOR )        : return 5;
        case ( RuleToken::GREATER_THAN_OPERATOR ) : return 5;
        case ( RuleToken::LESS_THAN_OPERATOR )    : return 5;
        case ( RuleToken::AND_OPERATOR )          : return 4;
        case ( RuleToken::OR_OPERATOR )           : return 3;
        case ( RuleToken::LEFT_PARENTHESIS )      : return 0;
        case ( RuleToken::RIGHT_PARENTHESIS )     : return 0;
        default: cerr << "invalid token type" << endl; exit(1);
    }
}

inline bool isRightAssociative(const RuleToken& token) {
    return (token.type == RuleToken::NOT_OPERATOR ||
            token.type == RuleToken::LEFT_PARENTHESIS);
}

inline bool isLeftAssociative(const RuleToken& token) {
    return !isRightAssociative(token);
}

inline bool isLeftParenthesis(const RuleToken& token) {
    return ( token.type == RuleToken::LEFT_PARENTHESIS );
}

inline bool isRightParenthesis(const RuleToken& token) {
    return ( token.type == RuleToken::RIGHT_PARENTHESIS );
}

inline bool isOperand(const RuleToken& token) {
    return ( token.type == RuleToken::OPERAND || 
             token.type == RuleToken::NUMBER ||
             token.type == RuleToken::NUMERIC_VARIABLE ||
             token.type == RuleToken::STRING_VARIABLE ||
             token.type == RuleToken::BOOLEAN_VARIABLE
           );
}

inline bool isOperator(const RuleToken& token) {
    return ( token.type == RuleToken::AND_OPERATOR ||
             token.type == RuleToken::OR_OPERATOR  ||
             token.type == RuleToken::NOT_OPERATOR ||
             token.type == RuleToken::EQUAL_OPERATOR ||
             token.type == RuleToken::GREATER_THAN_OPERATOR ||
             token.type == RuleToken::LESS_THAN_OPERATOR ||
             token.type == RuleToken::MULTIPLY_OPERATOR ||
             token.type == RuleToken::DIVIDE_OPERATOR ||
             token.type == RuleToken::ADD_OPERATOR ||
             token.type == RuleToken::SUBTRACT_OPERATOR
             );
}

inline bool isOperatorChar(const char& c) {
    return (c == '!' ||
            c == '&' ||
            c == '|' ||
            c == '=' ||
            c == '>' ||
            c == '<' ||
            c == '*' ||
            c == '/' ||
            c == '+' ||
            c == '-');
}

inline bool isParanChar(const char& c) {
    return (c == '(' || c == ')');
}

inline bool isNumeric(const RuleToken& token) {
    return token.type == RuleToken::NUMERIC_VARIABLE;
}

inline bool isString(const RuleToken& token) {
    return token.type == RuleToken::STRING_VARIABLE;
}

inline bool isBoolean(const RuleToken& token) {
    return token.type == RuleToken::BOOLEAN_VARIABLE;
}

inline bool isVariable(const RuleToken& token) {
    return isNumeric(token) || isString(token) || isBoolean(token);
}

void tokenizeFilterSpec(std::string& filterspec, stack<RuleToken>& tokens, map<std::string, VariantFieldType>& variables);


class VariantFilter {

public:

    enum VariantFilterType { SAMPLE = 0,
                             RECORD };

    std::string spec;
    std::queue<RuleToken> tokens; // tokens, infix notation
    std::queue<RuleToken> rules;  // tokens, prefix notation
    VariantFilterType type;
    VariantFilter(std::string filterspec, VariantFilterType filtertype, std::map<std::string, VariantFieldType>& variables);
    bool passes(Variant& var, std::string& sample); // all alts pass
    bool passes(Variant& var, std::string& sample, std::string& allele);
    void removeFilteredGenotypes(Variant& var, bool keepInfo);

};


// genotype manipulation

// TODO
//std::map<std::string, int> decomposeGenotype(std::string& genotype);

std::vector<int> decomposePhasedGenotype(const std::string& genotype);
std::map<int, int> decomposeGenotype(const std::string& genotype);

std::string genotypeToString(const std::map<int, int>& genotype);

std::string phasedGenotypeToString(const std::vector<int>& genotype);

bool isHet(const std::map<int, int>& genotype);

bool isHom(const std::map<int, int>& genotype);

bool hasNonRef(const std::map<int, int>& genotype);

bool isHomRef(const std::map<int, int>& genotype);

bool isHomNonRef(const std::map<int, int>& genotype);

bool isNull(const std::map<int, int>& genotype);

int ploidy(const std::map<int, int>& genotype);

std::string unionInfoHeaderLines(std::string& s1, std::string& s2);

// genotype likelihood ordering

std::list<std::list<int> > glorder(int ploidy, int alleles);
std::list<std::list<int> > _glorder(int ploidy, int alleles);
std::list<int> glsWithAlt(int alt, int ploidy, int numalts);
std::map<int, int> glReorder(int ploidy, int numalts, std::map<int, int>& alleleIndexMapping, std::vector<int>& altsToRemove);

std::vector<std::string>& unique(std::vector<std::string>& strings);

std::string varCigar(std::vector<VariantAllele>& vav, bool xForMismatch = false);
std::string mergeCigar(const std::string& c1, const std::string& c2);
std::vector<std::pair<int, std::string> > splitCigar(const std::string& cigarStr);
std::list<std::pair<int, std::string> > splitCigarList(const std::string& cigarStr);
int cigarRefLen(const std::vector<std::pair<int, char> >& cigar);
int cigarRefLen(const std::vector<std::pair<int, std::string> >& cigar);
std::vector<std::pair<int, std::string> > cleanCigar(const std::vector<std::pair<int, std::string> >& cigar);
std::string joinCigar(const std::vector<std::pair<int, std::string> >& cigar);
std::string joinCigar(const std::vector<std::pair<int, char> >& cigar);
std::string joinCigarList(const std::list<std::pair<int, std::string> >& cigar);
bool isEmptyCigarElement(const std::pair<int, std::string>& elem);

// for sorting, generating maps ordered by chromosome name
class ChromNameCompare {
public:
    bool operator()(const std::string& a, const std::string& b) const {
        return (filevercmp(a.c_str(), b.c_str()) < 0);
    }
};

class VCFHeader
{
public:
    VCFHeader();
    ~VCFHeader() {}

    /*
     * Adds header_column to this->header_columns if
     * it doesn't already exits.
     */
    void addHeaderColumn(const std::string& header_column);

    /*
     * Adds meta_line to either header_lines or header_lists.
     *
     * We parse out the ##_type_ from meta_line
     * - If the meta_line ##_type_ is a key in header_lines then meta_line is added to header_lines
     * - If the meta_line ##_type_ is a key in header_lists then meta_line is added to header_lists[##_type_] vector<std::string>
     *    Unless that header_lists[##_type_] vector already contains the ID that is in meta_line, in that case it is not added
     */
    void addMetaInformationLine(const std::string& meta_line);

    /*
     * Converts header_lines, header_lists and header_columns to a proper VCF header
     */
    std::string getHeaderString();

private:
    VCFHeader(const VCFHeader& vcfHeader); // Do not implement the copy constructor, there is no reason to add this functionality
    VCFHeader& operator=(const VCFHeader& vcfHeader); // Do not implement operator=, there is no reason to add this functionality

    /*
     * This is a helper function that determines if the ID substd::string contained in meta_line
     * exists as a ID substd::string within the vector<std::string> meta_lines. Returns true if
     * the ID exists within the vector and false otherwise.
     */
    bool metaInfoIdExistsInVector(const std::string& meta_line, std::vector<std::string>& meta_lines);

    /*
     * header_line_names_ordered contains all the header lines that
     * are available and in the expected order for a valid VCF file
     */
    std::vector<std::string> header_line_names_ordered;
    /*
     * header_list_names_ordered contains all the header lists that
     * are available and in the expected order for a valid VCF file
     */
    std::vector<std::string> header_list_names_ordered;

    /*
     * header_columns is set by the constructor to contain the 8 manditory VCF fields.
     * Also, unique header_columns for each of the vcf files are added as well.
     * Duplicates are not allowed, to prevent duplicates use addHeaderColumn when adding header columns
     */
    std::vector<std::string> header_columns;

    /* 
     * the maps we're going to be using will be case-insensitive
     * so that "fileFormat" and "fileformat" hash to the same item.
     */
    struct stringcasecmp : binary_function<std::string, std::string, bool> {
        struct charcasecmp : public std::binary_function<unsigned char, unsigned char, bool> {
            bool operator() (const unsigned char& c1, const unsigned char& c2) const {
                return tolower (c1) < tolower (c2); 
            }
        };
        bool operator() (const std::string & s1, const std::string & s2) const {
            return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end(), charcasecmp());
        }
    };

    // contains all the ##_types_ as keys, the value is either empty or a VCF file has set it
    std::map<std::string, std::string, stringcasecmp> header_lines; 

    // contains all the ##_types_ as keys, the value is a vector of ##_type_ (since there can be duplicate #INFO for example, duplicate ids are not allowed)
    std::map<std::string, vector<std::string>, stringcasecmp> header_lists; 

};

} // end namespace VCF

#endif
