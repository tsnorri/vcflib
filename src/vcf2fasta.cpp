#include "Variant.h"
#include "convert.h"
#include "join.h"
#include "split.h"
#include <set>
#include <getopt.h>
#include "Fasta.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace vcflib;

#define ALLELE_NULL -1


class SampleFastaFile {

public:

    ofstream fastafile;
    long int pos;
    string linebuffer;
    string filename;
    string seqname;
    int linewidth;

    void write(string sequence) {
        linebuffer += sequence;
        while (linebuffer.length() > linewidth) {
            fastafile << linebuffer.substr(0, linewidth) << endl;
            linebuffer = linebuffer.substr(linewidth);
        }
    }

    SampleFastaFile(void) { }

    void open(string& m_filename, string& m_seqname, int m_linewidth = 80) {
        filename = m_filename;
        seqname = m_seqname;
        pos = 0;
        linewidth = m_linewidth;
        if (fastafile.is_open()) fastafile.close();
        fastafile.open(filename.c_str());
        if (!fastafile.is_open()) {
            cerr << "could not open " << filename << " for writing, exiting" << endl;
            exit(1);
        }
        fastafile << ">" << seqname << endl;
    }

    ~SampleFastaFile(void) {
        if (fastafile.is_open()) {
            write(""); // flush
            fastafile << linebuffer << endl;
            fastafile.close();
        }
    }

};

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [file]" << endl
         << endl
         << "options:" << endl
         << "    -f, --reference REF     Use this reference when decomposing samples." << endl
         << "    -p, --prefix PREFIX     Affix this output prefix to each file, none by default" << endl
         << "    -P, --default-ploidy N  Set a default ploidy for samples which do not have information in the first record (2)." << endl
         << endl
         << "Outputs sample_seq:N.fa for each sample, reference sequence, and chromosomal copy N in [0,1... ploidy]." << endl;
        //<< "Impossible regions of haplotypes are noted with an error message.  The corresponding" << endl
        //<< "regions of the output FASTA files will be marked as N." << endl
    exit(0);
}

map<string, int>& getPloidies(Variant& var, map<string, int>& ploidies, int defaultPloidy=2) {
    for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
        int p = ploidy(decomposeGenotype(var.getGenotype(*s)));
        if (p == 0) ploidies[*s] = defaultPloidy;
        else        ploidies[*s] = p;
    }
    return ploidies;
}

void closeOutputs(map<string, map<int, SampleFastaFile*> >& outputs) {
    for (auto &vPair : outputs)
    {
        for (auto &mPair : vPair.second)
            delete mPair.second;
    }
    outputs.clear();
}

template <typename T>
void initOutputs(map<string, map<int, SampleFastaFile*> >& outputs, T const &sampleNames, string& seqName, map<string, int>& ploidies, string& prefix) {
    closeOutputs(outputs);
    for (auto const &sample : sampleNames)
    {
        map<int, SampleFastaFile*>& outs = outputs[sample];
        int p = ploidies[sample];
        for (int i = 0; i < p; ++i) {
            string name = prefix + sample + "_" + seqName + ":" + convert(i) + ".fasta";
            if (!outs[i]) {
                SampleFastaFile* fp = new SampleFastaFile;
                outs[i] = fp;
            }
            SampleFastaFile& f = *outs[i];
            f.open(name, seqName);
        }
    }
}

void vcf2fasta(VariantCallFile& variantFile, FastaReference& reference, string& outputPrefix, int defaultPloidy, string& nullAlleleString) {
    string lastSeq;
    map<string, map<int, SampleFastaFile*> > outputs;
    Variant var(variantFile);
    map<string, int> lastPloidies;
    
    // Handle samples in chunks.
    auto const &allSamples = variantFile.sampleNames;
    size_t chunkSize = 200;
    set<string> currentSamples;
    for (size_t i = 0, count = allSamples.size(); i < count; i += chunkSize)
    {
        long int lastPos=0, lastEnd=0;
        currentSamples.clear();
        lastSeq.clear();
        auto begin = allSamples.cbegin() + i;
        auto end = allSamples.cbegin() + std::min(i + chunkSize, count);
        std::copy(begin, end, std::inserter(currentSamples, currentSamples.end()));

        variantFile.reset();
        while (variantFile.getNextVariant(var)) {
            if (!var.isPhased()) {
                cerr << "variant " << var.sequenceName << ":" << var.position << " is not phased, cannot convert to fasta" << endl;
                exit(1);
            }
            map<string, int> ploidies;
            getPloidies(var, ploidies, defaultPloidy);
            if (var.sequenceName != lastSeq || lastSeq.empty()) {
                if (!lastSeq.empty()) {
                    string ref5prime = reference.getSubSequence(lastSeq, lastEnd, reference.sequenceLength(lastSeq)-lastEnd);
                    for (map<string, map<int, SampleFastaFile*> >::iterator s = outputs.begin(); s != outputs.end(); ++s) {
                        map<int, SampleFastaFile*>& f = s->second;
                        for (map<int, SampleFastaFile*>::iterator o = f.begin(); o != f.end(); ++o) {
                            o->second->write(ref5prime);
                        }
                    }
                }
                initOutputs(outputs, currentSamples, var.sequenceName, ploidies, outputPrefix);
                lastSeq = var.sequenceName;
                lastPos = 0;
            } else if (!lastPloidies.empty() && lastPloidies != ploidies) {
                cerr << "cannot handle mid-sequence change of ploidy" << endl;
                // in principle it should be possible...
                // it's a matter of representation, GFASTA anyone?
                exit(1);
            }
            lastPloidies = ploidies;
            if (var.position < lastEnd) {
                cerr << var.position << " vs " << lastEnd << endl;
                cerr << "overlapping or out-of-order variants at " << var.sequenceName << ":" << var.position << endl;
                exit(1);
            }
            // get reference sequences implied by last->current variant
            string ref5prime;
            if (var.position - 1 - lastEnd > 0) {
                ref5prime = reference.getSubSequence(var.sequenceName, lastEnd, var.position - 1 - lastEnd);
            }
            // write alt/ref seqs for current variant based on phased genotypes
            for (auto const &sample : currentSamples) {
                vector<int> gt = decomposePhasedGenotype(var.getGenotype(sample));
                // assume no-call == ref?
                if (gt.empty()) {
                    cerr << "empty genotype for sample " << sample << " at " << var.sequenceName << ":" << var.position << endl;
                    exit(1);
                }
                int i = 0;
                for (vector<int>::iterator g = gt.begin(); g != gt.end(); ++g, ++i) {
                    // @TCC handle uncalled genotypes (*g == NULL_ALLELE)
                    if( *g == NULL_ALLELE ){
                        if( nullAlleleString == "" ){
                            cerr << "empty genotype call for sample " << sample << " at " << var.sequenceName << ":" << var.position << endl;
                            cerr << "use -n option to set value to output for missing calls" << endl; 
                            exit(1);
                        }else{
                            outputs[sample].at(i)->write(nullAlleleString);
                        }
                    }else{
                        outputs[sample].at(i)->write(ref5prime+var.alleles.at(*g));
                    }
                }
            }
            lastPos = var.position - 1;
            lastEnd = lastPos + var.ref.size();
        }
        // write last sequences
        {
            string ref5prime = reference.getSubSequence(lastSeq, lastEnd, reference.sequenceLength(lastSeq)-lastEnd);
            for (map<string, map<int, SampleFastaFile*> >::iterator s = outputs.begin(); s != outputs.end(); ++s) {
                map<int, SampleFastaFile*>& f = s->second;
                for (map<int, SampleFastaFile*>::iterator o = f.begin(); o != f.end(); ++o) {
                    o->second->write(ref5prime);
                }
            }
        }
        closeOutputs(outputs);
        // outputs are closed by ~SampleFastaFile
    }
}

int main(int argc, char** argv) {

    VariantCallFile variantFile;
    string fastaFileName;
    int defaultPloidy;
    string outputPrefix;
    string nullAlleleString;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"reference", required_argument, 0, 'f'},
                {"prefix", required_argument, 0, 'p'},
                {"default-ploidy", required_argument, 0, 'P'},
                {"no-call-string", required_argument, 0, 'n'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hmf:p:P:n:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'f':
            fastaFileName = optarg;
            break;

        case 'h':
            printSummary(argv);
            break;

	    case 'p':
            outputPrefix = optarg;
            break;

        case 'P':
            defaultPloidy = atoi(optarg);
            break;

	    case 'n':
            nullAlleleString = optarg;
            break;

        case '?':
            printSummary(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    FastaReference reference;
    if (fastaFileName.empty()) {
        cerr << "a reference is required for haplotype allele generation" << endl;
        printSummary(argv);
        exit(1);
    }
    reference.open(fastaFileName);

    if (optind < argc) {
        string filename = argv[optind];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    vcf2fasta(variantFile, reference, outputPrefix, defaultPloidy, nullAlleleString);

    return 0;

}

