#ifndef BEDREADER_H
#define BEDREADER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>
#include "IntervalTree.h"
#include "split.h"

std::string strip(std::string const& str, char const* separators = " \t") {
    std::string::size_type const first = str.find_first_not_of(separators);
    return (first == std::string::npos) ? string()
        : str.substr(first, str.find_last_not_of(separators) - first + 1);
}

void parseRegion(
    std::string& region,
    std::string& startSeq,
    int& startPos,
    int& stopPos) {

    size_t foundFirstColon = region.find(":");

    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        startSeq = region;
        startPos = 0;
        stopPos = -1;
    } else {
        startSeq = region.substr(0, foundFirstColon);
        string sep = "..";
        size_t foundRangeSep = region.find(sep, foundFirstColon);
        if (foundRangeSep == string::npos) {
            sep = "-";
            foundRangeSep = region.find("-", foundFirstColon);
        }
        if (foundRangeSep == string::npos) {
            startPos = atoi(region.substr(foundFirstColon + 1).c_str());
            // differ from bamtools in this regard, in that we process only
            // the specified position if a range isn't given
            stopPos = startPos + 1;
        } else {
            startPos = atoi(region.substr(foundFirstColon + 1, foundRangeSep - foundFirstColon).c_str());
            // if we have range sep specified, but no second number, read to the end of sequence
            if (foundRangeSep + sep.size() != region.size()) {
                stopPos = atoi(region.substr(foundRangeSep + sep.size()).c_str()); // end-exclusive, bed-format
            } else {
                //stopPos = reference.sequenceLength(startSeq);
                stopPos = -1;
            }
        }
    }
}

// stores the posiitional information of a bed target entry
class BedTarget {

public:

    std::string seq;  // sequence name
    int left;    // left position
    int right;   // right position, adjusted to 0-base
    std::string desc; // descriptive information, target name typically

    BedTarget(std::string s) {
        parseRegion(s, seq, left, right); 
    }

    BedTarget(std::string s, int l, int r, std::string d = "")
        : seq(s)
        , left(l)
        , right(r)
        , desc(d)
    { }

};


class BedReader {

    bool _isOpen;
    std::ifstream file;

public:

    bool isOpen(void) { return _isOpen; }

    std::vector<BedTarget> targets;
    std::map<std::string, IntervalTree<BedTarget*> > intervals; // intervals by reference sequence

    std::vector<BedTarget> entries(void) {

        std::vector<BedTarget> entries;

        if (!isOpen()) {
            cerr << "bed targets file is not open" << endl;
            exit(1);
        }

        std::string line;
        while (std::getline(file, line)) {
            std::vector<std::string> fields = split(line, " \t");
            BedTarget entry(strip(fields[0]),
                            atoi(strip(fields[1]).c_str()),
                            atoi(strip(fields[2]).c_str()),
                            (fields.size() >= 4) ? strip(fields[3]) : "");
            entries.push_back(entry);
        }

        return entries;

    }

    std::vector<BedTarget*> targetsContained(BedTarget& target) {
        std::vector<Interval<BedTarget*> > results;
        intervals[target.seq].findContained(target.left, target.right, results);
        std::vector<BedTarget*> contained;
        for (std::vector<Interval<BedTarget*> >::iterator r = results.begin(); r != results.end(); ++r) {
            contained.push_back(r->value);
        }
        return contained;
    }

    std::vector<BedTarget*> targetsOverlapping(BedTarget& target) {
        std::vector<Interval<BedTarget*> > results;
        intervals[target.seq].findOverlapping(target.left, target.right, results);
        std::vector<BedTarget*> overlapping;
        for (std::vector<Interval<BedTarget*> >::iterator r = results.begin(); r != results.end(); ++r) {
            overlapping.push_back(r->value);
        }
        return overlapping;
    }

BedReader(void)
	: _isOpen(false)
    { }

BedReader(std::string& fname)
	: _isOpen(false) {
        open(fname);
    }

    void addTargets(std::vector<BedTarget>& targets) {
        std::map<std::string, std::vector<Interval<BedTarget*> > > intervalsBySeq;
        for (std::vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t) {
            intervalsBySeq[t->seq].push_back(Interval<BedTarget*>(1 + t->left, t->right, &*t));
        }
        for (std::map<std::string, std::vector<Interval<BedTarget*> > >::iterator s = intervalsBySeq.begin(); s != intervalsBySeq.end(); ++s) {
            intervals[s->first] = IntervalTree<BedTarget*>(s->second);
        }
    }

    void open(const std::string& fname) {
        file.open(fname.c_str());
        _isOpen = true;
        targets = entries();
        std::map<std::string, std::vector<Interval<BedTarget*> > > intervalsBySeq;
        for (std::vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t) {
            intervalsBySeq[t->seq].push_back(Interval<BedTarget*>(1 + t->left, t->right, &*t));
        }
        for (std::map<std::string, std::vector<Interval<BedTarget*> > >::iterator s = intervalsBySeq.begin(); s != intervalsBySeq.end(); ++s) {
            intervals[s->first] = IntervalTree<BedTarget*>(s->second);
        }
    }

};

#endif

