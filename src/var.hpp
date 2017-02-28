// not to complicate the issue but I need a different variant object to handle populations. 

#ifndef __VAR_H
#define __VAR_H

#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>      
#include <stdlib.h>
#include "split.h"

class zvar{
public:

  std::string name;

  int npop;
  
  std::string seqid;
  long int pos;

  double nalt ;
  double nref ;
  double af   ;

  double alpha; 
  double beta ;

  virtual void loadPop(std::vector< std::map< std::string, std::vector<std::string> > >& group, std::string seqid, long int position) = 0;
  virtual void estimatePosterior() = 0 ;
  virtual ~zvar() = 0;
  void setPopName(std::string  popName);
  
};

class genotype : public zvar {
  
public:

  double nhomr;
  double nhoma;
  double nhet ;
  double ngeno;
  double fis  ;
  double hfrq ;
  
  std::vector<int> genoIndex;
  std::vector<std::string> gts ;
  std::vector< std::vector < double > > genoLikelihoods;
  std::vector< std::vector < double > > genoLikelihoodsCDF;

  virtual double unphred(std::map< std::string, std::vector<std::string> > & geno, int index) = 0; 
  virtual void loadPop(std::vector< std::map< std::string, std::vector<std::string> > >& group, std::string seqid, long int position);
  virtual ~genotype() = 0;
  void estimatePosterior();
  

};

class pooled : public zvar{
public:

  double ntot  ;
  double afsum ; 
    
  std::vector<double> nalts;
  std::vector<double> nrefs;
  std::vector<double> afs  ; 

  void loadPop(std::vector< std::map< std::string, std::vector<std::string> > >& group, std::string seqid, long int position);
  void estimatePosterior();

  ~pooled();

  double bound(double v);

  pooled(void);
  
};

class gt : public genotype{
public:
  gt(void);
  double unphred(std::map< std::string, std::vector<std::string> > & geno, int index);
  ~gt();
};

class gl : public genotype{
public:
  gl(void);
  double unphred(std::map< std::string, std::vector<std::string> > & geno, int index);
  ~gl();
};

class gp : public genotype{
public:
  gp(void);
  double unphred(std::map< std::string, std::vector<std::string> > & geno, int index);
  ~gp();
};


class pl : public genotype{
public:
  pl(void);
  double unphred(std::map< std::string, std::vector<std::string> > & geno, int index);
  ~pl();
}; 

#endif 
