#include "popphyl.h"


#ifndef DATA_H
#define DATA_H

class Data{

 private:

  Switch      readNexus(string filename);
  Switch      readPhylip(string filename);
  Switch      readFasta(string filename);
  int**       data;
  int         Nstate;
  int         Nsite;
  int         Ntaxa;
  string*     SPname;
  Alphabet    alphabet;
  const char* charalphabet;
  string      filename;
  string* 	  listAcc;
  bool 		  takeAcc;

 public:

  bool verbose;

  Data();
  Data(Switch s);
  ~Data();
  Switch         initFrom(string file);
  Switch         initFrom(string file, datatype type);
  void           writePhylip(string oname);
  void           writeNexus(string oname);
  void           writeFasta(string oname);


  //*************// inline accessor //*************//
  const int      getNstate(){return Nstate;}
  const int      getNsite(){return Nsite;}
  const int      getNtaxa(){return Ntaxa;}
  Alphabet       getAlphabetName(){return alphabet;}
  const char*    getAlphabet(){return charalphabet;}
  string         getFileName(){return filename;}
  const string*  getTaxaList();
  int            getTaxaIndex(string name);
  string	     getAcc(int index);
  string         getPrimarySequence(int taxa);
  double*        getComposition();
  int            getNextATG(int seq, int pos);
  void           reverseComplement();
  void           changeName(int index, string name);
  int            getColumnOfStateN(string taxa, int staten);
  int            getStateNofColumn(string taxa, int col);
  void           removeTaxa(int index);
  void           removeSite(int index);
  void           concatenate(Data *ali);

  int operator() (int taxa, int site) const
  {
    if (taxa >= Ntaxa || site >= Nsite)
      throw ERROR_DAD_DATA_INDEX("Data::operator(): data index out of range\n");
    return data[taxa][site];
  };

  string getTaxaName(int taxa){
    if (taxa >= Ntaxa || taxa <0)
      throw ERROR_DAD_DATA_INDEX("Data::getTaxaName(): taxa index out of range\n");
    return SPname[taxa];
  }

  //*************// error message //*************//
  class ERROR_BAD_DATA_FORMAT  {   };
  class ERROR_DAD_DATA_INDEX {
  public:
    string msg;
    ERROR_DAD_DATA_INDEX(string m) : msg(m) {};
  };
};


#endif
