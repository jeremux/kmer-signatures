#include "classData.h"




Data::Data()
{
  // no warning for unused static function !
  printSwitch(getSwitch(0));
  verbose=dataVerbose;
  if(verbose) {
      cout << "Data::Data()\n";
      cout.flush();
  }  /**********************/
  Ntaxa=0;
  Nsite=0;
  data=NULL;
  SPname=NULL;
  filename="null";
  seqLength=NULL;
  takeAcc=No;
  if(verbose){
      cout << "end  Data::Data()\n";
      cout.flush();
  }
}

Data::Data(Switch s)
{
  printSwitch(getSwitch(0));
  verbose=dataVerbose;
  if(verbose) {
      cout << "Data::Data()\n";
      cout.flush();
  }  /**********************/
  Ntaxa=0;
  Nsite=0;
  data=NULL;
  SPname=NULL;
  filename="null";
  if(verbose){
      cout << "end  Data::Data()\n";
      cout.flush();
  }
  takeAcc = s;
  seqLength=NULL;
  listAcc=NULL;
}


Data::~Data()
{
  if(verbose){
      cout << "Data::~Data()\n";
      cout.flush();
  }  /**********************/
  if(Ntaxa>0)
    {
      for(int i=0;i<Ntaxa;i++)
        delete[] data[i];
      delete[] data;
      delete[] SPname;
      delete[] seqLength;
      if(takeAcc==Yes)
        delete[] listAcc;
    }

  if(verbose){
      cout << "end Data::~Data()\n";
      cout.flush();
  }
}


void Data::writeNexus(string oname)
{
  if(verbose){
      cout << "Data::writeNexus(" << oname <<")\n";
      cout.flush();
  }  /**********************/
  ofstream nex_os((oname).c_str());
  nex_os << "#NEXUS\n\n[]\nBegin data;\n    dimensions ntax=" << Ntaxa << " nchar=" << Nsite << ";\n    format datatype=";
  if(Nstate==4)
    nex_os << "dna";
  else
    nex_os << "protein";
  nex_os << " gap=-;\n    Matrix\n";
  for(int j=0;j<Ntaxa;j++)
    {
      nex_os << SPname[j] << " ";
      for(int i=0;i<Nsite;i++)
        {
          if(alphabet==dna)
            nex_os << DNAletters[data[j][i]];
          else
            nex_os << AminoAcids[data[j][i]];
        }
      nex_os << "\n";
    }
  nex_os << ";\nend;\n";
  nex_os.close();
  if(verbose){
      cout << "end Data::writeNexus()\n";
      cout.flush();
  }
}


void Data::writePhylip(string oname)
{
  if(verbose){
      cout << "Data::writePhylip(" << oname <<")\n";
      cout.flush();
  }  /**********************/
  ofstream nex_os(oname.c_str());
  nex_os << Ntaxa << " " << Nsite << "\n";
  for(int j=0;j<Ntaxa;j++)
    {
      nex_os << SPname[j] << " ";
      for(int i=0;i<Nsite;i++)
        {
          if(alphabet==dna)
            nex_os << DNAletters[data[j][i]];
          else
            nex_os << AminoAcids[data[j][i]];
        }
      nex_os << "\n";
    }
  nex_os.close();
  if(verbose){
      cout << "end Data::writePhylip()\n";
      cout.flush();
  }
}



void Data::writeFasta(string oname)
{
  if(verbose){
      cout << "Data::writeFasta(" << oname <<")\n";
      cout.flush();
  }  /**********************/
  ofstream nex_os(oname.c_str());
  for(int j=0;j<Ntaxa;j++)
    {
      nex_os << ">" << SPname[j] << "\n";
      for(int i=0;i<Nsite;i++)
        {
          if(alphabet==dna)
            nex_os << DNAletters[data[j][i]];
          else
            nex_os << AminoAcids[data[j][i]];
          if((i+1)%60==0)
            nex_os << "\n";
        }
      nex_os << "\n";
    }
  nex_os.close();
  if(verbose){
      cout << "end Data::writeFasta()\n";
      cout.flush();
  }
}


Switch Data::readNexus(string filename)
{
  if(verbose){
      cout << "Data::readNexus(" << filename <<")\n";
      cout.flush();
  } /**********************/
  string str, tmp;
  char line[MAXSITE];
  Switch ok=Yes;
  //unsigned x;
  int x;
  try {
      fstream map_is(filename.c_str(),ios::in);
      map_is.getline(line,MAXSITE);
      str=line;
      if(str.compare("#NEXUS")!=0)
        {
          throw "bad NEXUS format";
        } else {
            ok=No;
            while(!map_is.eof())
              {
                map_is.getline(line,MAXSITE);
                str=line;
                if(str.compare("Begin data;")==0)
                  {
                    ok=Yes;
                    break;
                  }
              }
            if(ok==No)
              {
                cerr << "reach EOF while waiting for \"Begin data;\"";
                cerr.flush();
                throw "bad NEXUS format";
              }
            map_is.getline(line,MAXSITE);
            str=line;
            x=str.find("=");
            while(x!=-1)
              {
                str.replace(x,1," = ");
                x=str.find("=",x+2);
              }
            x=str.find(";");
            str.replace(x,1," ;");
            istringstream str_is(str);
            str_is >> tmp; // dimensions ntax=20 nchar=1243;
            str_is >> tmp; // ntax
            str_is >> tmp; // =
            str_is >> Ntaxa;
            str_is >> tmp; // nchar
            str_is >> tmp; // =
            str_is >> Nsite;
            data=new int*[Ntaxa];
            SPname=new string[Ntaxa];
            for(int i=0;i<Ntaxa;i++)
              data[i]=new int[Nsite];
            map_is.getline(line,MAXSITE);
            str=line;
            x=str.find("=");
            while(x!=-1)
              {
                str.replace(x,1," = ");
                x=str.find("=",x+2);
              }
            x=str.find(";");
            str.replace(x,1," ;");
            istringstream str_is2(str);
            str_is2 >> tmp; //format datatype=protein gap=-;
            str_is2 >> tmp;
            str_is2 >> tmp;
            str_is2 >> tmp;
            if(tmp=="protein")
              {
                alphabet=protein;
                charalphabet=AminoAcids;
                Nstate=20;
              } else {
                  if(tmp=="dna")
                    {
                      alphabet=dna;
                      charalphabet=DNAletters;
                      Nstate=4;
                    } else {
                        cerr << "ERROR: unrocognize datatype \"" << tmp << "\" while reading nexus file " << filename << "\n";
                        cerr.flush();
                        throw "badformat";
                    }
              }
            map_is.getline(line,MAXSITE);
            str=line;
            istringstream str_is3(str);
            str_is3 >> tmp;
            if(tmp.compare("Matrix")!=0)
              {
                cerr << "expected to read \"Matrix\"\n";
                cerr.flush();
                throw "badformat";
              }
            for(int i=0;i<Ntaxa;i++)
              {
                map_is.getline(line,MAXSITE);
                str=line;
                istringstream str_is3(str);
                str_is3 >> SPname[i];
                str_is3 >> tmp;
                for(int j=0;j<Nsite;j++)
                  {
                    int k;
                    Switch seen=No;
                    if(alphabet==dna)
                      {
                        for(k=0;k<=Nstate;k++)
                          {
                            if(tmp.at(j)==DNAletters[k])
                              {
                                seen=Yes;
                                break;
                              }
                          }
                      }
                    if(alphabet==protein)
                      {
                        for(k=0;k<=Nstate;k++)
                          {
                            if(tmp.at(j)==AminoAcids[k])
                              {
                                seen=Yes;
                                break;
                              }
                          }
                      }
                    if(seen==No)
                      {
                        k--;
                        if(verbose){
                            cout << "WARNING: expected alphabet ";
                            if(alphabet==dna)
                              cout << "DNA";
                            else
                              cout << "Protein";
                            cout << " and read unknown state: " << tmp.at(j) << " at taxa " << SPname[i] <<"; state replaced by: " << charalphabet[k] << "\n";
                        }
                        ok=No;
                      }
                    data[i][j]=k;
                  }
              }
            map_is.close();
        }
  } catch( ... )      {
      if(verbose){
          cerr << "failed to read Nexus file: " << filename <<"\n";
          cerr.flush();
      }
      throw;
  }
  if(verbose){
      cout << "Data::readNexus()\n";
      cout.flush();
  }
  initLengthSeq();
  return ok;
}

Switch Data::readPhylip(string filename)
{
  if(verbose){
      cout << "Data::readPhylip(" << filename <<")\n";
      cout.flush();
  } /**********************/
  Switch ok=Yes;
  char line[MAXSITE];
  string tmp, str;
  int isdna=0;
  int notdna=0;
  try {
      ifstream tmp_is(filename.c_str());
      tmp_is >> tmp;
      Ntaxa=atoi(tmp.c_str());
      if(Ntaxa==0)
        throw "badformat";
      tmp_is >> tmp;
      Nsite=atoi(tmp.c_str());
      if(Nsite==0)
        throw "badformat";
      data=new int*[Ntaxa];
      SPname=new string[Ntaxa];
      for(int i=0;i<Ntaxa;i++)
        data[i]=new int[Nsite];
      Switch seen;
      tmp_is.getline(line,MAXSITE);
      for(int i=0;i<Ntaxa;i++)
        {
          tmp_is.getline(line,MAXSITE);
          str=line;
          istringstream str_is3(str);
          str_is3 >> SPname[i];
          str_is3 >> tmp;
          for(int j=0;j<Nsite;j++)
            {
              seen=No;
              for(int k=0;k<=4;k++)
                {
                  if(tmp.at(j)==DNAletters[k])
                    {
                      seen=Yes;
                      break;
                    }
                }
              if(seen==Yes)
                isdna++;
              else
                notdna++;
            }
        }
      if((isdna/double(isdna+notdna))>0.9)
        {
          alphabet=dna;
          charalphabet=DNAletters;
          Nstate=4;
        } else {
            alphabet=protein;
            charalphabet=AminoAcids;
            Nstate=20;
        }
      ifstream map_is(filename.c_str());
      map_is.getline(line,MAXSITE);
      for(int i=0;i<Ntaxa;i++)
        {
          map_is.getline(line,MAXSITE);
          str=line;
          istringstream str_is3(str);
          str_is3 >> tmp;
          str_is3 >> tmp;
          for(int j=0;j<Nsite;j++)
            {
              int k;
              seen=No;
              if(alphabet==dna)
                {
                  for(k=0;k<=Nstate;k++)
                    {
                      if(tmp.at(j)==DNAletters[k])
                        {
                          seen=Yes;
                          break;
                        }
                    }
                }
              if(alphabet==protein)
                {
                  for(k=0;k<=Nstate;k++)
                    {
                      if(tmp.at(j)==AminoAcids[k])
                        {
                          seen=Yes;
                          break;
                        }
                    }
                }
              if(seen==No)
                {
                  k--;
                  if(verbose){
                      cout << "WARNING: expected alphabet ";
                      if(alphabet==dna)
                        cout << "DNA";
                      else
                        cout << "Protein";
                      cout << " and read unknown state: " << tmp.at(j) << " at taxa " << SPname[i] <<"; state replaced by: " << charalphabet[k] << "\n";
                  }
                  ok=No;
                }
              data[i][j]=k;
            }
        }
      map_is.close();
  } catch(...)      {
      if(verbose){
          cerr << "failed to read phylip file: " << filename <<"\n";
          cerr.flush();
      }
      throw;
  }
  if(verbose){
      cout << "end Data::readPhylip()\n";
      cout.flush();
  }

  initLengthSeq();
  return ok;
}




Switch Data::readFasta(string filename)
{
  if(verbose)
  {
      cout << "Data::readFasta(" << filename <<")\n";
      cout.flush();
  }
  /**********************/
  Switch ok=Yes;
  const char * fileChar = filename.c_str();
  Switch seen;
  char line[MAXSITE];
  string str, tmp, field, spname;
  int nbbloc, nbsp=0, nbsite=0, maxnbsite=0, tmpnsite=0, sizebloc=0,indexAcc=0;
  int isdna=0,notdna=0, tmpsize;


  cout.flush();
  try {

	  cerr << "debut fstream " << fileChar << "\n";
      fstream map_is0(fileChar,ios::in);
      cerr << "fin fstream \n";


      map_is0.getline(line,MAXSITE);

      str=line;

      cout.flush();

      if(str.at(0)!='>')
        throw "baformat";


      cout.flush();

      map_is0.close();
      fstream map_is(filename.c_str(),ios::in);
      // compter nomber d especes et de sites


      while(!map_is.eof())
      {
    	  map_is.getline(line,MAXSITE);
    	  str=line;
    	  if(str.length()>0)
    	  {
    		  tmp=str.at(0);
    		  if(tmp==">")
    		  {
    			  if(tmpnsite!=0)
    			  {
    				  nbsite=tmpnsite;
    				  if(nbsite > maxnbsite)
    					  maxnbsite=nbsite;
    				  nbsite=0;
    				  tmpnsite=0;
    			  }
    			  nbsp++;
    		  } else {
    			  if(nbsite==0)
    			  {
    				  if(sizebloc==0)
    					  sizebloc=str.length();
    				  tmpnsite+=str.length();
    			  }
    			  tmpsize=str.length();
    			  for(int j=0;j<tmpsize;j++)
    			  {
    				  seen=No;
    				  for(int k=0;k<=4;k++)
    				  {
    					  if(str.at(j)==DNAletters[k])
    					  {
    						  seen=Yes;
    						  break;
    					  }
    					  if(str.at(j)==dnaletters[k])
    					  {
    						  seen=Yes;
    						  break;
    					  }
    				  }
    				  if(seen==Yes)
    					  isdna++;
    				  else
    					  notdna++;
    			  }
    		  }
    	  }
      }

      if((isdna/double(isdna+notdna))>0.9)
      {
    	  alphabet=dna;
    	  charalphabet=DNAletters;
    	  Nstate=4;
      } else {
    	  alphabet=protein;
    	  charalphabet=AminoAcids;
    	  Nstate=20;
      }
      if(tmpnsite > maxnbsite)
    	  maxnbsite=tmpnsite;
      //if(nbsp==1)
      //maxnbsite=tmpnsite;

      Ntaxa=nbsp;
      Nsite=maxnbsite;
//      cout << "Ntaxa " << Ntaxa << " " << Nsite  <<"(" << nbsite << " " << tmpnsite << ")\n";
//      cout.flush();

//      cout << "début alloc data et SPname \n ";
      data=new int*[Ntaxa];
      SPname=new string[Ntaxa];
//      cout << "fin alloc data et SPname \n ";

      if(takeAcc==Yes)
    	  listAcc=new string[Ntaxa];

//      cout << "début alloc data\n ";
      for(int i=0;i<Ntaxa;i++)
      {
//    	  cout << "alloc data[" << i << "] = new int[" << Nsite << "]\n";
    	  data[i]=new int[Nsite];
    	  for(int j=0;j<Nsite;j++)
    		  data[i][j]=Nstate;
      }
//      cout << "fin alloc data\n ";
      map_is.close();
      fstream map_is2(filename.c_str(),ios::in);
      int length;
      nbsp=-1;

//      cout << "debut lire sequence et espece\n";
//      cout.flush();
      // lire les sequences et nom d especes
      while(!map_is2.eof())
      {
    	  map_is2.getline(line,MAXSITE);
    	  str=line;
    	  length=str.length();
    	  if(length>0)
    	  {

    		  tmp=str.at(0);
    		  if(tmp!=">")
    		  {
    			  for(int i=0;i<length;i++)
    			  {
    				  int k;
    				  seen=No;
    				  if(alphabet==dna)
    				  {
    					  for(k=0;k<=(Nstate+1);k++)
    					  {
    						  if(str.at(i)==DNAletters[k])
    						  {
    							  seen=Yes;
    							  break;
    						  }
    						  if(str.at(i)==dnaletters[k])
    						  {
    							  seen=Yes;
    							  break;
    						  }
    					  }
    				  }
    				  if(alphabet==protein)
    				  {
    					  for(k=0;k<=(Nstate+1);k++)
    					  {
    						  if(str.at(i)==AminoAcids[k])
    						  {
    							  seen=Yes;
    							  break;
    						  }
    					  }
    				  }
    				  if(seen==No)
    				  {
    					  k--;
    					  if(verbose){
    						  cout << "WARNING: expected alphabet ";
    						  if(alphabet==dna)
    							  cout << "DNA";
    						  else
    							  cout << "Protein";
    						  cout << " and read unknown state: " << str.at(i) << " at taxa " << SPname[nbsp] <<"; state replaced by: " << charalphabet[k] << "\n";
    					  }
    					  ok=No;
    				  }
    				  data[nbsp][(sizebloc*nbbloc)+i]=k;
    			  }
    			  nbbloc++;
    		  }
    		  else
    		  {
    			  nbsp++;
    			  nbbloc=0;
    			  SPname[nbsp]=str.substr(1,str.length());
//			  cout << SPname[nbsp] << "\n";
//			  cout.flush();
    			  if(takeAcc==Yes)
    			  {
    				  int x = str.find( "__" ) + 2 ;
    				  listAcc[indexAcc++] = str.substr(x,str.length());

    			  }
    		  }
    	  }
      }
  } catch( ... )      {
	  if(verbose){
		  cerr << "failed to read Fasta file: " << filename <<"\n";
		  cerr.flush();
	  }
	  throw;
  }
  if(verbose){
	  cout << "end Data::readFasta()\n";
	  cout.flush();
  }


  //	for (int var = 0; var < Ntaxa; var++)
  //	{
  //		cout << "\n";
  //		for (int var2 = 0; var2 < Nsite; var2++)
  //		{
  //			cout << data[var][var2] << "\t";
  //		}
  //	}
  //	cout << "\n";
  initLengthSeq();
  return ok;
}


void Data::initLengthSeq()
{
  seqLength = new int[Ntaxa];
  for(int k=0;k<Ntaxa;k++)
    seqLength[k]=0;
  for(int i=0;i<Ntaxa;i++)
  {
      for(int j=0;j<Nsite;j++)
          if(data[i][j]<4 && data[i][j]>=0)
                  seqLength[i]+=1; /* TODO: changer lol */
  }
}

int Data::getLengthSeq(int i)
{
  return seqLength[i];
}

Switch Data::initFrom(string file, datatype type)
{
  if(verbose){
      cout << "Data::initfrom(" << file<< ",";
      if(type==Nexus)
        cout << "Nexus";
      else if(type==Phylip)
        cout << "Phylip";
      else if(type==Fasta)
        cout << "Fasta";
      else
        cout << "Unknown Format";
      cout <<")\n";
      cout.flush();
  } /**********************/
  Switch ok;
  switch(type)
  {
  case 0:
    ok=readNexus(file);
    break;
  case 1 :
    ok=readPhylip(file);
    break;
  case 2 :
    ok=readFasta(file);
    break;
  default:
    ok=No;
  }
  if(verbose){
      cout << "end Data::initfrom()\n";
      cout.flush();
  }
  return ok;
}

Switch Data::initFrom(string file)
{
  if(verbose){
      cout << "Data::initfrom(" << file<< ")\n";
      cout.flush();
  } /**********************/
  Switch ok;
  filename=file;
  try {
      if(verbose){
          cout << "Trying Nexus\n";
          cout.flush();
      }
      ok=readNexus(file);
  } catch(...)      {
      try {
          if(verbose){
              cout << "Trying Fasta\n";
              cout.flush();
          }
          ok=readFasta(file);
      } catch(...)      {
          try {
              if(verbose){
                  cout << "Trying Phylip\n";
                  cout.flush();
              }
              ok=readPhylip(file);
          } catch ( ... ) {
              cerr << "ERROR in Data::initfrom(" << file<< "), fail to read datafile in Nexus, Phylip or Fasta format\n";
              cerr.flush();
              throw ERROR_BAD_DATA_FORMAT();
          }
      }
  }
  if(verbose){
      cout << "end Data::initfrom()\n";
      cout.flush();
  }
  return ok;
}

const string* Data::getTaxaList()
{
  return SPname;
}

double* Data::getComposition()
{
  double *compo=new double[getNstate()];
  for(int i=0;i<getNstate();i++)
    compo[i]=0;
  int nstate=0;
  for(int i=0;i<getNtaxa();i++)
    {
      for(int j=0;j<getNsite();j++)
        {
          if((*this)(i,j)!=getNstate())
            {
              compo[(*this)(i,j)]+=1;
              nstate++;
            }
        }
    }
  for(int i=0;i<getNstate();i++)
    {
      compo[i]/=(double)nstate;
      //cout << (getAlphabet())[i] << "\t" << compo[i] << "\n";
    }
  return compo;
}

void Data::reverseComplement()
{
  if(getAlphabetName()!=dna)
    {
      // skip
    } else {
        for(int i=0;i<getNtaxa();i++)
          {
            int *reverse=new int[getNsite()];
            for(int j=0;j<getNsite();j++)
              {
                if(data[i][j]==4)
                  reverse[getNsite()-1-j]=4;
                if(data[i][j]==0)
                  reverse[getNsite()-1-j]=3;
                if(data[i][j]==3)
                  reverse[getNsite()-1-j]=0;
                if(data[i][j]==1)
                  reverse[getNsite()-1-j]=2;
                if(data[i][j]==2)
                  reverse[getNsite()-1-j]=1;

              }
            delete[] data[i];
            data[i]=reverse;
          }
    }

}


void Data::changeName(int index, string name)
{
  if(index<Ntaxa)
    SPname[index]=name;
}


int Data::getTaxaIndex(string name)
{
  int index=-1;
  for(int i=0;i<getNtaxa();i++)
    {
      if(getTaxaName(i)==name)
        {
          index=i;
          break;
        }
    }
  return index;
}

string Data::getPrimarySequence(int taxa)
{
  string seq="";
  for(int i=0;i<getNsite();i++)
    {
      if((*this)(taxa,i)!=getNstate())
        seq+=getAlphabet()[(*this)(taxa,i)];
    }
  return seq;
}

int Data::getNextATG(int seq, int pos)
{
  int nextpos=-1;
  if(getNstate()!=4)
    {
      cerr << "WARNING in Data::getNextATG, cannot find ATG codon in non DNA dataset " << getFileName() << "\n";
      cerr.flush();
      return -1;
    }
  int posend=getNsite();
  if(pos+2<getNsite())
    {
      //cout << "gogogo\n";
      //cout.flush();
      int codon[3];
      int lastpos;
      while(pos < posend)
        {
          lastpos=pos;
          //cout << "? look at " << pos << "\n";
          //cout.flush();
          codon[0]=4;
          codon[1]=4;
          codon[2]=4;
          Switch complete=Yes;
          while((*this)(seq,pos)==getNstate())
            {
              pos++;
              if(pos>=getNsite())
                {
                  pos=posend+1;
                  break;
                }
            }
          if((pos < getNsite()) && (pos < posend))
            {
              codon[0]=(*this)(seq,pos);
              pos++;
              if(pos<getNsite())
                {
                  while((*this)(seq,pos)==getNstate())
                    {
                      pos++;
                      if(pos>=getNsite())
                        {
                          pos=posend+1;
                          break;
                        }
                    }
                  if((pos < getNsite()) && (pos < posend))
                    {
                      codon[1]=(*this)(seq,pos);
                      pos++;
                      if(pos<getNsite())
                        {
                          while((*this)(seq,pos)==getNstate())
                            {
                              pos++;
                              if(pos>=getNsite())
                                {
                                  pos=posend+1;
                                  break;
                                }
                            }
                          if((pos < getNsite()) && (pos < posend))
                            {
                              codon[2]=(*this)(seq,pos);
                              pos++;
                            } else {
                                //cout << "!!! 1\n";

                                complete=No;
                                codon[2]=4;
                            } // end got third codon position
                        } else {
                            //cout << "!!! 2\n";
                            complete=No;
                        } // end trunc at codon position 3
                    } else {
                        //cout << "!!! 3\n";
                        complete=No;
                    } // end got second codon position
                } else {
                    //cout << "!!! 4\n";
                    complete=No;
                } //end trunc at codon position 2
            } else {
                //cout << "!!! 5\n";
                complete=No;
            } // end got first codon position
          //cout << "read "<< DNAletters[codon[0]] << DNAletters[codon[1]] << DNAletters[codon[2]] << "\n";
          cout.flush();


          /************* Modif Jeremy *************/
          /* complete not use warning: donc je l'utilise */
          /***************************************/

          if (complete==Yes)
            {
              complete=Yes;
            }
          if(codon[0]==0 && codon[1]==3 && codon[2]==2)
            {
              //cout << "atg at "<< lastpos <<"!!!\n";
              //cout.flush();
              nextpos=lastpos;
              break;
            } else {
                pos=lastpos+1;
            }
        }

    } else {
        return -1;
    }

  return nextpos;
}

int Data::getColumnOfStateN(string taxa, int staten)
{
  int taxaid=getTaxaIndex(taxa);
  int seenstate=0;
  int pos=0;
  for(int i=0;i<getNsite();i++)
    {
      if((*this)(taxaid,i)!=getNstate())
        seenstate++;
      if(seenstate==staten)
        {
          pos=i;
          break;
        }
    }

  return pos;
}

int Data::getStateNofColumn(string taxa, int col)
{
  //int taxaid=getTaxaIndex(taxa);

  return 1;
}

void Data::removeTaxa(int index)
{
  int ntaxa=getNtaxa();
  if(index >= 0 && index<ntaxa)
    {
      delete[] data[index];
      data[index]=NULL;
      string *newsp=new string[ntaxa-1];
      int **newdata=new int*[ntaxa-1];
      for(int i=0;i<index;i++)
        {
          newsp[i]=SPname[i];
          newdata[i]=data[i];
        }
      for(int i=index;i<(ntaxa-1);i++)
        {
          newsp[i]=SPname[i+1];
          newdata[i]=data[i+1];
        }
      delete[] data;
      delete[] SPname;
      data=newdata;
      SPname=newsp;
      Ntaxa-=1;
    }
}


void Data::removeSite(int index)
{
  if(index<getNsite())
    {
      for(int i=0;i<getNtaxa();i++)
        {
          int *tmp=new int[getNsite()-1];
          for(int j=0;j<index;j++)
            tmp[j]=data[i][j];
          for(int j=(index+1);j<getNsite();j++)
            tmp[j-1]=data[i][j];
          delete[] data[i];
          data[i]=tmp;

        }
      Nsite-=1;
    }
}

void Data::concatenate(Data *ali)
{
  if(getAlphabetName()!=ali->getAlphabetName())
    {
      cerr << "WARNING in Data::concatenate, can not concatenate alignment "<< getFileName() << " and "<< ali->getFileName() <<" which have different alphabet\n";
      cerr.flush();
      exit(0);
    } else {
        for(int i=0;i<getNtaxa();i++)
          {
            int totalsite=getNsite()+ali->getNsite();
            int* totalseq=new int[totalsite];
            for(int j=0;j<getNsite();j++)
              totalseq[j]=data[i][j];
            int aliid=ali->getTaxaIndex(SPname[i]);
            if(aliid!=-1)
              {
                // taxa i in this dataset is found at index aliid in dataset ali
                for(int j=0;j<ali->getNsite();j++)
                  totalseq[getNsite()+j]=(*ali)(aliid,j);
              } else {
                  // taxa not found in ali, fill with gaps
                  for(int j=0;j<ali->getNsite();j++)
                    totalseq[getNsite()+j]=getNstate();
                  cerr << "WARNING in Data::concatenate(), sequence for taxa " <<getTaxaName(i) << " is missing in file " << ali->getFileName() << "\n";
                  cerr.flush();
              }
            delete[] data[i];
            data[i]=totalseq;
          }
        for(int i=0;i<ali->getNtaxa();i++)
          {
            int id=getTaxaIndex(ali->getTaxaName(i));
            if(id==-1)
              {
                // found a new sequence name in ali
                cerr << "WARNING in Data::concatenate(), sequence for taxa " <<ali->getTaxaName(i) << " is missing in file " << getFileName() << "\n";
                cerr.flush();
                int totalsite=getNsite()+ali->getNsite();
                int* totalseq=new int[totalsite];
                for(int j=0;j<getNsite();j++)
                  totalseq[j]=getNstate();
                for(int j=0;j<ali->getNsite();j++)
                  totalseq[getNsite()+j]=(*ali)(i,j);
                int **tmpdata=new int*[getNtaxa()+1];
                for(int k=0;k<getNtaxa();k++)
                  tmpdata[k]=data[k];
                tmpdata[getNtaxa()]=totalseq;
                delete[] data;
                data=tmpdata;
                string *tmpname=new string[getNtaxa()+1];
                for(int k=0;k<getNtaxa();k++)
                  tmpname[k]=getTaxaName(k);
                tmpname[getNtaxa()]=ali->getTaxaName(i);
                delete[] SPname;
                SPname=tmpname;
                Ntaxa+=1;
              }
          }
        Nsite=getNsite()+ali->getNsite();

    }
}

string Data::getAcc(int i)
{
  string res = "";
  if (takeAcc==Yes)
    {
      res = listAcc[i];
    }
  return res;
}
