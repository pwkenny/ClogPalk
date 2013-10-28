//  ClogPalk  Version 1.0
//  Peter W. Kenny, IQSC-USP, 20-May-2013
//  This source code is provided as supporting information for:
//  PW Kenny, CA Montanari, IM Prokopcyk (2013) ClogPalk: A method for predicting alkane/water
//  partition coefficient. JCAMD 27:*-*  (If citing please consult journal for page numbers) 
  

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include "openeye.h"
#include "oechem.h"
#include "oespicoli.h"


using namespace OEChem;
using namespace OESystem;
using namespace std;
using namespace OESpicoli;

const char *InterfaceData = 
"!BRIEF UsingOEInterfaceHelp [-o] <output> [-i] <input> [-p] <parameters> [-v] <vb> [-d] <debug> \n"
"!PARAMETER -input\n"
"  !ALIAS -i\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF Input file\n"
"  !DETAIL\n"
"         Input file of molecules\n"
"!END\n"
"!PARAMETER -output\n"
"  !ALIAS -o\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF Output file\n"
"  !DETAIL\n"
"         Output file of predicted ClogPalk values\n"
"!END\n"
"!PARAMETER -param\n"
"  !ALIAS -p\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF SMARTS\n"
"  !DETAIL\n"
"         Parameters\n"
"!END\n"
"!PARAMETER -vectorbind\n"
"  !ALIAS -v\n"
"  !TYPE string\n"
"  !REQUIRED false\n"
"  !BRIEF Vector bindings\n"
"  !DETAIL\n"
"         Vector bindings (optional)\n"
"!END\n"
"!PARAMETER -debug\n"
"  !ALIAS -d\n"
"  !TYPE bool\n"
"  !REQUIRED false\n"
"  !BRIEF Debug mode\n"
"  !DETAIL\n"
"         Specfied debugging output\n"
"!END\n" ;



int main(int argc , char** argv)

{

 OEInterface itf(InterfaceData, argc, argv);

 OEGraphMol mol;
 OESurface surf;
 OESubSearch sss;
 vector<OESubSearch> ss_fgrp , ss_interact;

 double dlogp , dlogp_interc , dlogp_slope , int_dlogp , ref_eqn_interc , ref_eqn_slope , logp , msa ;

 string vbname , vbdef , line , smarts , dummy , molname ;
 vector< pair<string, string> > definitions;
 vector<string> fgrp_smt , interact_smt;
 vector<double> fgrp_dlogp_interc , fgrp_dlogp_slope , interact_dlogp;
 unsigned int i , j ;
 bool debug;

 oemolistream infil;
 ifstream parfil;
 ifstream vbfil;
 ofstream outfil;



 //Kommand line stuff
 infil.open(itf.Get<std::string>("-i"));
 outfil.open(itf.Get<std::string>("-o").c_str());
 parfil.open(itf.Get<std::string>("-p").c_str());

 if (itf.Has<bool>("-d")) debug = itf.Get<bool>("-d");

 if (itf.Has<string>("-v"))
   {
     vbfil.open(itf.Get<std::string>("-v").c_str());
     if (!vbfil) OEThrow.Error("Unable to open vector binding file\n",itf.Get<std::string>("-v").c_str());
   }

 if (!infil) OEThrow.Error("Unable to open input file\n",itf.Get<std::string>("-i").c_str());
 if (!parfil) OEThrow.Error("Unable to open SMARTS file\n",itf.Get<std::string>("-s").c_str());

  // Read vektor bindings

 if (vbfil)
    {
      while(getline(vbfil,line))
        {
          if (!vbfil.eof())
            {  
               if ( line.find("#") != 0) 
                 {
                    istringstream iss(line);
                    iss >> vbname >> vbdef; 
                    definitions.push_back(pair<string,string>( vbname, vbdef ));
                 }
            }
        }
    }

//  Parse parameter file
 
 while(getline(parfil,line))
   {
     if (!parfil.eof())
        {  
           if ( line.find("#") != 0) 
              {
		if ( line.find("f") == 0 || line.find("F") == 0 )
		  { 
                    istringstream iss(line);
                    if ( iss >> dummy >> smarts >> dlogp_slope >> dlogp_interc )
		      {
                        fgrp_smt.push_back(smarts);
                        fgrp_dlogp_slope.push_back(dlogp_slope);
                        fgrp_dlogp_interc.push_back(dlogp_interc);
                      }
		  }
		else if ( line.find("i") == 0 || line.find("I") == 0 )
		  { 
                    istringstream iss(line);
                    if ( iss >> dummy >> smarts >> int_dlogp )
		      {
                        interact_smt.push_back(smarts);
                        interact_dlogp.push_back(int_dlogp);
                      }
		  }
                else if ( line.find("r") == 0 || line.find("R") == 0 )
		  { 
                    istringstream iss(line);
                    if  (!( iss >> dummy >> ref_eqn_slope >> ref_eqn_interc ))
		      {
                        exit(0);
                      }
		  }
                else
                  {
		    cout << "Failed to parse line in parameter file: " << line << endl;
                  }

               }
         }
   }

 if (debug)
   {
     cout << "intercept: "  << ref_eqn_interc << " slope: " <<  ref_eqn_slope << endl;
     cout << "fgroup array size: " << fgrp_smt.size()  << endl;
     for ( i = 0 ; i < fgrp_smt.size() ; i++ )
        cout << "SMARTS: " << fgrp_smt[i] << " " << fgrp_dlogp_interc[i] << " " << fgrp_dlogp_slope[i]  << endl;
        cout << "interaction array size: " << interact_smt.size()  << endl;
     for ( i = 0 ; i < interact_smt.size() ; i++ )
     cout << "SMARTS2: " << interact_smt[i] << " " << interact_dlogp[i] << endl;
   }

 // Creat substructural search objects from functional group and interaction SMARTS

 for ( i = 0 ; i < fgrp_smt.size() ; i++ )
   {
     OESmartsLexReplace(fgrp_smt[i], definitions);
     if (!sss.Init(fgrp_smt[i].c_str()))
       {
         OEThrow.Error("Unable to parse %s",fgrp_smt[i].c_str());
       }
     else
       {
         ss_fgrp.push_back(sss);
        }
    }

 for ( i = 0 ; i < interact_smt.size() ; i++ )
   {
     OESmartsLexReplace(interact_smt[i], definitions);
     if (!sss.Init(interact_smt[i].c_str()))
       {
         OEThrow.Error("Unable to parse: %s",interact_smt[i].c_str());
       }
     else
       {
         ss_interact.push_back(sss);
        }
    }


 //  Write header for output file
 outfil << "Name                                ClogPalk     MSA/Ang**2      logPcorr\n";   
 

// Read input molecules

 while (OEReadMolecule(infil, mol))
   {
     // Calculate logP for the saturated hydrocarbon reference molecule
     OEAssignBondiVdWRadii(mol);
     OEMakeMolecularSurface(surf, mol);
     msa = OESurfaceArea(surf);
     logp = ref_eqn_interc + ref_eqn_slope*msa;
     dlogp = 0;
     molname = mol.GetTitle();
     if ( debug)
          cout << "Name: " << mol.GetTitle()<< " Area/Ang**2: " << OESurfaceArea(surf) << "\n";

     //  Initialisation
     double *atm_dlogp = new double[mol.GetMaxAtomIdx()];
     for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom ; ++atom ) atm_dlogp[atom->GetIdx()] = 0;
     vector<char> fgrp_flag( ss_fgrp.size() , 'n');

     for ( i = 0 ; i < ss_fgrp.size() ; i++) 
	     {
               
	       for (OEIter<OEMatchBase> match = ss_fgrp[i].Match(mol); match ; ++match )
                 {
                   OEIter<OEMatchPair<OEAtomBase> > apr;
                   j = 0;
                   for (apr = match->GetAtoms();apr;++apr,++j)
		       {
                         fgrp_flag[i] = 'y';
                         if ( j == 0 )
                            atm_dlogp[apr->target->GetIdx()] = fgrp_dlogp_slope[i];
                        }
                  }
	     }

     // Apply atomic logP increments
     for ( i = 0 ; i < ss_fgrp.size() ; i++) 
             {
               if ( fgrp_flag[i] == 'y' )
                  dlogp = dlogp + fgrp_dlogp_interc[i]; 
             }   

     for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) 
          dlogp = dlogp + atm_dlogp[atom->GetIdx()];

     //   Check
     if (debug) 
         cout << "Checking the output " << molname << " " <<  logp << " " << dlogp << "\n";

     // Apply interaction logP increments
         
     for ( i = 0 ; i < ss_interact.size() ; i++) 
	     {
	       for (OEIter<OEMatchBase> match = ss_interact[i].Match(mol); match ; ++match )
                 {
		   // OEIter<OEMatchPair<OEAtomBase> > apr;
                   dlogp = dlogp + interact_dlogp[i];
                   if (debug)
                       cout << "Interact " << i << " " <<  interact_dlogp[i] << endl;
                  }
	     }
     
     // Output predicted logP 
     logp = logp - dlogp;

     outfil.setf(ios::left); 
     outfil << setw(20) << molname;
     outfil.unsetf(ios::left);
     outfil.setf(ios::right);
     outfil.setf(ios::fixed, ios::floatfield );
     outfil.precision(3);
 
     outfil << setw(20)  << logp; 
     outfil << setw(20)  << msa;
     outfil << setw(20) << dlogp << "\n";
     outfil.unsetf(ios::right);

     // Check

     if ( debug )
       {     
          for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom)
          cout <<  atom->GetIdx() <<  " " << atom->GetAtomicNum() <<  " " << atm_dlogp[atom->GetIdx()] << 
          " " << atom->GetRadius() <<std::endl ;
       }
    mol.Clear();
   }
 infil.close(); 

}







