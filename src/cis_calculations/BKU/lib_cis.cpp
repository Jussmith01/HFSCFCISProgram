#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <vector>

//#include "../utils_cpp/lib_includes.h"
//#include "../classes/classes.h"
//#include "../scf_main/scf_classes.h"
#include "lib_cis.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//            Two Electron Integral Spin Basis Conversion Class
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 3/30/2015  ----
                   ----Modified By:               ----
*/
class TEIsb
{
	int N;//Number of spin orbitals
	vector<vector<vector<vector<double>>>> vec4d;

	public:

        //********************************************//
        //            Class Initialization            //
        //********************************************//
	TEIsb (int N)
	{
		this->N = N;
		vec4d.resize(N);

		for (int i=0;i<N;++i)
		{
			vec4d[i].resize(N);
		}

                for (int i=0;i<N;++i)
                	for (int j=0;j<N;++j)
			{
                        	vec4d[i][j].resize(N);
			}

		for (int i=0;i<N;++i)
                        for (int j=0;j<N;++j)
				for (int k=0;k<N;++k)
                        	{
                               		vec4d[i][j][k].resize(N);
                        	}
	};

        //********************************************//
        //          	  Save Element                //
        //********************************************//
	void SaveElement (int i,int j,int k,int l,double value) 
	{
		vec4d[i][j][k][l] = value;
	};

        //********************************************//
        //                 Get Element                //
        //********************************************//
	double GetElement (int i,int j,int k,int l)
	{
		return vec4d[i][j][k][l];	
	};

        //********************************************//
        //   Build Two Electron Integral Spin Basis   //
        //********************************************//
	void BuildTEIsb(scf_data *scf)
	{
		/* Translate integrals to spin-orbital basis */
		cout << "Number of TEI: " << N*N*N*N << " Num Spatial Orbitals: " << scf->orbitals << " NumSpin: " << N <<  endl;
  		for(int p=0; p < N; p++)
    		   for(int q=0; q < N; q++)
      		      for(int r=0; r < N; r++)
        		for(int s=0; s < N; s++) 
			{
        		  int i = p/2;
			  int j = q/2;
        		  int k = r/2;
			  int l = s/2;
			  //cout << "TEST1\n";
			  //cout << "i: " << i << " j: " << j << " k: " << k << " l: " << l << endl;
        		  double value1 = scf->FEI.get_val(p,q,r,s) * (p%2 == q%2) * (r%2 == s%2);
			  if(value1 < 1.0E-14) {value1 == 0.0f;}
			  //cout << "TEST2\n";
        		  //double value2 = scf->FEI.get_val(i,l,j,k) * (p%2 == s%2) * (q%2 == r%2);
			  //if(value2 < 1.0E-14) {value2 == 0.0f;}
			  //cout << "TEST3\n";
        		  SaveElement(p,q,r,s,value1);
			  //cout << "p: " << p << " q: " << q << " r: " << r << " s: " << s << " Val: " << GetElement(p,q,r,s) << endl;
        		}
	};
};

//________________________________________________________________________//
//      *************************************************************     //
//             	      Hcore Spin Basis Conversion Class
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 3/30/2015  ----
                   ----Modified By:               ----
*/
class Hcoresb
{
	int N;
	vector<vector<double>> vec2d;

        public:
 
        Hcoresb (int N)
        {
                this->N = N;
                vec2d.resize(N);

                for (int i=0;i<N;++i)
                {
                        vec2d[i].resize(N);
                }
	};

        void SaveElement (int i,int j,double value)
        {
                vec2d[i][j] = value;
        };

        double GetElement (int i,int j)
        {
                return vec2d[i][j];
        };

        void BuildHcoresb(scf_data *scf)
        {
                  /* Translate integrals to spin-orbital basis */
                for(int p=0; p < N; p++)
                   for(int q=0; q < N; q++)
                   {
                          int i = p/2;
                          int j = q/2;
                          double value1 = scf->Hcore[j + i * scf->orbitals] * ((p == q) ? 1.0f : 0.0f);
			  if(value1 < 1.0E-14) {value1 == 0.0f;}
                          SaveElement(p,q,value1);
                          //cout << "p: " << p << " q: " << q << " i: " << i << " j: " << j <<  " Val: " << GetElement(p,q) << endl;
                   }
        };

	void PrintMatrix(dataOutput *optfile)
	{
		double *MP = new double [N * N];
                for(int p=0; p < N; p++)
                   for(int q=0; q < N; q++)
                   {
                          MP[q+p*N] = GetElement(p,q);
                          //cout << "p: " << p << " q: " << q << " i: " << i << " j: " << j <<  " Val: " << GetElement(p,q) << endl;
                   }
		optfile->Matrix_printer("Spin-Orbital Hcore Matrix",MP,N);
		delete [] MP; MP = NULL;

	}
};

//________________________________________________________________________//
//      *************************************************************     //
//                             CIS MAIN CALC
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 3/26/2015  ----
                   ----Modified By:               ----
*/

extern void cis_energy_calc(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	if (data_mem->ipt_parms.cis_flag==1)
	{
        	optfile->ofile << "|------------------------------------------------------|\n";
		optfile->ofile << "|-----------------Begin CIS Calculation----------------|\n";
        	optfile->ofile << "|------------------------------------------------------|\n";

		int N = scf->orbitals;
		optfile->Matrix_printer("SCF Fock Matrix: ",scf->F,N);

		int Nocc = data_mem->num_elec;//Number of occupied orbitals
		int Nvir = 2*N-Nocc;//Number of virtual orbitals
		int NSingles = Nvir * Nocc;//Number of virtual orbitals
		optfile->ofile << "Number of occupied spin orbitals: " << Nocc << " Number of virtual: " << Nvir << endl;

                //********************************************//
                //   	  Buld Spin Orbital Form of TEI       //
                //********************************************//
		optfile->ofile << "\nBuilding Two-Electron Spin-Orbital form of TEI..." << endl;

		TEIsb TEI(N*2);//Allocate Class Memory with N*2 spin orbitals
		TEI.BuildTEIsb(scf);//Build TEI spin basis 4d matrix

                //********************************************//
                //       Buld Spin Orbital Form of Hcore      //
                //********************************************//
		optfile->ofile << "\nBuilding Hcore Spin-Orbital form of TEI..." << endl;
	
		Hcoresb Hcore(N*2);//Allocate Class Memory
		Hcore.BuildHcoresb(scf);//Build TEI spin basis 2d matrix
		Hcore.PrintMatrix(optfile);

		optfile->ofile << "\nTotal number of possible single excitations: " << NSingles << endl;
		optfile->ofile << "\nBuilding CIS Hamitonian Matrix...\n";

		//Note, Excluding the ground state since it does not couple with the single excited states.
		double *Fcis = new double [N*2 * N*2];

	        int NO = data_mem->num_elec / 2;

	        //********************************************//
	        //   Find Minimum Energy Orbital Occupancy    //
	        //********************************************//
	        int *min = new int [NO];
		find_mins(min,scf->E,N,NO);

		vector<int> minsb;
		for (int i=0;i<N*2;++i)
		{
			int idx = i/2;
			minsb.push_back(min[idx]*2 + i%2);
		}
		delete [] min; min = NULL;
		

                //********************************************//
                //      Calculate SpinOrbit Fock Matrix       //
                //********************************************//
		for (int p=0;p<N*2;++p)
		{	
			for (int q=0;q<N*2;++q)
			{
				double Fpq = Hcore.GetElement(p,q);
				cout << "p: " << p << " q: " << q;
				for (int j=0;j<Nocc;++j)
				{
					int m = minsb[j];
					cout << " m[" << j << "]: " << m;
					Fpq += TEI.GetElement(p,m,q,m);
				}
				cout << endl;
				Fcis[q + p * N * 2] = Fpq;
			}
		}
		optfile->Matrix_printer("Spin-Orbital Fock Matrix",Fcis,N*2);

                //********************************************//
                //             Obtain Unocc. Orbs.            //
                //********************************************//
		double *Hcis = new double [NSingles * NSingles];

		vector<int> maxsb;
		for (int i=0;i<N*2;++i)
		{
			bool idxfound = false;
			for (int j=0;j<Nocc;++j)
			{
				if (i == minsb[j])
					{idxfound=true;}
			}

			if (!idxfound)
			{maxsb.push_back(i);cout << "Unocc. Orb. Found: " << i << endl ;}
		}

                //********************************************//
                //        Build CIS Hamiltonian Matrix        //
                //********************************************//
		//cout << "MinSB Size: " << (int)minsb.size() << " MaxSB Size: " << (int)maxsb.size() << endl;
		int Nuocc = (int)maxsb.size();
		double Escf = data_mem->Eelec;
		for (int i=0;i<Nocc;++i)
                for (int a=0;a<Nuocc;++a)
                {
		   	for (int j=0;j<Nocc;++j)
		   	for (int b=0;b<Nuocc;++b)
		   	{
				int iidx = minsb[i];
				int aidx = maxsb[a];
				int jidx = minsb[j];
				int bidx = maxsb[b];

				cout << "i: " << iidx << " a: " << aidx << " j: " << jidx << " b: " << bidx; 

				double val = 0;
				//val += (((iidx == jidx) && (aidx == bidx)) ? Escf : 0);
				val += ((iidx == jidx) ? Fcis[bidx + aidx * N * 2] : 0);
				val -= ((aidx == bidx) ? Fcis[jidx + iidx * N * 2] : 0);
				val += TEI.GetElement(aidx,jidx,iidx,bidx);

                           	Hcis[(j*Nuocc+b) + (i*Nuocc+a) * NSingles] = val;

				cout << " " << (j*Nuocc+b) + (i*Nuocc+a) * NSingles << " \n";
		   	}
			//cout << endl;
                }

		optfile->Matrix_printer("CIS Hamiltonian",Hcis,NSingles);

                //********************************************//
                //    Diagonalize CIS Hamiltonian Matrix      //
                //********************************************//
		cout << "TEST6\n";
		double *Evec = new double [NSingles*NSingles];
		double *Eval = new double [NSingles];

                //Find eigenvectors and eigenvalues
                optfile->ofile << "Running Jacobi Algorithm... \n";

                double *convergence_data;
                convergence_data = new double [3];

		cout << "TEST7\n";
                convergence_data[0] = data_mem->ipt_parms.j_conv;
                convergence_data[1] = 8.0 * NSingles * NSingles;
                Jacobi_Algorithm(Evec,Eval,Hcis,NSingles,convergence_data);//Diagonlize H CIS and get eigenvectors and values
                optfile->ofile << "Jacobi Iterations: " << int_to_str(convergence_data[2]) <<" of " << int_to_str(convergence_data[1]) <<  "\n";
                delete [] convergence_data; convergence_data = NULL;

		optfile->Matrix_printer("CIS Eigen Vectors",Evec,NSingles);
                optfile->Array_printer("CIS Eigen Values",Eval,NSingles,4);

		//********************************************//
                //    		   Clean Up                   //
                //********************************************//
		cout << "TEST8\n";
		delete [] Evec; Evec = NULL;
		delete [] Eval; Eval = NULL;
		delete [] Hcis; Hcis = NULL;
		delete [] Fcis; Fcis = NULL;

		cout << "TEST9\n";
                optfile->ofile << "|------------------------------------------------------|\n";
                optfile->ofile << "|------------------End CIS Calculation-----------------|\n";
                optfile->ofile << "|------------------------------------------------------|\n";
	} else {
        	optfile->ofile << "CIS Calculation Disabled... moving on.\n\n";
	}
}
