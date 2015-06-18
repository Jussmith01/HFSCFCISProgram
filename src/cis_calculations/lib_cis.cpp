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

#include "../utils_cpp/lib_includes.h"
//#include "../classes/classes.h"
//#include "../scf_main/scf_classes.h"
#include "lib_cis.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//            Two Electron Integral MO Basis Conversion Class
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/21/2015  ----
                   ----Modified By:               ----
*/
class TEIMO
{
	int NM;//Number of spacial orbitals
	vector<vector<vector<vector<double>>>> vec4dMO;

	public:

        //********************************************//
        //            Class Initialization            //
        //********************************************//
	TEIMO (int N)
	{
		this->NM = N;
		alloc4d(vec4dMO,N);
	};

        //********************************************//
        //                 4D Allocator               //
        //********************************************//
        void alloc4d(vector<vector<vector<vector<double>>>> &arr,int N)
        {
                arr.resize(N);

                for (int i=0;i<N;++i)
                        {arr[i].resize(N);}

                for (int i=0;i<N;++i)
                for (int j=0;j<N;++j)
                        {arr[i][j].resize(N);}

                for (int i=0;i<N;++i)
                for (int j=0;j<N;++j)
                for (int k=0;k<N;++k)
                        {arr[i][j][k].resize(N,0.0);}
        };

        //********************************************//
        //          	  Save Element                //
        //********************************************//
	void SaveElement (int i,int j,int k,int l,double value) 
		{vec4dMO[i][j][k][l] = value;};

        //********************************************//
        //                 Get Element                //
        //********************************************//
	double GetElement (int i,int j,int k,int l)
		{return vec4dMO[i][j][k][l];};

        //********************************************//
        //   	    TEI MO TRansformation	      //
        //********************************************//
	void TEIMOTransform (scf_data *scf,dataOutput *optfile)
	{
		timer trans_timer;
		optfile->ofile << "\nTranforming Two-Electron Integrals to MO Basis..." << endl;
        	trans_timer.set_timer();
        	//optfile->ofile << "Transforming Two Electron Integrals to MO Basis: " << endl;
		optfile->ofile.setf( std::ios::fixed, std::ios::floatfield );

        	int nao = NM;
		int NTEI = NM*NM*NM*NM;
		int step = NTEI / (float)20;
		if(step == 0) {step=1;}
		int cntr=0;

		//cout << "TEST1!!" << endl;

		vector<vector<vector<vector<double>>>> tmp4d;
		alloc4d(tmp4d,NM);

		for(int i=0; i < nao; i++) {
                for(int j=0; j < nao; j++) {
                for(int k=0; k < nao; k++) {
                for(int l=0; l < nao; l++) {
			tmp4d[i][j][k][l] = scf->FEI.get_val(i,j,k,l);
		}}}}

		//cout << "TEST2!!" << endl;
       		for(int i=0; i < nao; i++) {
       		for(int j=0; j < nao; j++) {
       		for(int k=0; k < nao; k++) {
        	for(int l=0; l < nao; l++) {
        	    //cout << " i: " << i << " j: " << j << " k: " << k << " l: " << l;

		    double value = 0.0;		

		    for(int p=0; p < nao; p++) {
		    for(int q=0; q < nao; q++) {
		    for(int r=0; r < nao; r++) {
		    for(int s=0; s < nao; s++) {
		    	value += scf->C[i+p*nao] * scf->C[j+q*nao] * scf->C[k+r*nao] * scf->C[l+s*nao] * tmp4d[p][q][r][s];
		    	//cout << "  -p: " << p << " q: " << q << " r: " << r << " s: " << s << "\n";
		    }}}}

		    vec4dMO[i][j][k][l] = value;

		    if (cntr%step==0) {cout << (cntr/(float)NTEI)*100.0 << endl;}

		    //cout << cntr << " of " << NTEI << endl;
		    ++cntr;
		    //optfile->ofile << setprecision(6) << "   VALUE: " << value << endl;
		}}}}
		cout << "100%" << endl;
        	trans_timer.end_timer();
        	string message = "TEI MO Transform Clock Time: ";
        	trans_timer.print_clock_time(message,optfile->ofile);
	};
};

//________________________________________________________________________//
//      *************************************************************     //
//             	      Hcore MO Basis Conversion Class
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/21/2015  ----
                   ----Modified By:               ----
*/
class HcoreMO
{
	int NM;
	vector<vector<double>> vec2d;

        public:
 
        HcoreMO (int N)
        {
                this->NM = N;
                vec2d.resize(NM);

                for (int i=0;i<NM;++i)
                	{vec2d[i].resize(NM);}
	};

        void SaveElement (int i,int j,double value)
        	{vec2d[i][j] = value;};

        double GetElement (int i,int j)
        	{return vec2d[i][j];};

        void MOTransform(scf_data *scf)
        {
                for(int i=0; i < NM; i++) {
                for(int j=0; j < NM; j++) {
			double value = 0.0;

                	for(int p=0; p < NM; p++) {
                	for(int q=0; q < NM; q++) {
				value += scf->C[i+p*NM] * scf->C[j+q*NM] * scf->Hcore[p+q*NM];
			}}

                        SaveElement(i,j,value);
                }}
        };

	void PrintMatrix(dataOutput *optfile)
	{
		double *MP = new double [NM * NM];
                for(int p=0; p < NM; p++) {
                   for(int q=0; q < NM; q++) {
                          MP[p+q*NM] = GetElement(p,q);
                          //cout << "p: " << p << " q: " << q << " i: " << i << " j: " << j <<  " Val: " << GetElement(p,q) << endl;
                   }}
		optfile->Matrix_printer("MO-Basis Hcore Matrix",MP,NM);
		delete [] MP; MP = NULL;

	}
};

//________________________________________________________________________//
//      *************************************************************     //
//                         CIS Spin-Adapted Class
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/22/2015  ----
                   ----Modified By:               ----
*/
class CIS_tools
{
	int NM;//Number of Molecular orbitals
	int NO;//Number of Occupied molecular orbitals
	int NV;//Number of Virtual molecular orbitals
	int ST;//Number of unique singlets and triplets

	vector<double> Eps;//Contains the orbital energies

	int *min;//Contains the occupied orbital index
	int *vir;//Contains the virtual orbital index

	double *SH;
	double *TH;

        double *ES;
        double *EVS;
        double *ET;
        double *EVT;

	//Transition Dipole Moment Data
	int td_flag;
	int NumE;
	double *TDSM;//Transition Dipole Square Matrix
	double *TDx;//Transition Dipole x Matrix
	double *TDy;//Transition Dipole y Matrix
	double *TDz;//Transition Dipole z Matrix
	double *OSM;//Oscillator Strength Matrix

	string optapp;

	void AllocateHamiltonians() 
	{
		SH = new double [ST*ST];
		TH = new double [ST*ST];
	};
	
	double GetSH (int i,int j) {return SH[j+i*ST];};
	double PutSH (int i,int j,double value) {SH[j+i*ST]=value;};
	double GetTH (int i,int j) {return TH[j+i*ST];};
	double PutTH (int i,int j,double value) {TH[j+i*ST]=value;};
	double GetTDSM (int i) {return TDSM[i];};
	double PutTDSM (int i,double value) {TDSM[i]=value;};
        double GetTDx (int i,int j) {return TDx[j+i*NM];};
        double PutTDx (int i,int j,double value) {TDx[j+i*NM]=value;};
        double GetTDy (int i,int j) {return TDy[j+i*NM];};
        double PutTDy (int i,int j,double value) {TDy[j+i*NM]=value;};
        double GetTDz (int i,int j) {return TDz[j+i*NM];};
        double PutTDz (int i,int j,double value) {TDz[j+i*NM]=value;};
        double GetOSM (int i,int j) {return OSM[j+i*NM];};
        double PutOSM (int i,int j,double value) {OSM[j+i*NM]=value;};


	public:
        //***********************************
        //	  Class Constructor
        //***********************************
	CIS_tools (int N,int td_flag,string optapp)
	{
		this->td_flag = td_flag;
		this->optapp = optapp;
		this->NM = N;
		Eps.resize(N,0.0);
	};

        //***********************************
        // Transition Dipole Square Matrix
        //***********************************
        void MOTransformTDM(scf_data *scf,dataOutput *optfile)
        {
	   if(td_flag==1)
	   {
		TDSM = new double [ST];
		TDx = new double [NM*NM];
		TDy = new double [NM*NM];
		TDz = new double [NM*NM];
		OSM = new double [NM*NM];

		double isq2 = 1/(double)sqrt(2);

                for(int i=0; i < NM; i++) {
                for(int j=0; j < NM; j++) {
                        double valuexp = 0.0,valuexn = 0.0;
                        double valueyp = 0.0,valueyn = 0.0;
                        double valuezp = 0.0,valuezn = 0.0;

                        for(int p=0; p < NM; p++) {
                        for(int q=0; q < NM; q++) {
				double tvalx = scf->C[i+p*NM] * scf->C[j+q*NM] * scf->TDx[p+q*NM];
				double tvaly = scf->C[i+p*NM] * scf->C[j+q*NM] * scf->TDy[p+q*NM];
				double tvalz = scf->C[i+p*NM] * scf->C[j+q*NM] * scf->TDz[p+q*NM];

                                if (tvalx >= 0) {valuexp += tvalx;} else {valuexn += abs(tvalx);}
                                if (tvaly >= 0) {valueyp += tvaly;} else {valueyn += abs(tvaly);}
                                if (tvalz >= 0) {valuezp += tvalz;} else {valuezn += abs(tvalz);}
                        }}

			double valuex = valuexp - valuexn;
			double valuey = valueyp - valueyn;
			double valuez = valuezp - valuezn;

		  	double tmp = ( valuex*valuex + valuey*valuey + valuez*valuez );
			PutTDx(i,j,-valuex);
			PutTDy(i,j,-valuey);
			PutTDz(i,j,-valuez);
			//PutTDSM(i,j,tmp);
			//PutTDSM(i,j,tmp);
                        //cout << "i: " << i << " j: " << j << " VEC: [" << valuex << "," << valuey << "," << valuez << "] VAL: " << tmp << endl;
                }}

		optfile->Matrix_printer("Transition Dipole x Matrix (MO Basis)",TDx,NM);
		optfile->Matrix_printer("Transition Dipole y Matrix (MO Basis)",TDy,NM);
		optfile->Matrix_printer("Transition Dipole z Matrix (MO Basis)",TDz,NM);
		//optfile->Matrix_printer("Transition Dipole Squared Matrix (MO Basis)",TDSM,NM);

		optfile->ofile << endl;

		int cntr = 0;

                ofstream plot;
                stringstream name;
                name << optapp.c_str() << "lplot.dat";
                plot.open (name.str().c_str());

		optfile->ofile << "State Transitions:\n";
		for(int a=0; a < NO; a++) {
                for(int r=0; r < NV; r++) {

			int aidx = min[a];
			int ridx = vir[r];

                	double SE = ES[cntr];
                	double SEeV = SE * 27.21138505;

                	double hc = 1.23984193;

                	double SEnm = (hc/SEeV)*1000;

			double Dx=0;
			double Dy=0;
			double Dz=0;
			double DS=0;

                	plot << SEnm << "   " << 0.6666666 * ES[cntr] * DS << "\n";

			int idx = 0;

			cout << aidx << "->" << ridx << " ";
			for(int b=0; b < NO; b++) {
	                for(int s=0; s < NV; s++) {

				int bidx = min[b];
				int sidx = vir[s];

				double Ci = isq2 * EVS[cntr + idx * ST]; 

				//cout << "  " << bidx << "->" << sidx << " idx: " << cntr + idx * ST << " VAL: " << Ci << endl;
	
				Dx += Ci * GetTDx(bidx,sidx);
				Dy += Ci * GetTDy(bidx,sidx);
				Dz += Ci * GetTDz(bidx,sidx);

				//DS += Dx*Dx + Dy*Dy + Dz*Dz;

				//cout << "  Dx: " << Dx << " Dy: " << Dy << " Dz: " << Dz << " DS: " << DS << endl;

				++idx;
			}}

			Dx = 2.0*Dx;
			Dy = 2.0*Dy;
			Dz = 2.0*Dz;

			DS = Dx*Dx+Dy*Dy+Dz*Dz;
	
			cout << "Dx: " << Dx << " Dy: " << Dy << " Dz: " << Dz << " DS: " << DS << endl;
			TDSM[cntr] = DS;
			
                	optfile->ofile << "  " << min[a] << " -> " << vir[r] 
				<< " ENERGY(AU/eV/nm): " << SE << "/" << SEeV << "/" << SEnm 
				<<  " VEC: [" << -Dx << "," << -Dy << "," << -Dz 
				<<  "] Dip. S.: " << DS << " f: " 
				<< 0.6666666 * ES[cntr] * DS << endl;
           		++cntr;
		}}

		NumE = cntr;

		optfile->ofile << endl;
		plot.close();
	   }
        };

        //***********************************
        // Transition Dipole Square Matrix
        //***********************************
        void TDCalcTest(scf_data *scf,dataOutput *optfile)
        {
           if(td_flag==1)
           {
                double *EVg = new double [NM * NO];
                double *EVe = new double [NM * NO];
                double *Rho = new double [NM * NM];

                double tol = 1.0E-14;

		//ofstream plot;
  		//plot.open ("specplot.dat");

		cout << "Ground State Composed of Orbitals: ";
                for (int o = 0; o < NO; ++o)
                {
                        for (int i = 0; i < NM; ++i)
                        {EVg[o + i * NO] = scf->C[min[o] + i * NM];}
			cout << " g: " << min[o];
                }
		cout << "\n";

		int cntr = 0;

		cout << "State Transition: \n";
		optfile->ofile << "State Transition: \n";

                for(int a=0; a < NO; a++) {
                for(int r=0; r < NV; r++) {

			cout << "Excited State Composed of Orbitals: ";

			for (int o = 0; o < NO; ++o)
                	{
                                int Orb = min[o];
                                if (Orb == min[a])
                                	{Orb = vir[r];cout << " e: " << Orb;}
				else
					{cout << " g: " << Orb;}

                        	for (int i = 0; i < NM; ++i)
                        	{
					EVe[o + i * NO] = scf->C[Orb + i * NM];
				}
                	}

			cout << "\n";

        		for (int i = 0; i < NM; ++i) {
        		for (int j = 0; j < NM; ++j) {
       		                	double pos = 0;
			        	double neg = 0;
	
        	                for (int o = 0; o < NO; ++o)
        	                {
        	                	double val = EVg[o + i * NO] * EVe[o + j * NO];
	
		                        if (val > 0 && val > tol)
	                                	{pos += val;}
	                               	else if (val < 0 && abs(val) > tol)
	                                       	{neg += abs(val);}
	                               	//cout << "C[" << i << "," << o <<  "]*C[" << j << ","  << o << "] pos = " << pos << " neg = " << neg <<  " | ";
	                        }
	
	                        Rho[j + i * NM] = (double) (pos - neg);
		         }}

			 stringstream ss;
                         ss << "DENSITY TEST: " << min[a] << " -> " << vir[r] << " ";
                         optfile->Matrix_printer(ss.str().c_str(),Rho,NM);

                        double valuex = 0.0;
                        double valuey = 0.0;
                        double valuez = 0.0;

                        for(int p=0; p < NM; p++) {
                        for(int q=0; q < NM; q++) {
                                valuex += Rho[p+q*NM] * scf->TDx[q+p*NM];
                                valuey += Rho[p+p*NM] * scf->TDy[q+p*NM];
                                valuez += Rho[p+p*NM] * scf->TDz[q+p*NM];
                        }}

                        double tmp = valuex*valuex + valuey*valuey + valuez*valuez;
                        //cout << "State Transition: " << min[a] << " -> " << vir[r] << " ENERGY: " << ES[cntr] <<  " VEC: [" << -valuex << "," << -valuey << "," << -valuez << "] VAL: " << tmp << " f: " << 0.6666666 * ES[cntr] * tmp << endl;
                        //optfile->ofile << "State Transition: " << min[a] << " -> " << vir[r] << " ENERGY: " << ES[cntr] <<  " VEC: [" << -valuex << "," << -valuey << "," << -valuez << "] VAL: " << tmp << " f: " << 0.6666666 * ES[cntr] * tmp << endl;

			double SE = ES[cntr];
			double SEeV = SE * 27.21138505;
			
			double hc = 1.23984193;

			double SEnm = (hc/SEeV)*1000;

			//plot << SEnm << "   " << 0.6666666 * ES[cntr] * tmp << "\n";

                        cout << " " << min[a] << " -> " << vir[r] << " ENERGY(AU/eV/nm): " << SE << "/" << SEeV << "/" << SEnm <<  " VEC: [" << -valuex << "," << -valuey << "," << -valuez <<  "] Dip. S.: " << tmp << " f: " << 0.6666666 * ES[cntr] * tmp << endl;
                        optfile->ofile << " " << min[a] << " -> " << vir[r] << " ENERGY(AU/eV/nm): " << SE << "/" << SEeV << "/" << SEnm <<  " VEC: [" << -valuex << "," << -valuey << "," << -valuez <<  "] Dip. S.: " << tmp << " f: " << 0.6666666 * ES[cntr] * tmp << endl;
			++cntr;
		}}

  		//plot.close();
		delete [] EVg;  EVg = NULL;
		delete [] EVe;  EVe = NULL;
		delete [] Rho;  Rho = NULL;
           }
        };

        //***********************************
        //Find the Virtual Molecular Orbitals
        //***********************************
	void find_vir(int *min,dataOutput *optfile)
	{
		int vcnt = 0;
		int ocnt = 0;
		vir = new int [NV];
		for (int i=0;i<NM;++i)
		{
			int jcnt = -1;
			for (int j=0;j<NO;++j)
			{
				if (min[j]==i)
				{
					optfile->ofile << "Occupied Orbital[" << ocnt << "] Index: " << min[j] << endl;
					++ocnt;
					jcnt = j;
					break;
				}
			}

			if (jcnt == -1)
			{
				vir[vcnt]=i;++vcnt;
				optfile->ofile << "Virtual Orbital[" << vcnt-1 << "] Index: " << vir[vcnt-1] << endl;
			}
		}
	};

	//***********************************
        //      	FindMax
        //***********************************
        double DetermineMaxVal(vector<double> &array)
        {
		int N = array.size();
		double tmpMax = 0;

		for(int i=0; i < N; i++) 
		{
			double val = array[i];

			if (val > tmpMax)
			{tmpMax=val;}
		}

		return tmpMax;
        };

        //***********************************
        //         Epsilon Function
        //***********************************
        double EpsilonFunction(double std,double f,double nui,double nu)
        {
		double C1 = 1.3062974E+8;
		double C2 = (1.0E+7/(double)3099.6);
		//double C2 = (1.0E+7/(double)600.6);
		//double C3 = 1/(double)600.6;
		double C3 = 1/(double)3099.6;
		double C4 = 1/(double)nu;
		double C5 = 1/(double)nui;

		//cout << "f: " << f << " nui: " << nui;
		//cout << " nu: " << nu << " C1: " << C1 << " C2: " << C2 << " C3: " << C3 << " C4: " << C4 << " C5: " << C5 << endl;

                return C1*(f/C2)*exp(-1.0*((C4-C5)/C3)*((C4-C5)/C3));
        };

        //***********************************
        //	Print Epsilon Graph
        //***********************************
        void PrintEpsilonGraph(dataOutput *optfile)
        {
		if(td_flag==1)
		{
 	                //--------------------------------------------------
 	                //     	        Convert and Save Data
 	                //--------------------------------------------------
                        double hc = 1.23984193;
			vector<double> EE; // nu_i
			vector<double> OS; // f
			int cntr=0;

			for(int a=0; a < NO; a++) {
	                for(int r=0; r < NV; r++) {
                        	int aidx = min[a];
                        	int ridx = vir[r];

                        	double SE = ES[cntr];
                        	double SEeV = SE * 27.21138505;

                        	EE.push_back((hc/SEeV)*1000);
                        	OS.push_back(0.666666*SE*GetTDSM(cntr));

				++cntr;
			}}

			//--------------------------------------------------
                	//                 Setup Graphing    
                        //--------------------------------------------------
                        int gpoints = 5000;
			double std = 0.2;
			double min = 0.0;
                        double max = DetermineMaxVal(EE);
			double dEE = (max + 150.0)/(double)gpoints;

			ofstream plot;
			stringstream name;
			name << optapp.c_str() << "wplot.dat";
                	plot.open (name.str().c_str());

			//--------------------------------------------------
                        //                Begin Graphing    
                        //--------------------------------------------------
			for(int i=0; i < gpoints; i++)
			{
				double sum = 0.0;
				double nu = i*dEE;

				for(int j=0; j < cntr; j++)
				{sum += EpsilonFunction(std,OS[j],EE[j],nu);}

				plot << nu << "   " << sum << endl;
			}
                        //--------------------------------------------------
                        //                Close Graphing    
                        //--------------------------------------------------
			plot.close();
		}
        };

	//**********************************************************
	//Determine Occupied/Virtual Indicies and Test MO Transforms
	//**********************************************************
	void ProduceOrbitalEnergies(scf_data *scf,MemHandler *data_mem,HcoreMO &H,TEIMO &T,dataOutput *optfile)
	{
		double totE = 0;	

        	NO = data_mem->num_elec / 2;

                //--------------------------------------------------
                //     Determine occupied and virtual orbitals
                //--------------------------------------------------
		min = new int [NO];
        	find_mins(min,scf->E,NM,NO);//This function is from the scf functions file

		NV = NM-NO;
		find_vir(min,optfile);//this is a class function

                //--------------------------------------------------
                //Test MO Transforms by calculating orbital energies
                //--------------------------------------------------
        	optfile->ofile << "Orbital Energies: \n\n";
		for (int i=0;i<NM;++i)
		{
			Eps[i] = H.GetElement(i,i);
			//optfile->ofile << "ii = " << i << i << " = " << Eps[i] << endl; 
			for (int m=0;m<NO;++m)
			{
				int Om = min[m];
				Eps[i] += 2.0 * T.GetElement(i,i,Om,Om) - T.GetElement(i,Om,Om,i);
				//optfile->ofile << " iimm = " << i << i << Om << Om << " = " << T.GetElement(i,i,Om,Om) << " immi = "<< i << Om << Om << i << " = " << T.GetElement(i,Om,Om,i) << endl;
			}

			optfile->ofile << "Eps(" << i << "): " << Eps[i] << "\n";
		}

		//----------------------------------------------
        	//Test MO Transforms by calculating total energy
        	//----------------------------------------------
                double Etot = 0.0;
                for (int a=0;a<NO;++a)
                {
			int ai = min[a];
                        Etot += H.GetElement(ai,ai);
                }

                Etot = 2.0 * Etot;

                for (int a=0;a<NO;++a) {
                for (int b=0;b<NO;++b) {
			int ai = min[a];
			int bi = min[b];
                        Etot += 2.0 * T.GetElement(ai,ai,bi,bi) - T.GetElement(ai,bi,bi,ai);	
                }}

		optfile->ofile << "MO Transformed Test - Eelec: " << Etot << endl;
	};

        //**********************************************************
        //	  Build required Spin-Adapted Hamiltonians
        //**********************************************************
	void BuildSpinAdaptedHamiltonians(dataOutput *optfile,TEIMO &T)
	{
		ST = NO*NV;
		optfile->ofile << "Number of Unique Singlets/Triplets: " << ST << endl;
		AllocateHamiltonians();
		
		                //double Escf = data_mem->Eelec;
		int icnt = 0;
		int jcnt = 0;

                for (int a=0;a<NO;++a)
                for (int r=0;r<NV;++r)
                {
			jcnt = 0;
                        for (int b=0;b<NO;++b)
                        for (int s=0;s<NV;++s)
                        {
                                int aidx = min[a];
                                int ridx = vir[r];
                                int bidx = min[b];
                                int sidx = vir[s];

                                optfile->ofile << "a: " << aidx << " r: " << ridx << " b: " << bidx << " s: " << sidx << "| ";

                                double Sval = 0;
                                //Sval += (((iidx == jidx) && (aidx == bidx)) ? Escf : 0);
				if ((ridx == sidx) && (aidx == bidx))
				{Sval += Eps[ridx] - Eps[aidx];}
                                //Sval += ((ridx == sidx) ? Eps[ridx] : 0);
                                //Sval -= ((aidx == bidx) ? Eps[bidx] : 0);
                                Sval += 2*T.GetElement(ridx,aidx,bidx,sidx) - T.GetElement(ridx,sidx,bidx,aidx);

                                PutSH(icnt,jcnt,Sval);

                                double Tval = 0;
                                if ((ridx == sidx) && (aidx == bidx))
                                {Tval += Eps[ridx] - Eps[aidx];}
                                //Tval += ((iidx == jidx) ? Eps[bidx + aidx * NM] : 0);
                                //Tval -= ((aidx == bidx) ? Eps[jidx + iidx * NM] : 0);
                                Tval -= T.GetElement(ridx,sidx,bidx,aidx);

                                PutTH(icnt,jcnt,Tval);

                                //cout << " " <<  << " \n";
				++jcnt;
                        }

                        optfile->ofile << endl;
			++icnt;
                }

		//cout << "TEST1!!\n";
                optfile->Matrix_printer("Singlet Hamiltonian",SH,ST);
                optfile->Matrix_printer("Triplet Hamiltonian",TH,ST);
		//cout << "TEST2!!\n";
	};

        //**********************************************************
        //          Diagonalize Spin-Adapted Hamiltonians
        //**********************************************************
	void DiagHamiltonians (MemHandler *data_mem,dataOutput *optfile)
	{
		if (ST == 1)
		{SetValstoHamiltonians(data_mem,optfile);}
		else
		{DiagonalizeHamiltonians(data_mem,optfile);}
	}

	void DiagonalizeHamiltonians(MemHandler *data_mem,dataOutput *optfile)
        {
                //----------------------------------------------
                //      Diagonalize the Singlet Hamiltonian
                //----------------------------------------------
		ES = new double [ST];
		EVS = new double [ST*ST];

               //Find eigenvectors and eigenvalues
                optfile->ofile << "\nRunning Jacobi Algorithm for Singlet Hamiltonian... \n";

                double *convergence_data;
                convergence_data = new double [3];

                convergence_data[0] = data_mem->ipt_parms.j_conv;
                convergence_data[1] = 6.0 * ST * ST;
                Jacobi_Algorithm(EVS,ES,SH,ST,convergence_data);//Diagonlize H CIS and get eigenvectors and values
                optfile->ofile << "Jacobi Iterations: " << int_to_str(convergence_data[2]) <<" of " << int_to_str(convergence_data[1]) <<  "\n";

                optfile->Matrix_printer("Singlet CIS Eigenvectors",EVS,ST);
                optfile->Array_printer("Singlet CIS Eigenvalues(Hartree)",ES,ST,4);
                optfile->Array_printer_harTOeV("Singlet CIS Eigenvalues(eV)",ES,ST,4);

                //----------------------------------------------
                //     Diagonalize the Triplet Hamiltonian
                //----------------------------------------------
		ET = new double [ST];
                EVT = new double [ST*ST];

                //Find eigenvectors and eigenvalues
                optfile->ofile << "\nRunning Jacobi Algorithm for Triplet Hamiltonian... \n";

                convergence_data[0] = data_mem->ipt_parms.j_conv;
                convergence_data[1] = 6.0 * ST * ST;
                Jacobi_Algorithm(EVT,ET,TH,ST,convergence_data);//Diagonlize H CIS and get eigenvectors and values
                optfile->ofile << "Jacobi Iterations: " << int_to_str(convergence_data[2]) <<" of " << int_to_str(convergence_data[1]) <<  "\n";
                delete [] convergence_data; convergence_data = NULL;

                optfile->Matrix_printer("Triplet CIS Eigenvectors",EVT,ST);
                optfile->Array_printer("Triplet CIS Eigenvalues(Hartree)",ET,ST,4);
                optfile->Array_printer_harTOeV("Triplet CIS Eigenvalues(eV)",ET,ST,4);
	};

        void SetValstoHamiltonians(MemHandler *data_mem,dataOutput *optfile)
        {
                //----------------------------------------------
                //      Diagonalize the Singlet Hamiltonian
                //----------------------------------------------
                ES = new double [ST];
                EVS = new double [ST*ST];

                optfile->ofile << "\nSingle Element Singlets Matrix... \n";

		ES[0] = SH[0];
		EVS[0] = 1.0;

                optfile->Matrix_printer("Singlet CIS Eigenvectors",EVS,ST);
                optfile->Array_printer("Singlet CIS Eigenvalues(Hartree)",ES,ST,4);
                optfile->Array_printer_harTOeV("Singlet CIS Eigenvalues(eV)",ES,ST,4);

                //----------------------------------------------
                //     Diagonalize the Triplet Hamiltonian
                //----------------------------------------------
                ET = new double [ST];
                EVT = new double [ST*ST];

                //Find eigenvectors and eigenvalues
                optfile->ofile << "\nSingle Element Triplets Matrix... \n";

                ET[0] = TH[0];
                EVT[0] = 1.0;

                optfile->Matrix_printer("Triplet CIS Eigenvectors",EVT,ST);
                optfile->Array_printer("Triplet CIS Eigenvalues(Hartree)",ET,ST,4);
                optfile->Array_printer_harTOeV("Triplet CIS Eigenvalues(eV)",ET,ST,4);
        };

        //***********************************
        //      Free Allocated Memory
        //***********************************
	void freeMem ()
	{
        	if(td_flag==1)
        	{
			delete [] TDSM; TDSM = NULL;
        		delete [] TDx; TDx = NULL;
        		delete [] TDy; TDy = NULL;
        		delete [] TDz; TDz = NULL;
        		delete [] OSM; OSM = NULL;
		}

                delete [] ES; ES = NULL;
                delete [] EVS; EVS = NULL;
                delete [] ET; ET = NULL;
                delete [] EVT; EVT = NULL;

		delete [] SH; SH = NULL;
		delete [] TH; TH = NULL;
		delete [] vir; vir = NULL;
		delete [] min; min = NULL;
	};
};

//________________________________________________________________________//
//      *************************************************************     //
//                             CIS MAIN CALC
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/21/2015   ----
                   ----Modified By:               ----
*/

extern void cis_energy_calc(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	if (data_mem->ipt_parms.cis_flag==1)
	{
		timer cis_timer;
	        cis_timer.set_timer();

		optfile->ofile.setf( std::ios::fixed, std::ios::floatfield );
        	optfile->ofile << setprecision(6) << "|------------------------------------------------------|\n";
		optfile->ofile << "|-----------------Begin CIS Calculation----------------|\n";
        	optfile->ofile << "|------------------------------------------------------|\n";

		cout << "CIS Calculation Beginning...\n";

		int N = scf->orbitals;

		int Nocc = data_mem->num_elec;//Number of occupied orbitals
		int Nvir = 2*N-Nocc;//Number of virtual orbitals
		int NSingles = Nvir * Nocc;//Number of virtual orbitals
		optfile->ofile << "Number of occupied spin orbitals: " << Nocc << " Number of virtual: " << Nvir << endl;

                //********************************************//
                //   	   	Buld MO Form of TEI	      //
                //********************************************//
		cout << "Transforming TEI...\n";

		TEIMO TEI(N);//Allocate Class Memory with N*2 spin orbitals
		TEI.TEIMOTransform(scf,optfile);//Build TEI spin basis 4d matrix

                //********************************************//
                //             Buld MO Form of Hcore          //
                //********************************************//
		optfile->ofile << "\nTranforming Hcore Integrals to MO Basis..." << endl;
		cout << "Transforming Hcore...\n";
	
		HcoreMO Hcore(N);//Allocate Class Memory
		Hcore.MOTransform(scf);
		Hcore.PrintMatrix(optfile);

                optfile->ofile << "\nCalculating Orbital Energies..." << endl;
		cout << "Calculating Orbital Energies...\n";

                //******************************************************//
                //    Build and Diagonalize CIS Hamiltonian Matrix      //
                //******************************************************//
		cout << "Cuilding and Diagonalizing the CIS Hamiltonian...\n";
                CIS_tools CIS(N,scf->td_flag,data_mem->ipt_parms.dataout);//Allocate Class Memory
                CIS.ProduceOrbitalEnergies(scf,data_mem,Hcore,TEI,optfile);
		CIS.BuildSpinAdaptedHamiltonians(optfile,TEI);
		CIS.DiagHamiltonians(data_mem,optfile);

		//******************************************************//
                //  Transform the Transition Dipole Matrix to MO Space  //
		//	Also Calculates the Sqare values into TDSM	//
                //******************************************************//
		cout << "Calculate Transition Dipoles...\n";
		CIS.MOTransformTDM(scf,optfile);
		CIS.PrintEpsilonGraph(optfile);
		//CIS.TDCalcTest(scf,optfile);

		//********************************************//
                //             Free the Memory 	              //
                //********************************************//
		CIS.freeMem();

		cis_timer.end_timer();
        	string message = "CIS Clock Time: ";
        	cis_timer.print_clock_time(message,optfile->ofile);

		cout << "End of CIS Calculation\n";
                optfile->ofile << "|------------------------------------------------------|\n";
                optfile->ofile << "|------------------End CIS Calculation-----------------|\n";
                optfile->ofile << "|------------------------------------------------------|\n";
	} else {
        	optfile->ofile << "CIS Calculation Disabled... moving on.\n\n";
	}
}
