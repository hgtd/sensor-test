 	TString fname = "data/0716/";
 	TString fname = "data/0716/";
void pltdoping()
{ 
// #include <iostream> 

// 	TString fname = "data/0716/HPK-W28-LGAD-L5P6-CV";
// 	TString fname = "data/0716/HPK-W33-LGAD-L5P6-CV";
// 	TString fname = "data/0716/HPK-W33-PIN-L5P9-CV";
// 	TString fname = "data/0716/HPK-W37-LGAD-L5P6-CV";
// 	TString fname = "data/0716/HPK-W37-PIN-L5P9-CV";
// 	TString fname = "data/0713/HPK-W43-LGAD-L5P6-CV";
// 	TString fname = "data/0713/HPK-W43-PIN-L5P9-CV";
// 	TString fname = "data/0713/HPK-W43-PIN-L5P9-CV";
// 	TString fname = "data/0728/HPK-W37-LGAD-L5P5-CV";
// 	TString fname = "data/0728/HPK-W33-LGAD-L5P5-CV";
 	TString fname = "data/0728/HPK-W28-LGAD-L5P5-CV";
	TGraph *g_cap = new TGraph( Form("%s.csv",fname.Data()), "%lf,%*lf,%lf");
//	TGraph *g_cap = new TGraph("28lgad.csv", "%lf,%*lf,%lf");
    int nvol = g_cap->GetN();
	std::cout<<nvol<<std::endl;
    double *vol = g_cap->GetX(); // [V]
    double *cap = g_cap->GetY(); // [pF]
    double step = fabs(vol[1]-vol[0]);

    double A = 0.13*0.13; // [cm^2] // HPK
    double esi = 11.7;
    double e0 = 8.854E-14; // [F/cm]
    double q = 1.602E-19; // [C]

    double *invcap2 = new double[nvol]; // 1/C^2
    double *invcap2_dv = new double[nvol-1]; // d(1/C^2)/dV
    double *doppf = new double[nvol-1];
    double *depth = new double[nvol-1];

    // double *invcap2 = new double[nvol]; // 1/C^2
    // double *invcap2_dv = new double[nvol]; // d(1/C^2)/dV
    // double *doppf = new double[nvol];
    // double *depth = new double[nvol];


    for(int i=0; i<nvol; ++i){
        vol[i] = -vol[i]; // convert to positive
        cap[i] = cap[i]*1E-12; // [pF] to [F]
        invcap2[i] = 1.0/cap[i]/cap[i]; // [1/F^2]
        g_cap->SetPoint(i, vol[i], cap[i]);
    }

    ofstream outFile;
    outFile.open(Form("%s_doping.csv",fname.Data()),  ios::out);
	// TGraph *g_cap = new TGraph( Form("%s.csv",fname.Data()), "%lf,%*lf,%lf");
    outFile << "Voltage" << "," << "depth" << "," << "dop" << endl;

    for(int i=0; i<nvol-1; ++i){ // calculate derivatives and doping profile
        invcap2_dv[i] = (invcap2[i+1]-invcap2[i])/step; // forward difference
        doppf[i] = 2.0/q/esi/e0/A/A/invcap2_dv[i]; // [cm^{-3}]
        if(doppf[i] < 0.){doppf[i] = 0.1;}
        depth[i] = A*esi*e0*(1.0/cap[i])*1E+4; // [cm] to [um]
        outFile << vol[i] << "," << depth[i] << "," << doppf[i] << endl;
//        cout<<vol[i]<<"V "<<cap[i]<<"F "<<invcap2[i]<<"F^{-2} "<<depth[i]<<"um "<<doppf[i]<<"cm^{-3} "<<endl; // for debugging
    }

    outFile.close();
}

   
