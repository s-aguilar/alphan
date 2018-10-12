{
  gROOT->Reset();

  // Get input file
  char input_filename[128];

  cout << "Enter SRIM input file: ";
  cin >> input_filename;

  //ifstream *input_file = new ifstream(input_filename);
  ifstream *input_file = new ifstream("10C_HeCO2_400_Torr.txt");

  if(input_file->is_open()){
    cout << "File opened successfully!" << endl;
  }
  else{
    cout << "Could not open file" << endl;
    exit(-1);
  }

  // defines region of tabulated numbers
  int first_line;
  int last_line;

  cout << "Enter first line of data: ";
  cin >> first_line;

  cout << "Enter last line of data: ";
  cin >> last_line;

  int n_data = last_line - first_line + 1;

  // parameters read in from data
  double *energy = new double[n_data];
  double *stopping_elec = new double[n_data];
  double *stopping_nuclear = new double[n_data];

  // calculated parameters
  double *x = new double[n_data];
  double *stopping_total = new double[n_data];
  double *dxdE_total = new double[n_data];
  double *Energy_x = new double[n_data]; // E(x)

  // skip to first line of data
  char c_dummy[100];
  for(int i=0; i<(first_line - 1); i++){
    input_file->getline(c_dummy, 100);
  }

  // to convert to same energy units
  for(int i=0; i<n_data; i++){
    // for energy units
    char e_order[5];
    string keV("keV");
    string MeV("MeV");

    double  f_dummy; // double variable

    // read in values
    *input_file >> energy[i] >> e_order >> dec >> stopping_elec[i]
		>> stopping_nuclear[i] >> f_dummy >> c_dummy
		>> f_dummy >> c_dummy >> f_dummy >> c_dummy;

    // convert energies to MeV
    if(keV.compare(e_order)==0){
      energy[i] = energy[i]/1000.0;
    }
    else if(MeV.compare(e_order)==0){
      energy[i] = energy[i];
    }

    // calculate total dE/dx
    stopping_total[i] = stopping_elec[i] + stopping_nuclear[i];

    // calculate dx/dE
    dxdE_total[i] = 1.0/stopping_total[i];
  }


  // calculate x
  x[n_data-1] = 0.0; // starting point; remember indexing starts at 0

  //cout << n_data - 1<< "\t" << x[n_data-1] << endl;


  // "Stair step" of integration: may need to improve this
  for(int i=((n_data-1)-1); i>=0; i--){
    x[i] = x[i+1] + 0.5*(dxdE_total[i+1] + dxdE_total[i])*(energy[i+1] - energy[i]);
  }

  // print results
  for(int i=0; i<n_data; i++){
    cout << i << "\t" << energy[i] << "\t" <<  x[i] << "\t" << dxdE_total[i] << "\t" << stopping_total[i] << endl;
  }

  // plot results
  TCanvas *c1 = new TCanvas("c1");
  TGraph *dEdx_E = new TGraph(n_data, energy, stopping_total);
  dEdx_E->SetLineColor(2);
  dEdx_E->SetLineWidth(2);
  dEdx_E->SetMarkerStyle(8);
  dEdx_E->GetXaxis()->SetTitle("energy (keV)");
  dEdx_E->GetXaxis()->CenterTitle();
  dEdx_E->GetYaxis()->SetTitle("dE/dx (MeV/mm)");
  dEdx_E->GetYaxis()->CenterTitle();
  dEdx_E->Draw("ALP");

  TCanvas *c2 = new TCanvas("c2", "Stopping Power");
  TGraph *dEdx_x = new TGraph(n_data, x, stopping_total);
  dEdx_x->SetLineColor(2);
  dEdx_x->SetLineWidth(2);
  dEdx_x->SetMarkerStyle(8);
  dEdx_x->Draw("AP");
  dEdx_x->GetXaxis()->SetTitle("x (mm)");
  dEdx_x->GetXaxis()->CenterTitle();
  dEdx_x->GetYaxis()->SetTitle("dE/dx (MeV/mm)");
  dEdx_x->GetYaxis()->CenterTitle();

  // define spline knots
  const Int_t n_knots = 13;
  Int_t knots[n_knots] = {93,90,88,83,80,77,71,67,62,53,39,12,0};

  Double_t x_spline_data[n_knots];
  Double_t dEdx_spline_data[n_knots];
  Double_t E_spline_data[n_knots];
  for(Int_t i=0; i<n_knots; i++){
    x_spline_data[i] = x[knots[i]];
    dEdx_spline_data[i] = stopping_total[knots[i]];
    E_spline_data[i] = energy[knots[i]];
    cout << x_spline_data[i] << "\t" << dEdx_spline_data[i] << "\t" << E_spline_data[i] << endl;
  }

  TSpline3 *dEdx_spline = new TSpline3("dEdx_spline",x_spline_data,dEdx_spline_data,n_knots);
  dEdx_spline->SetLineColor(2);
  dEdx_spline->SetLineWidth(2);
  dEdx_spline->Draw("same");


  TCanvas *c3 = new TCanvas("c3");
  TGraph *E_x = new TGraph(n_data, x, energy);
  E_x->SetLineColor(2);
  E_x->SetLineWidth(2);
  E_x->SetMarkerStyle(8);
  E_x->Draw("AP");
  E_x->GetXaxis()->SetTitle("x (mm)");
  E_x->GetXaxis()->CenterTitle();
  E_x->GetYaxis()->SetTitle("E (MeV)");
  E_x->GetYaxis()->CenterTitle();

  TSpline3 *E_spline = new TSpline3("E_spline",x_spline_data,E_spline_data,n_knots);
  E_spline->SetLineColor(2);
  E_spline->SetLineWidth(2);
  E_spline->Draw("same");

  Double_t x_data;
  cout << "Enter x position: ";
  cin >> x_data;
  cout << "Energy at " << x_data << " is " << E_spline->Eval(x_data) << endl;

  const Double_t M1_M2 = 14.0; // sum of m1 + m2 ?

  // Integral of fit, E(x)
  //TF1 *E = new TF1("E", "[3]-([0]*x + [1]/2.0*pow(x,2) + [2]/3.0*pow(x,3))", 0.0, 500.0);
  //E->SetParameter(0, p0);

  TFile *E_output = new TFile("E_spline.root", "recreate");
  E_spline->Write();
  E_output->Write();
  E_output->Close();
}
