// Root script to take a data file and return a spline created from
// that data. The data file must be formatted with # followed by the
// number of knots (number of data points) on the first line,
// labels for the two columns separated by a space on the second line,
// and data on the remainin lines. For example:

// #3
// #E    xsec
// 0.0   0.0
// 0.5   1.2e-11
// 1.3   3.0e-11

// A pointer to the TSpline3 object made from the points will be returned.

TSpline3 txt2spline(TString file = "data.txt"){
  ifstream in;
  in.open("data.txt");

  // Get number of knots from the first line
  TString temp;
  in >> temp;
  temp = temp(1,temp.Length());
  Int_t nknots = temp.Atoi();
  //nknots = 20;

  // ignore the second line
  in >> temp;
  in >> temp;

  // put the data in two arrays of doubles
  TArrayD x(nknots);
  TArrayD y(nknots);
  for(int i=0; i<nknots; i++){
    in >> x[i];
    in >> y[i];
    cout << x[i] << ", " << y[i] << endl;
  }

  //TGraph tempGraph(nknots, x.GetArray(), y.GetArray());

  TSpline3 * spl = new TSpline3("title",x.GetArray(),y.GetArray(),nknots);

  return *spl;
}
