

//return (x, ii, abs(newNorm-oldNorm), res.real)
//A,x,f, do split_b, tolerance,my,lambda,rows,columns
std::tuple<std::vector<double>,int,double,double> doOptimize(std::vector<double>,std::vector<double> ,std::vector<double>, bool, double,double,double,int,int);

//def objectivefunctionDer(A,x, f, mu, lmbda, d, b, AtA, ftA,rows,columns);


double minObjFunc(std::vector<double>,std::vector<double>,std::vector<double>,double,double,std::vector<double>,std::vector<double>,std::vector<int>,std::vector<double>,int,int);

std::vector<double> minObjFuncDer(std::vector<double>,std::vector<double>,std::vector<double>,double,double,std::vector<double>,std::vector<double>,std::vector<int>,std::vector<double>,int,int);
