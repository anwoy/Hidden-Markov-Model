#ifndef GAURD_defs_1
#define GAURD_defs_1

#define PI 3.14159265359
#define LOG_PI 1.1447298858494002
#define LOG_2 0.69314718055994529

#define A_PRIOR 1.0000001
#define Q0 10000000.0


using namespace std;

vector<string> split(const string &line,const char &separator){
	vector<string> ans;
	string::const_iterator i(line.begin()),j(line.begin());
	while(i!=line.end()){
		j=i;
		while(*j!=separator && j!=line.end())
			++j;
		string temp(i,j);
		ans.push_back(temp);
		
		i=j;
		if(i!=line.end())
			++i;
		if(j!=line.end() && i==line.end())
			ans.push_back(string(""));
	}
	
	return ans;
}

vector<string> split(const string &line){
	vector<string> ans;
	string::const_iterator i(line.begin()),j(line.begin());
	while(i!=line.end()){
		while(isspace(*i))
			++i;
		j=i;
		while(j!=line.end() && !isspace(*j))
			++j;
		string temp(i,j);
		if(temp.size()>0)
			ans.push_back(temp);
		i=j;
	}
	return ans;
}



template<class T>
ostream& operator<<(ostream &out,const vector<T> &V){
	for(size_t i(0);i<V.size()-1;++i)
		out<<V[i]<<",";
	out<<V[V.size()-1];
	return out;
}

template<class T>
ostream& operator<<(ostream &out,const vector<vector<T> > &V){
	for(int i(0);i<V.size();++i)
		out<<V[i]<<"\n";
	return out;
}


vector<vector<double> > columnextractor(const vector<vector<double> > &V,const pair<size_t,size_t> &P){
	if(P.first > V[0].size()-1 || P.second > V[0].size())
		throw(domain_error("indices are out of range in columnextractor"));
	
	vector<vector<double> > ans(V.size(),vector<double>(P.second-P.first));
	for(size_t i(0);i<V.size();++i){
		for(size_t j(P.first);j<P.second;++j)
			ans[i][j-P.first] = V[i][j];
	}
	return ans;
}




vector<double> extract(const vector<double> &V,const vector<bool> &index){
	if(index.size()!=V.size())
		throw(domain_error("ERROR1 in extract"));
	size_t count(0);
	for(size_t i(0);i<index.size();++i)
		if(index[i])
			count += 1;
	vector<double> ans(count);
	count=0;
	for(size_t i(0);i<index.size();++i){
		if(index[i]){
			ans[count]=V[i];
			++count;
		}
	}
	return ans;
}


vector<vector<double> > extract(const vector<vector<double> > &V,const vector<bool> &index1,const vector<bool> &index2){
	if(V.size()!=index1.size() || V[0].size()!=index2.size())
		throw(domain_error("ERROR2 in extract"));
	size_t count1(0),count2(0);
	for(size_t i(0);i<index1.size();++i)
		if(index1[i])
			count1 += 1;
	for(size_t i(0);i<index2.size();++i)
		if(index2[i])
			count2 += 1;
		
	vector<vector<double> > ans(count1,vector<double>(count2));
	count1=0;
	for(size_t i(0);i<index1.size();++i){
		if(index1[i]){
			count2=0;
			for(size_t j(0);j<index2.size();++j){
				if(index2[j]){
					ans[count1][count2] = V[i][j];
					count2 += 1;
				}
			}
			count1 += 1;
		}
	}
	return ans;
}





vector<double> LOG_1D(const vector<double> &V){
	vector<double> ans(V.size());
	for(size_t i(0);i<V.size();++i)
		ans[i] = log(V[i]);
	return ans;
}

void LOG_1D_INPLACE(vector<double> &V){
	for(size_t i(0);i<V.size();++i)
		V[i] = log(V[i]);
}


vector<vector<double> > LOG_2D(const vector<vector<double> > &V){
	vector<vector<double> > ans(V.size(),vector<double>(V[0].size()));
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[0].size();++j)
			ans[i][j] = log( V[i][j] );
	return ans;
}


vector<double> EXP_1D(const vector<double> &V){
	vector<double> ans(V.size());
	for(size_t i(0);i<V.size();++i)
		ans[i] = exp(V[i]);
	return ans;
}


void EXP_1D_INPLACE(vector<double> &V){
	for(size_t i(0);i<V.size();++i)
		V[i] = exp(V[i]);
}

vector<vector<double> > EYE(const size_t &D){
	vector<vector<double> > ans(D,vector<double>(D,0));
	for(size_t i(0);i<D;++i)
		ans[i][i]=1;
	return ans;
}

int determinant_sign(const boost::numeric::ublas::permutation_matrix<std ::size_t> &pm){
	int pm_sign=1;
	size_t size = pm.size();
	for (size_t i = 0; i < size; ++i)
		if (i != pm(i))
			pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
	return pm_sign;
}
 
double determinant( boost::numeric::ublas::matrix<double> &m ){
	boost::numeric::ublas::permutation_matrix<size_t> pm(m.size1());
	double det = 1.0;
	if( boost::numeric::ublas::lu_factorize(m,pm) ) {
		det = 0.0;
	}
	else {
		for(int i = 0; i < m.size1(); i++)
			det *= m(i,i); // multiply by elements on diagonal
		det = det * determinant_sign( pm );
	}
	return det;
}

double determinant(const vector<vector<double> > &cov){
	if(cov.size() != cov[0].size())
		throw(domain_error("ERROR1 in determinant"));
	const size_t N(cov.size());
	boost::numeric::ublas::matrix<double> M(N,N);
	for(size_t i(0);i<N;++i)
		for(size_t j(0);j<N;++j)
			M(i,j) = cov[i][j];
	return determinant(M);
}

template<class T>
bool InvertMatrix(const boost::numeric::ublas::matrix<T> &input, boost::numeric::ublas::matrix<T> &inverse){
	
	// create a working copy of the input
	boost::numeric::ublas::matrix<T> A(input);
	
	// create a permutation matrix for the LU-factorization
	boost::numeric::ublas::permutation_matrix<size_t> pm(A.size1());
	
	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;
	
	// create identity matrix of "inverse"
	inverse.assign(boost::numeric::ublas::identity_matrix<T> (A.size1()));
	
	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
	
	return true;
}


vector<vector<double> > InvertMatrix(const vector<vector<double> > &V){
	if(V.size() != V[0].size())
		throw(domain_error("ERROR1 in InvertMatrix"));
	const size_t N(V.size());
	boost::numeric::ublas::matrix<double> M(N,N),inverse(N,N);
	for(size_t i(0);i<N;++i)
		for(size_t j(0);j<N;++j)
			M(i,j)=V[i][j];
	InvertMatrix(M,inverse);
	
	vector<vector<double> > ans(N,vector<double>(N));
	for(size_t i(0);i<N;++i)
		for(size_t j(0);j<N;++j)
			ans[i][j] = inverse(i,j);
	return ans;
}


void testdata(const vector<vector<vector<double> > > &allX){
	const size_t Ncols( allX[0][0].size() );
	int count(0);
	bool mark(1);
	for(size_t i(0);i<allX.size();++i){
		for(size_t j(0);j<allX[i].size();++j){
			++count;
			if(allX[i][j].size() != Ncols)
				mark=0;
		}
	}
	if(mark)
		cout<<"test passed!!\n";
	else
		cout<<"test failed\n";
	cout<<"Ncols="<<Ncols<<"\n";
	cout<<"count="<<count<<"\n";
}

vector<double> DOT(const vector<double> &V,const vector<vector<double> > &M){
	if(V.size()!=M.size())
		throw(domain_error("ERROR1 in DOT"));
	
	vector<double> ans(M[0].size());
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += V[j]*M[j][i];
	}
	return ans;
}

void DOT(const vector<double> &V,const vector<vector<double> > &M,vector<double> &ans){
	if(V.size()!=M.size())
		throw(domain_error("ERROR12 in DOT"));
	if(M[0].size() != ans.size())
		throw(domain_error("ERROR123 in DOT"));
	
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += V[j]*M[j][i];
	}
}


vector<double> DOT(const vector<vector<double> > &M,const vector<double> &V){
	if(M[0].size()!=V.size())
		throw(domain_error("ERROR2 in DOT"));
	vector<double> ans(M.size());
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += M[i][j]*V[j];
	}
	return ans;
}


void DOT(const vector<vector<double> > &M,const vector<double> &V,vector<double> &ans){
	if(M[0].size()!=V.size())
		throw(domain_error("ERROR21 in DOT"));
	if(M.size() != ans.size())
		throw(domain_error("ERROR213 in DOT"));
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += M[i][j]*V[j];
	}
}



double DOT(const vector<double> &v1,const vector<double> &v2){
	if(v1.size() != v2.size())
		throw(domain_error("ERROR3 in DOT"));
	double ans(0);
	for(size_t i(0);i<v1.size();++i)
		ans += v1[i]*v2[i];
	return ans;
}
vector<vector<double> > DOT(const vector<vector<double> > &v1,const vector<vector<double> > &v2){
	if(v1[0].size()!=v2.size())
		throw(domain_error("ERROR4 in DOT"));
	vector<vector<double> > ans(v1.size(),vector<double>(v2[0].size()));
	for(size_t i(0);i<ans.size();++i){
		for(size_t j(0);j<ans[0].size();++j){
			ans[i][j] = 0;
			for(size_t k(0);k<v1[0].size();++k)
				ans[i][j] += v1[i][k]*v2[k][j];
		}
	}
	return ans;
}


vector<double> operator*(const double &S,const vector<double> &V){
	vector<double> ans(V.size());
	for(size_t i(0);i<ans.size();++i)
		ans[i] = S*V[i];
	return ans;
}

vector<vector<double> > operator*(const double &S,const vector<vector<double> > &V){
	vector<vector<double> > ans(V);
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[0].size();++j)
			ans[i][j] = S*V[i][j];
	return ans;
}

vector<double> operator+(const vector<double> &v1,const vector<double> &v2){
	if(v1.size() != v2.size())
		throw(domain_error("ERROR1 in operator+"));
	vector<double> ans(v1.size());
	for(size_t i(0);i<ans.size();++i)
		ans[i] = v1[i]+v2[i];
	return ans;
}

vector<vector<double> > operator+(const vector<vector<double> > &v1,const vector<vector<double> > &v2){
	if(v1.size()!=v2.size() || v1[0].size()!=v2[0].size())
		throw(domain_error("ERROR in vector<vector<double> > operator+"));
	
	vector<vector<double> > ans(v1.size(),vector<double>(v1[0].size()));
	for(size_t k1(0);k1<v1.size();++k1)
		for(size_t k2(0);k2<v1[0].size();++k2)
			ans[k1][k2] = v1[k1][k2]+v2[k1][k2];
	
	return ans;
}

vector<double> operator+(const vector<double> &V,const double &S){
	vector<double> ans(V.size());
	for(size_t i(0);i<V.size();++i)
		ans[i] = V[i]+S;
	return ans;
}


void operator-=(vector<double> &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		V[i] = V[i]-S;
}

void operator+=(vector<double> &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		V[i] = V[i]+S;
}

void operator*=(vector<double> &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		V[i] = V[i]*S;
}


double MAX(const vector<double> &V){
	double ans(V[0]);
	for(size_t i(0);i<V.size();++i)
		if(V[i]>ans)
			ans=V[i];
	return ans;
}

double MAX_2D(const vector<vector<double> > &V){
	double ans(V[0][0]);
	for(size_t i(0);i<V.size();++i){
		for(size_t j(0);j<V[0].size();++j){
			if(V[i][j]>ans)
				ans = V[i][j];
		}
	}
	return ans;
}


size_t ARGMAX(const vector<double> &V){
	size_t ans(0);
	for(size_t i(0);i<V.size();++i)
		if(V[i]>V[ans])
			ans=i;
	return ans;
}


double SUM(const vector<double> &V){
	double ans(0);
	for(size_t i(0);i<V.size();++i)
		ans += V[i];
	return ans;
}

double SUM_2D(const vector<vector<double> > &V){
	double ans(0);
	for(size_t i(0);i<V.size();++i)
		ans += SUM(V[i]);
	return ans;
}

void savevectors(const vector<double> &V,ofstream &fout){
	const size_t d1( V.size() );
	fout.write((char*)&d1,sizeof(size_t));
	for(size_t i(0);i<d1;++i)
		fout.write((char*)&V[i],sizeof(double));
}

void savevectors(const vector<vector<double> > &V,ofstream &fout){
	const size_t d1( V.size() ),d2( V[0].size() );
	fout.write((char*)&d1,sizeof(size_t));
	fout.write((char*)&d2,sizeof(size_t));
	for(size_t i(0);i<d1;++i)
		for(size_t j(0);j<d2;++j)
			fout.write((char*)&V[i][j],sizeof(double));
}

void savevectors(const vector<vector<vector<double> > > &V,ofstream &fout){
	const size_t d1( V.size() ),d2( V[0].size() ),d3( V[0][0].size() );
	fout.write((char*)&d1,sizeof(size_t));
	fout.write((char*)&d2,sizeof(size_t));
	fout.write((char*)&d3,sizeof(size_t));
	for(size_t i(0);i<d1;++i)
		for(size_t j(0);j<d2;++j)
			for(size_t k(0);k<d3;++k)
				fout.write((char*)&V[i][j][k],sizeof(double));
}

void loadvectors(vector<double> &V,ifstream &fin){
	size_t d1;
	fin.read((char*)&d1,sizeof(size_t));
	V = vector<double>(d1);
	for(size_t i(0);i<d1;++i)
		fin.read((char*)&V[i],sizeof(double));
}

void loadvectors(vector<vector<double> > &V,ifstream &fin){
	size_t d1,d2;
	fin.read((char*)&d1,sizeof(size_t));
	fin.read((char*)&d2,sizeof(size_t));
	V = vector<vector<double> >(d1,vector<double>(d2));
	for(size_t i(0);i<d1;++i)
		for(size_t j(0);j<d2;++j)
			fin.read((char*)&V[i][j],sizeof(double));
}

void loadvectors(vector<vector<vector<double> > > &V,ifstream &fin){
	size_t d1,d2,d3;
	fin.read((char*)&d1,sizeof(size_t));
	fin.read((char*)&d2,sizeof(size_t));
	fin.read((char*)&d3,sizeof(size_t));
	V = vector<vector<vector<double> > >(d1,vector<vector<double> >(d2,vector<double>(d3)));
	for(size_t i(0);i<d1;++i)
		for(size_t j(0);j<d2;++j)
			for(size_t k(0);k<d3;++k)
				fin.read((char*)&V[i][j][k],sizeof(double));
}


class probties{
	public:
	size_t Nstates,D;
	vector<vector<double> > transition_probties;
	vector<vector<double> > transition_probties_log;
	vector<vector<double> > mu;
	vector<vector<vector<double> > > cov;
	vector<vector<vector<double> > > covinv;
	vector<double> detcov;
	
	probties() {}
	probties(const string&);
	probties(const size_t &nstates,const size_t &d){ initialize(nstates,d); }
	
	void initialize(const size_t &,const size_t &);
	void transition_probties_initialize();
	
	void savebinary(ofstream &) const;
	void loadbinary(ifstream &);
};

bool operator == (const probties &O1,const probties &O2) {
	if(O1.Nstates == O2.Nstates and O1.D == O2.D)
		if(O1.transition_probties == O2.transition_probties)
			if(O1.transition_probties_log == O2.transition_probties_log)
				if(O1.mu == O2.mu and O1.cov == O2.cov and O1.covinv == O2.covinv and O1.detcov == O2.detcov)
					return true;
	else
		return false;
}

void probties::loadbinary(ifstream &fin){
	//ifstream fin(filename.c_str() , ios::in | ios::binary);
	fin.read((char*)&Nstates,sizeof(size_t));
	fin.read((char*)&D,sizeof(size_t));
	loadvectors(transition_probties,fin);
	loadvectors(transition_probties_log,fin);
	loadvectors(mu,fin);
	loadvectors(cov,fin);
	loadvectors(covinv,fin);
	loadvectors(detcov,fin);
	fin.close();
}

void probties::savebinary(ofstream &fout) const {
	//ofstream fout(filename.c_str() , ios::out | ios::binary);
	fout.write((char*)&Nstates,sizeof(size_t));
	fout.write((char*)&D,sizeof(size_t));
	savevectors(transition_probties,fout);
	savevectors(transition_probties_log,fout);
	savevectors(mu,fout);
	savevectors(cov,fout);
	savevectors(covinv,fout);
	savevectors(detcov,fout);
	fout.close();
}

void probties::transition_probties_initialize(){
	srand(time(NULL));
	transition_probties = vector<vector<double> >(Nstates,vector<double>(Nstates));
	vector<double> colsums(Nstates,0);
	for(size_t i(0);i<Nstates;++i){
		for(size_t j(0);j<Nstates;++j){
			transition_probties[i][j] = rand()*1.0/RAND_MAX;
			colsums[i] += transition_probties[i][j];
		}
	}
	
	for(size_t i(0);i<Nstates;++i)
		for(size_t j(0);j<Nstates;++j)
			transition_probties[i][j] /= colsums[i];
}

void probties::initialize(const size_t &nstates,const size_t &d){
	Nstates = nstates;
	D = d;
	transition_probties_initialize();
	transition_probties_log = LOG_2D( transition_probties );
	mu = vector<vector<double> >(Nstates,vector<double>(D,0));
	cov = vector<vector<vector<double> > >(Nstates);
	for(size_t i(0);i<Nstates;++i)
		cov[i] = EYE(D);
	covinv = vector<vector<vector<double> > >(Nstates);
	for(size_t i(0);i<Nstates;++i)
		covinv[i] = EYE(D);
	detcov = vector<double>(Nstates,1);
}


probties::probties(const string &filename){
	ifstream fin(filename.c_str());
	string line;
	getline(fin,line);
	vector<string> split_line( split(line) );
	Nstates = atoi( split_line[0].c_str() );
	
	getline(fin,line);
	split_line = split(line);
	D = atoi( split_line[0].c_str() );
	
	getline(fin,line);
	getline(fin,line);
	for(int i(0);i<Nstates;++i){
		getline(fin,line);
		split_line = split(line,',');
		if(split_line.size() != Nstates)
			throw(domain_error("ERROR1 in probties::probties"));
		vector<double> row;
		for(int j(0);j<Nstates;++j)
			row.push_back( atof(split_line[j].c_str()) );
		transition_probties.push_back(row);
	}
	
	getline(fin,line);
	getline(fin,line);
	for(int i(0);i<Nstates;++i){
		getline(fin,line);
		split_line = split(line,',');
		if(split_line.size() != D)
			throw(domain_error("ERROR2 in probties::probties"));
		vector<double> row;
		for(int j(0);j<D;++j)
			row.push_back( atof(split_line[j].c_str()) );
		mu.push_back(row);
	}
	
	getline(fin,line);
	getline(fin,line);
	for(int i(0);i<Nstates;++i){
		vector<vector<double> > covval,covvalinv;
		for(int j(0);j<D;++j){
			getline(fin,line);
			split_line = split(line,',');
			if(split_line.size() != D)
				throw(domain_error("ERROR3 in probties::probties"));
			vector<double> row;
			for(int k(0);k<D;++k)
				row.push_back(atof(split_line[k].c_str()));
			covval.push_back(row);
		}
		cov.push_back(covval);
		covvalinv=InvertMatrix(covval);
		covinv.push_back(covvalinv);
		detcov.push_back(determinant(covval));
		
		getline(fin,line);
	}
	
	transition_probties_log = LOG_2D(transition_probties);
	
	
}

ostream& operator<<(ostream &out,const probties &O){
	out<<"Nstates,D:\n";
	out<<O.Nstates<<"\t"<<O.D<<"\n";
	out<<"transition_probties:\n";
	out<<O.transition_probties;
	out<<"transition_probties_log:\n";
	out<<O.transition_probties_log;
	out<<"mu:\n";
	out<<O.mu;
	out<<"cov:\n";
	for(int i(0);i<O.cov.size();++i)
		out<<O.cov[i]<<"\n";
	out<<"covinv:\n";
	for(int i(0);i<O.covinv.size();++i)
		out<<O.covinv[i]<<"\n";
	out<<"detcov:\n";
	out<<O.detcov<<"\n";
	
	out<<"matrix multiplications:\n";
	for(size_t i(0);i<O.cov.size();++i)
		cout<< DOT(O.cov[i],O.covinv[i]) <<"\n";
}



void forward_backward_aux(const vector<double> &row,const probties &O,vector<double> &mess){
	if(mess.size()!=O.Nstates)
		throw(domain_error("ERROR1 in forward_backward_aux"));
	
	vector<bool> nonnanindex(row.size(),0);
	size_t nancount(0);
	for(size_t i(0);i<row.size();++i){
		if(!isnan(row[i]))
			nonnanindex[i]=1;
		else{
			nonnanindex[i]=0;
			nancount += 1;
		}
	}
	
	
	if(nancount==0){
		for(size_t i(0);i<O.Nstates;++i){
			const vector<double> temp(row+((-1)*O.mu[i]));
			mess[i] = -int(O.D)/2.0*(LOG_2+LOG_PI) - 0.5*log(O.detcov[i]) - 0.5*DOT(DOT(temp,O.covinv[i]),temp);
		}
	}
	else{
		if(nancount==O.D){
			for(size_t i(0);i<O.Nstates;++i)
				mess[i]=0;
		}
		else{
			for(size_t i(0);i<O.Nstates;++i){
				const vector<double> newmu( extract(O.mu[i],nonnanindex) );
				const vector<vector<double> > newcov( extract(O.cov[i],nonnanindex,nonnanindex) );
				const vector<vector<double> > newcovinv( InvertMatrix(newcov) );
				const double newdetcov( determinant(newcov) );
				const vector<double> newrow( extract(row,nonnanindex) );
				const vector<double> temp( newrow+((-1)*newmu) );
				const size_t newD( O.D - nancount );
				mess[i] = -int(newD)/2.0*(LOG_2+LOG_PI) - 0.5*log(newdetcov) -0.5*DOT(DOT(temp,newcovinv),temp);
			}
		}
	}
}


void ADD(const vector<double> &v1,const vector<double> &v2,vector<double> &ans){
	if(v1.size()!=v2.size() || v1.size()!=ans.size())
		throw(domain_error("ERROR in ADD"));
	for(size_t i(0);i<ans.size();++i)
		ans[i] = v1[i]+v2[i];
}

/*
-0->   -2->
<-1-   <-3-
	^
	| |
	4 5
	| |
*/

void forward_backward(const probties &O,const vector<vector<double> > &X,vector<vector<vector<double> > > &Mess){
	if(Mess.size()!=X.size())
		throw(domain_error("ERROR1 in forward_backward"));
	
	const size_t N( X.size() );
	
	for(size_t i(0);i<O.Nstates;++i)
		Mess[0][0][i] = log(1.0/O.Nstates);
	
	for(size_t i(0);i<N;++i){
		if(i!=0){
			vector<double> temp( Mess[i-1][2] );
			const double norm( MAX(temp) );
			temp -= norm;
			EXP_1D_INPLACE(temp);
			DOT(O.transition_probties,temp,Mess[i][0]);
			LOG_1D_INPLACE(Mess[i][0]);
			Mess[i][0] += norm;
		}
		forward_backward_aux(X[i],O,Mess[i][4]);
		ADD(Mess[i][0],Mess[i][4],Mess[i][2]);
	}
	
	Mess[N-1][3] = vector<double>(O.Nstates,0);
	for(size_t i(N-1);i>0;--i){
		if(i!=N-1){
			vector<double> temp(Mess[i+1][1]);
			const double norm( MAX(temp) );
			temp -= norm;
			EXP_1D_INPLACE(temp);
			DOT(temp,O.transition_probties,Mess[i][3]);
			LOG_1D_INPLACE(Mess[i][3]);
			Mess[i][3] += norm;
		}
		ADD(Mess[i][0],Mess[i][3],Mess[i][5]);
		ADD(Mess[i][4],Mess[i][3],Mess[i][1]);
	}
	{
		const size_t i(0);
		if(i!=N-1){
			vector<double> temp(Mess[i+1][1]);
			const double norm( MAX(temp) );
			temp -= norm;
			EXP_1D_INPLACE(temp);
			DOT(temp,O.transition_probties,Mess[i][3]);
			LOG_1D_INPLACE(Mess[i][3]);
			Mess[i][3] += norm;
		}
		ADD(Mess[i][0],Mess[i][3],Mess[i][5]);
		ADD(Mess[i][4],Mess[i][3],Mess[i][1]);
	}
}

struct wrapper{
	vector<vector<vector<double> > > allX;
	
	void readdata_csv(const string &);
	
	double getD() const;
};

double wrapper::getD() const{
	const vector<vector<double> > &X(allX[0]);
	const vector<double> &row(X[0]);
	return row.size();
}


void wrapper::readdata_csv(const string &filename){
	ifstream fin(filename.c_str());
	vector<vector<double> > X(0);
	while(!fin.eof()){
		string line;
		getline(fin,line);
		if(line.size()>0){
			vector<string> split_line(split(line,','));
			vector<double> row(split_line.size());
			for(size_t i(0);i<row.size();++i)
				row[i] = atof( split_line[i].c_str() );
			X.push_back(row);
		}
		else{
			if(X.size()>0){
				allX.push_back(X);
				X = vector<vector<double> >(0);
			}
		}
	}
}


struct Messdata{
	vector<vector<vector<vector<double> > > > allMess;
	
	Messdata(const vector<vector<vector<double> > > &allX,const probties &O);
	void FB(const probties &,const wrapper &);
	void predstates(const string &) const;
};

/*
-0->   -2->
<-1-   <-3-
	^
	| |
	4 5
	| |
*/

void Messdata::predstates(const string &filename) const {
	ofstream fout(filename.c_str());
	for(size_t i(0);i<allMess.size();++i){
		const vector<vector<vector<double> > > &Mess( allMess[i] );
		for(size_t k(0);k<Mess.size();++k)
			fout<< ARGMAX( Mess[k][0] + Mess[k][4] + Mess[k][3] ) <<"\n";
		fout<<"\n";
	}
	fout.close();
}


void Messdata::FB(const probties &O,const wrapper &WO){
	if(allMess.size()!=WO.allX.size())
		throw(domain_error("ERROR in Messdata::forward_backward"));
	
	for(size_t i(0);i<allMess.size();++i)
		forward_backward(O,WO.allX[i],allMess[i]);
}


ostream& operator<<(ostream &out,const Messdata &MDO){
	cout<<MDO.allMess.size()<<"\n";
	for(size_t i(0);i<MDO.allMess.size();++i){
		const vector<vector<vector<double> > > &Mess(MDO.allMess[i]);
		for(size_t j(0);j<Mess.size();++j){
			cout<<j<<")\n";
			cout<<"***********\n";
			cout<<Mess[j]<<"\n";
		}
	}
}



double prior_contribution(const probties &O){
	double ans(0);
	for(size_t i(0);i<O.Nstates;++i)
		for(size_t j(0);j<O.Nstates;++j)
			ans += (A_PRIOR-1)*log(O.transition_probties[i][j]);
	
	const double V0( O.D+1 );
	
	for(size_t i(0);i<O.Nstates;++i){
		ans += -0.5*log(O.detcov[i])-0.5*DOT(DOT(O.mu[i],O.covinv[i]),O.mu[i])/Q0;
		ans += -(V0+O.D+1)/2.0*log(O.detcov[i]);
		for(size_t j(0);j<O.D;++j)
			ans += -0.5*O.covinv[i][j][j];
	}
	
	return ans;
}




Messdata::Messdata(const vector<vector<vector<double> > > &allX,const probties &O){
	allMess = vector<vector<vector<vector<double> > > >(allX.size());
	for(size_t i(0);i<allX.size();++i){
		const vector<vector<double> > &X(allX[i]);
		const vector<vector<vector<double> > > Mess(X.size(),vector<vector<double> >(6,vector<double>(O.Nstates)));
		allMess[i] = Mess;
	}
}

vector<double> mucalc_aux(const vector<double> &row,const probties &O,const size_t STATE){
	if(row.size()!=O.D)
		throw(domain_error("ERROR in mucalc_aux"));
	
	vector<bool> nonnanindex(O.D),nanindex(O.D);
	size_t nancount(0);
	for(size_t i(0);i<row.size();++i){
		if(isnan(row[i])){
			nancount += 1;
			nanindex[i] = 1;
			nonnanindex[i] = 0;
		}
		else{
			nanindex[i] = 0;
			nonnanindex[i] = 1;
		}
	}
	
	if(nancount==0)
		return row;
	else
		if(nancount==O.D)
			return O.mu[STATE];
		else{
				
			const vector<double> mu_a( extract(O.mu[STATE],nanindex) );
			const vector<vector<double> > covab( extract(O.cov[STATE],nanindex,nonnanindex) );
			const vector<vector<double> > covbb( extract(O.cov[STATE],nonnanindex,nonnanindex) );
			const vector<vector<double> > covbbinv( InvertMatrix(covbb) );
			const vector<double> rowb( extract(row,nonnanindex) );
			const vector<double> mu_b( extract(O.mu[STATE],nonnanindex) );
			const vector<double> newval( mu_a+DOT(covab,DOT(covbbinv,rowb+((-1)*mu_b))) );
			
			vector<double> ans(row);
			size_t indexer(0);
			for(size_t i(0);i<ans.size();++i){
				if(isnan(ans[i])){
					ans[i]=newval[indexer];
					indexer += 1;
				}
			}
			
			return ans;
			
		}
}


vector<vector<double> > covcalc_aux(const vector<double> &row,const probties &O,const vector<vector<double> > &MU,const size_t STATE){
	if(row.size()!=O.D)
		throw(domain_error("ERROR in covcalc_aux"));
	
	if(MU.size()!=O.Nstates || MU[0].size()!=O.D)
		throw(domain_error("ERROR2 in covcalc_aux"));
	
	vector<bool> nonnanindex(O.D),nanindex(O.D);
	size_t nancount(0);
	for(size_t i(0);i<row.size();++i){
		if(isnan(row[i])){
			nancount += 1;
			nanindex[i] = 1;
			nonnanindex[i] = 0;
		}
		else{
			nanindex[i] = 0;
			nonnanindex[i] = 1;
		}
	}
	
	
	if(nancount==0){
		const vector<double> temp(row+((-1)*MU[STATE]));
		vector<vector<double> > ans(O.D,vector<double>(O.D));
		for(size_t k1(0);k1<O.D;++k1)
			for(size_t k2(0);k2<O.D;++k2)
				ans[k1][k2] = temp[k1]*temp[k2];
		return ans;
	}
	else
		if(nancount==O.D){
			vector<vector<double> > ans(O.D,vector<double>(O.D));
			for(size_t k1(0);k1<O.D;++k1)
				for(size_t k2(0);k2<O.D;++k2)
					ans[k1][k2] = O.cov[STATE][k1][k2]+O.mu[STATE][k1]*O.mu[STATE][k2]-O.mu[STATE][k1]*MU[STATE][k2]-MU[STATE][k1]*O.mu[STATE][k2]+MU[STATE][k1]*MU[STATE][k2];
			
			return ans;
		}
		else{
			const vector<double> mu_a( extract(O.mu[STATE],nanindex) );
			const vector<vector<double> > covab( extract(O.cov[STATE],nanindex,nonnanindex) );
			const vector<vector<double> > covba( extract(O.cov[STATE],nonnanindex,nanindex) );
			const vector<vector<double> > covbb( extract(O.cov[STATE],nonnanindex,nonnanindex) );
			const vector<vector<double> > covbbinv( InvertMatrix(covbb) );
			const vector<vector<double> > covaa( extract(O.cov[STATE],nanindex,nanindex) );
			const vector<double> rowb( extract(row,nonnanindex) );
			const vector<double> mu_b( extract(O.mu[STATE],nonnanindex) );
			
			const vector<double> newval( mu_a+DOT(covab,DOT(covbbinv,rowb+((-1)*mu_b))) );
			const vector<vector<double> > covnew( covaa+((-1)*DOT(DOT(covab,covbbinv),covba)) );
			
			
			vector<vector<double> > ans(O.D,vector<double>(O.D,0));
			size_t I1(0),I2(0);
			for(size_t k1(0);k1<O.D;++k1){
				for(size_t k2(0);k2<O.D;++k2){
					
					if(nanindex[k1]==0 && nanindex[k2]==0){
						ans[k1][k2] += (row[k1]-MU[STATE][k1])*(row[k2]-MU[STATE][k2]);
						//cout<<"case0\t"<<k1<<"\t"<<k2<<"\t"<<I1<<"\t"<<I2<<"\n";
					}
					if(nanindex[k1]==1 && nanindex[k2]==0){
						ans[k1][k2] += newval[I1]*row[k2]-newval[I1]*MU[STATE][k2]-MU[STATE][k1]*row[k2]+MU[STATE][k1]*MU[STATE][k2];
						//cout<<"case1\t"<<k1<<"\t"<<k2<<"\t"<<I1<<"\t"<<I2<<"\n";
					}
					if(nanindex[k1]==0 && nanindex[k2]==1){
						ans[k1][k2] += row[k1]*newval[I2]-row[k1]*MU[STATE][k2]-MU[STATE][k1]*newval[I2]+MU[STATE][k1]*MU[STATE][k2];
						//cout<<"case2\t"<<k1<<"\t"<<k2<<"\t"<<I1<<"\t"<<I2<<"\n";
					}
					if(nanindex[k1]==1 && nanindex[k2]==1){
						ans[k1][k2] += covnew[I1][I2]+newval[I1]*newval[I2]-newval[I1]*MU[STATE][k2]-MU[STATE][k1]*newval[I2]+MU[STATE][k1]*MU[STATE][k2];
						//cout<<"case3\t"<<k1<<"\t"<<k2<<"\t"<<I1<<"\t"<<I2<<"\n";
					}
					
					if(nanindex[k2]==1)
						I2 += 1;
					if(I2 == nancount)
						I2=0;
					
				}
				
				if(nanindex[k1]==1)
					I1 += 1;
				
			}
			return ans;
		}
}


probties calcnewprobties(const Messdata &MDO,const vector<vector<vector<double> > > &allX,const probties &O){
	// assume forward_backward is already run
	const double V0(O.D+1);
	if(MDO.allMess.size() != allX.size())
		throw(domain_error("ERROR1 in calcnewprobties"));
	
	vector<vector<double> > transition_probties(O.Nstates,vector<double>(O.Nstates,A_PRIOR-1));
	vector<vector<double> > mu(O.Nstates,vector<double>(O.D,0));
	vector<vector<double> > munum(O.Nstates,vector<double>(O.D,0));
	vector<double> muden(O.Nstates,1.0/Q0);
	for(size_t i(0);i<MDO.allMess.size();++i){
		const vector<vector<vector<double> > > &Mess( MDO.allMess[i] );
		const vector<vector<double> > &X( allX[i] );
		if(X.size() != Mess.size())
			throw(domain_error("ERROR2 in calcnewprobties"));
		
		for(size_t j(0);j<Mess.size();++j){
			if(j!=0){
				vector<vector<double> > temp(O.transition_probties_log);
				double norm;
				for(size_t k1(0);k1<O.Nstates;++k1){
					for(size_t k2(0);k2<O.Nstates;++k2){
						temp[k1][k2] += Mess[j][3][k1]+Mess[j][4][k1] + Mess[j-1][0][k2]+Mess[j-1][4][k2];
						if(k1==0 && k2==0)
							norm=temp[k1][k2];
						else
							if(temp[k1][k2] > norm)
								norm=temp[k1][k2];
					}
				}
				double summ(0);
				for(size_t k1(0);k1<O.Nstates;++k1)
					for(size_t k2(0);k2<O.Nstates;++k2){
						temp[k1][k2] -= norm;
						temp[k1][k2] = exp( temp[k1][k2] );
						summ += temp[k1][k2];
					}
				for(size_t k1(0);k1<O.Nstates;++k1)
					for(size_t k2(0);k2<O.Nstates;++k2){
						temp[k1][k2] /= summ;
						transition_probties[k1][k2] += temp[k1][k2];
					}
			}
			
			{
				vector<double> temp( Mess[j][0]+Mess[j][4]+Mess[j][3] );
				if(temp.size()!=O.Nstates)
					throw(domain_error("ERROR2 in calcnewprobties"));
				const double norm( MAX(temp) );
				temp -= norm;
				EXP_1D_INPLACE(temp);
				const double summ( SUM(temp) );
				temp *= 1.0/summ;
				for(size_t k1(0);k1<O.Nstates;++k1){
					const vector<double> newrow( mucalc_aux(X[j],O,k1) );
					for(size_t k2(0);k2<O.D;++k2)
						munum[k1][k2] += temp[k1]*newrow[k2];
					muden[k1] += temp[k1];
				}
			}
			
		}
	}
	
	for(size_t i(0);i<O.Nstates;++i){
		double summ(0);
		for(size_t j(0);j<O.Nstates;++j)
			summ += transition_probties[j][i];
		for(size_t j(0);j<O.Nstates;++j)
			transition_probties[j][i] /= summ;
	}
	
	for(size_t i(0);i<O.Nstates;++i)
		for(size_t j(0);j<O.D;++j)
			mu[i][j] = munum[i][j]/muden[i];
	
	const vector<double> covden(muden+(-1.0/Q0+V0-O.D));
	vector<vector<vector<double> > > covnum(O.Nstates,vector<vector<double> >(O.D,vector<double>(O.D)));
	for(size_t i(0);i<O.Nstates;++i){
		for(size_t k1(0);k1<O.D;++k1){
			for(size_t k2(0);k2<O.D;++k2){
				covnum[i][k1][k2] = 1.0/Q0*mu[i][k1]*mu[i][k2];
				if(k1==k2)
					covnum[i][k1][k2] += 1;
			}
		}
	}
	
	
	for(size_t i(0);i<MDO.allMess.size();++i){
		const vector<vector<vector<double> > > &Mess( MDO.allMess[i] );
		const vector<vector<double> > &X( allX[i] );
		if(X.size() != Mess.size())
			throw(domain_error("ERROR2 in calcnewprobties"));
		
		for(size_t j(0);j<Mess.size();++j){
			vector<double> temp( Mess[j][0]+Mess[j][4]+Mess[j][3] );
			if(temp.size()!=O.Nstates)
				throw(domain_error("ERROR2 in calcnewprobties"));
			const double norm( MAX(temp) );
			temp -= norm;
			EXP_1D_INPLACE(temp);
			const double summ( SUM(temp) );
			temp *= 1.0/summ;
			
			for(size_t KK(0);KK<O.Nstates;++KK){
				const vector<vector<double> > MAT( covcalc_aux(X[j],O,mu,KK) );
				for(size_t k1(0);k1<O.D;++k1){
					for(size_t k2(0);k2<O.D;++k2){
						covnum[KK][k1][k2] += temp[KK]*MAT[k1][k2];
					}
				}
			}
		}
	}
	
	vector<vector<vector<double> > > cov(covnum);
	
	for(size_t i(0);i<O.Nstates;++i)
		for(size_t k1(0);k1<O.D;++k1)
			for(size_t k2(0);k2<O.D;++k2)
				cov[i][k1][k2] = covnum[i][k1][k2]/covden[i];
	
	probties newO(O);
	newO.transition_probties = transition_probties;
	newO.transition_probties_log = LOG_2D( transition_probties );
	newO.mu = mu;
	newO.cov = cov;
	for(size_t i(0);i<newO.Nstates;++i){
		newO.covinv[i] = InvertMatrix( cov[i] );
		newO.detcov[i] = determinant( cov[i] );
	}
	
	return newO;
}




double lhoodcalc(const vector<vector<vector<double> > > &Mess,size_t pivot){
	if(pivot>=Mess.size())
		pivot = Mess.size()-1;
	
	vector<double> temp( Mess[pivot][0] + Mess[pivot][4] + Mess[pivot][3] );
	const double norm( MAX(temp) );
	temp -= norm;
	EXP_1D_INPLACE(temp);
	const double summ( SUM(temp) );
	
	return norm+log(summ);
}


struct EMhandler{
	vector<double> postprobties;
	probties O;
	Messdata MDO;
	
	EMhandler(const wrapper &WO,const size_t &Nstates): O(Nstates,WO.getD()) , MDO(WO.allX,O) {}
	
	void EMworker(const wrapper &);
};


void EMhandler::EMworker(const wrapper &WO){
	for(size_t i(0);i<4;++i){
		MDO.FB(O,WO);
		double lhood( prior_contribution(O) );
		for(size_t j(0);j<WO.allX.size();++j)
			lhood += lhoodcalc(MDO.allMess[j],10);
		postprobties.push_back(lhood);
		O = calcnewprobties(MDO,WO.allX,O);
	}
	
	vector<double>::reverse_iterator I(postprobties.rbegin());
	
	while(1){
		MDO.FB(O,WO);
		double lhood( prior_contribution(O) );
		for(size_t j(0);j<WO.allX.size();++j)
			lhood += lhoodcalc(MDO.allMess[j],10);
		if(abs((*I)-lhood) < 10e-7)
			break;
		postprobties.push_back(lhood);
		I = postprobties.rbegin();
		O = calcnewprobties(MDO,WO.allX,O);
	}
	
}

/////////////////

void learn(const string &datafile,const int Nstates,const string &paramsfile){
	wrapper WO;
	WO.readdata_csv(datafile);
	EMhandler EMHO(WO,Nstates);
	EMHO.EMworker(WO);
	
	ofstream fout(paramsfile.c_str() , ios::out | ios::binary);
	savevectors(EMHO.postprobties,fout);
	EMHO.O.savebinary(fout);
	fout.close();
}

void display(const string &paramsfile){
	vector<double> PP;
	probties O;
	
	ifstream fin(paramsfile.c_str() , ios::in | ios::binary);
	loadvectors(PP,fin);
	O.loadbinary(fin);
	fin.close();
	
	cout.precision(30);
	cout<<PP<<"\n";
	cout<<O;
}


void prediction(const string &paramsfile,const string &datafile,const string &outfile){
	wrapper WO;
	WO.readdata_csv(datafile);
	
	vector<double> PP;
	probties O;
	
	ifstream fin(paramsfile.c_str() , ios::in | ios::binary);
	loadvectors(PP,fin);
	O.loadbinary(fin);
	fin.close();
	
	Messdata MDO(WO.allX,O);
	MDO.FB(O,WO);
	MDO.predstates(outfile);
}




#endif














