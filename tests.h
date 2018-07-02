#ifndef GAURD_tests_1
#define GAURD_tests_1

using namespace std;


void tester_01(){
	wrapper WO;
	WO.readdata_csv("data.csv");
	probties O("probties.txt");
	Messdata MDO(WO.allX,O);
	calcnewprobties(MDO,WO.allX,O);
}

void tester_02(){
	wrapper WO;
	WO.readdata_csv("data.csv");
	probties O(6,4);
	Messdata MDO(WO.allX,O);
	for(size_t i(0);i<500;++i){
		MDO.FB(O,WO);
		double lhood( prior_contribution(O) );
		for(size_t j(0);j<WO.allX.size();++j)
			lhood += lhoodcalc(MDO.allMess[j],10);
		cout.precision(30);
		cout<<lhood<<",";
		cout.flush();
		O = calcnewprobties(MDO,WO.allX,O);
	}
	cout<<"\n";
	cout<<O;
	for(size_t i(0);i<WO.allX.size();++i)
		cout<< lhoodcalc(MDO.allMess[i],10) <<"\t"<< lhoodcalc(MDO.allMess[i],1000) <<"\n";
		
}

void tester_03(){
	wrapper WO;
	WO.readdata_csv("data.csv");
	EMhandler EMHO(WO,3);
	EMHO.EMworker(WO);
	cout.precision(30);
	cout<<EMHO.postprobties<<"\n";
	cout<<EMHO.O;
	{
		ofstream fout("MO.dat" , ios::out | ios::binary);
		savevectors(EMHO.postprobties,fout);
		EMHO.O.savebinary(fout);
		fout.close();
	}
	probties Onew;
	vector<double> PPnew;
	{
		ifstream fin("MO.dat" , ios::in | ios::binary);
		loadvectors(PPnew,fin);
		Onew.loadbinary(fin);
		fin.close();
	}
	if(!(Onew == EMHO.O and EMHO.postprobties == PPnew))
		throw(domain_error("ERROR1 in tester_03"));
	else
		cout<<"saving and loading works fine"<<"\n";
	
	// testing prediction of hidden states
	const probties AP("probties.txt");
	Messdata MDO(WO.allX,AP);
	MDO.FB(AP,WO);
	MDO.predstates("states_new.txt");
}

void tester_04(){
	learn("data.csv",3,"MO.dat");
	display("MO.dat");
	prediction("MO.dat","data.csv","states_new.txt");
}



#endif



















