#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "definitions.h"
#include "tests.h"

using namespace std;


int main(int argc,char *argv[]){
	if(argc > 1){
		if(string(argv[1]) == string("--learn")){
			if(argc == 5){
				const string datafile(argv[2]);
				const int Nstates( atoi( argv[3] ) );
				const string paramsfile( argv[4] );
				learn(datafile,Nstates,paramsfile);
				return 0;
			}
			else{
				cout<<"usage:\n";
				cout<<"hmm --learn <datafile> <Nstates> <paramsfile>\n";
				return 0;
			}
		}
		
		if(string(argv[1]) == string("--display")){
			if(argc == 3){
				const string paramsfile( argv[2] );
				display(paramsfile);
				return 0;
			}
			else{
				cout<<"usage:\n";
				cout<<"hmm --display <paramsfile>\n";
				return 0;
			}
		}
		
		if(string(argv[1]) == string("--predict")){
			if(argc == 5){
				const string paramsfile(argv[2]) , datafile(argv[3]) , outfile(argv[4]);
				prediction(paramsfile,datafile,outfile);
				return 0;
			}
			else{
				cout<<"usage:\n";
				cout<<"hmm --predict <paramsfile> <datafile> <outfile>\n";
				return 0;
			}
		}
		
		{
			cout<<"usage:\n";
			cout<<"hmm --learn <datafile> <Nstates> <paramsfile>\n";
			cout<<"hmm --display <paramsfile>\n";
			cout<<"hmm --predict <paramsfile> <datafile> <outfile>\n";
			return 0;
		}
	}
	else{
		cout<<"usage:\n";
		cout<<"hmm --learn <datafile> <Nstates> <paramsfile>\n";
		cout<<"hmm --display <paramsfile>\n";
		cout<<"hmm --predict <paramsfile> <datafile> <outfile>\n";
		return 0;
	}
	
	return 0;
}








