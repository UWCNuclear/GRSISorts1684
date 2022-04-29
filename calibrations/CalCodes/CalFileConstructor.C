// g++ CalFileConstructor.C -std=c++0x -o ConstructCalibrationFile 

#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip> 

using namespace std;

void ConstructCalFile(const char *inp = "Conf_File.txt", const char *out = "CalibrationFile.cal"){

	std::string tempadd, tempmnemonic,tempadd2,tempdig;
	double tempgain,tempoffset,tempnonlin;
	int tempchan;
	
	std::vector<std::string> ADDRESS,MNEMONIC,DIGITIZER;
	std::vector<int> channelnumber;
	std::vector<double> gain,offset,nonlin;

	int j=0;

	ifstream infile;
	infile.open(inp);
	if(infile.is_open()){
		printf("Conf file opened\n");
		while(!infile.eof()){
			infile >> tempchan >> tempadd >> tempmnemonic >> tempgain >> tempoffset >> tempnonlin >> tempdig;
			channelnumber.push_back(tempchan);
			ADDRESS.push_back(tempadd);
			MNEMONIC.push_back(tempmnemonic);
			gain.push_back(tempgain);
			offset.push_back(tempoffset);
			nonlin.push_back(tempnonlin);
			DIGITIZER.push_back(tempdig);
			j++;
		}	
	}
	else{
		printf("Conf file %s failed to open, abort!\n",inp);
		return;	
	}

	printf("Config file successfully read, constructing calibration file...\n");
		
	ofstream outfile(out,std::ofstream::out);
	
	
	for(int i=0;i<channelnumber.size();i++){
	
		if(MNEMONIC.at(i)[0] == 'T' && MNEMONIC.at(i)[1] == 'I'){
	
			outfile << MNEMONIC.at(i) << " " << ADDRESS.at(i) << " { \n";
			outfile << "Name:\t" << MNEMONIC.at(i) << " " << ADDRESS.at(i) << "\n";
			outfile << "Number:\t" << channelnumber.at(i) << "\n";
			outfile << "Address:\t" << ADDRESS.at(i) << "\n";
			outfile << "Digitizer:\t" << DIGITIZER.at(i) << "\n";
			outfile << "EngCoeff:\t" << offset.at(i) << " " << gain.at(i) << " " << nonlin.at(i) << "\n"; 
			outfile << "Integration:\t" << 125 << "\n";
			outfile << "ENGChi2:\t" << "\n";
			outfile << "FileInt:\t" << "true\n";
			outfile << "}\n";
			outfile << "\n";
			outfile << "//====================================//\n";
			
		}
		else if(MNEMONIC.at(i)[0] == 'S' && MNEMONIC.at(i)[1] == 'H'){
	
			outfile << MNEMONIC.at(i) << " " << ADDRESS.at(i) << " { \n";
			outfile << "Name:\t" << MNEMONIC.at(i) << " " << ADDRESS.at(i) << "\n";
			outfile << "Number:\t" << channelnumber.at(i) << "\n";
			outfile << "Address:\t" << ADDRESS.at(i) << "\n";
			outfile << "Digitizer:\t" << DIGITIZER.at(i) << "\n";
			outfile << "EngCoeff:\t" << offset.at(i) << " " << gain.at(i) << " " << nonlin.at(i) << "\n"; 
			outfile << "Integration:\t" << 25 << "\n";
			outfile << "ENGChi2:\t" << "\n";
			outfile << "FileInt:\t" << "true\n";
			outfile << "}\n";
			outfile << "\n";
			outfile << "//====================================//\n";	
			
		}
		else if(MNEMONIC.at(i)[0] == 'S' && MNEMONIC.at(i)[1] == 'P'){
		
			outfile << MNEMONIC.at(i) << " " << ADDRESS.at(i) << " { \n";
			outfile << "Name:\t" << MNEMONIC.at(i) << " " << ADDRESS.at(i) << "\n";
			outfile << "Number:\t" << channelnumber.at(i) << "\n";
			outfile << "Address:\t" << ADDRESS.at(i) << "\n";
			outfile << "Digitizer:\t" << DIGITIZER.at(i) << "\n";
			outfile << "EngCoeff:\t" << offset.at(i) << " " << gain.at(i) << " " << nonlin.at(i) << "\n"; 
			outfile << "Integration:\t" << 125 << "\n";
			outfile << "ENGChi2:\t" << "\n";
			outfile << "FileInt:\t" << "true\n";
			outfile << "}\n";
			outfile << "\n";
			outfile << "//====================================//\n";	
		
		}
		else {
			outfile << MNEMONIC.at(i) << " " << ADDRESS.at(i) << " { \n";
			outfile << "Name:\t" << MNEMONIC.at(i) << " " << ADDRESS.at(i) << "\n";
			outfile << "Number:\t" << channelnumber.at(i) << "\n";
			outfile << "Address:\t" << ADDRESS.at(i) << "\n";
			outfile << "Digitizer:\t" << DIGITIZER.at(i) << "\n";
			outfile << "EngCoeff:\t" << offset.at(i) << " " << gain.at(i) << " " << nonlin.at(i) << "\n"; 
			outfile << "Integration:\t" << 125 << "\n";
			outfile << "ENGChi2:\t" << "\n";
			outfile << "FileInt:\t" << "true\n";
			outfile << "}\n";
			outfile << "\n";
			outfile << "//====================================//\n";	
		}   

	}

	cout << "Construction complete!!!\n";

	outfile.close();
	
}

int main(int argc, char **argv){

	printf("Constructing CalFile\n");
	if(argc==1){
		printf("Assuming default .conf file name: Conf_File.conf\nFormat: Input conf file (default: Conf_File.conf),\tOutput calibration file (default: CalibrationFile.cal)\n");
		ConstructCalFile();
	}
	else if(argc==2)
		ConstructCalFile(argv[1]);
	else if(argc==3)
		ConstructCalFile(argv[1],argv[2]);
	
	return 0;

}
