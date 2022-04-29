//g++ Conf_Setup.C -std=c++0x -o SetupConfFile

#include "CalibrationParameters.h"
#include <stdio.h>
#include <string.h>
#include <fstream>

using namespace std;

bool CreateConfFile(const char * out = "Conf_File.txt",
  const char * inp = "NULL") {

  int ring[24] = {7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,23,22,21,20,19,18,17,16};


  double gains[960], offsets[960];

  const char * testval = "NULL";
  if (strcmp(inp, testval) != 0) {
    ifstream segpar;
    segpar.open(inp);
    if (segpar.is_open()) {
      printf("Segment parameter file: %s opened!\n", inp);
      int j = 0;
      while (!segpar.eof() && j < 960) {
        segpar >> offsets[j] >> gains[j];
        j++;
      }
    } else {
      printf("%s not opened\n", inp);
      return false;
    }
    printf("Segment parameters read in successfully\n");
  } else {
    printf("No segment parameter file declared, setting gains, offsets and nonlinear components to zero\n");
    for (int i = 0; i < 960; i++) {
      gains[i] = 1;
      offsets[i] = 0;
    }
  }

// Uncomment these lines to get a good EmptyGains.cal
/*  for (int i = 0; i < 112; i++){
	bam_gains[i] = 1;
	bam_offsets[i] = 0;
  }
   for (int i = 0; i < 64; i++){
	gain[i] = 1;
	offset[i] = 0;
	non_lin[i] = 0;
  }
*/
  printf("Creating conf file: %s\n", out);
  ofstream outfile;
  outfile.open(out);
  if (!outfile.is_open()) {
    printf("Output not opened\n");
    return false;
  }

  bool tig_flag = true;
  printf("Creating TIGNAMES.txt\n");
  ofstream tignames;
  tignames.open("TIGNAMES.txt");
  if (!tignames.is_open()) {
    printf("Failed to create TIGNAMES, will not be filled\n");
    tig_flag = false;
  }

  char line[128];

  char
  var [64];

  printf("Creating data\n");

  for (int i = 0; i < 1200; i++) {

    int slave = (i / 120) + 1;
    int port = (i % 120) / 10 + 1;
    int channel = i % 10;

    int dnum = i / 60;

    char electronicaddress[32];

 if(i>479 && i<544){
      	      slave = 5;
      port = (i - 480) / 10 + 1;
      channel = (i - 480) % 10;
      sprintf(electronicaddress, "0x%02x0%02x%02x", slave, port, channel);
	          if (i < 512) {
		     sprintf(var, "BAE02EP%02dx", i - 480);
		     outfile << i << "\t" << electronicaddress << "\t" <<
		     var << "\t" << bam_gains[i - 480] << "\t" << bam_offsets[i - 480] << "\t0" << "\tTig10\n";
	          }
		  else if (i > 519) {
        	     sprintf(var, "BAE02EN%02dx", ring[i - 520]);
                     outfile << i << "\t" << electronicaddress << "\t" <<
                     var << "\t" << bam_gains[i - 488] << "\t" << bam_offsets[i - 488] << "\t0" << "\tTig10\n";
      		  }
		  else {
		        int ringval = ring[i+1000];
	          	sprintf(var,"Empty");
	          	outfile << i << "\t" << electronicaddress << "\t" << var << "\t" << 0 << "\t" << 0 << "\t0" << "\tTig10\n"; 
		    }
		  }

  else if (i < 960) {
      int det_pos = i % 60;

      int det_num = (i / 60) + 1;

      sprintf(electronicaddress, "0x%02x0%02x%02x", slave, port, channel);

      if (false) {
        if (det_pos == 0) {
          sprintf(var, "RFL00XSx");
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t0\t0\t0\tTig10\n";
        } else {
          sprintf(var, "Empty");
          outfile << i << "\t0xffffffff\t" <<
            var << "\t0\t0\t0\tTig10\n";
          sprintf(electronicaddress, "0xffffffff", slave, port, channel);
        }
      }
      if (det_num == 9 || det_num == 10) {
        sprintf(var, "Empty");
        outfile << i << "\t0xffffffff\t" <<
          var << "\t0\t0\t0\tTig10\n";
	continue;
      }
	else {

	if(det_num==6)det_num=9;
        int crys_num = (det_pos % 15) + 1;

        bool core = false;
        bool segment = false;

        char colour[1];

        if (det_pos < 15) sprintf(colour, "B");
        else if (det_pos > 14 && det_pos < 30) sprintf(colour, "G");
        else if (det_pos > 29 && det_pos < 45) sprintf(colour, "R");
        else if (det_pos > 44 && det_pos < 60) sprintf(colour, "W");

        if (strcmp(colour, "B") == 0 || strcmp(colour, "R") == 0) {
          if (crys_num < 11) {
            if (det_num < 10) sprintf(var, "TIG%i%i%s", 0, det_num, colour);
            else sprintf(var, "TIG%i%s", det_num, colour);
          } else if (crys_num > 10) {
            if (det_num < 10) sprintf(var, "TIS%i%i%sN", 0, det_num, colour);
            else sprintf(var, "TIS%i%sN", det_num, colour);
          }

          if (crys_num == 1) {
            strcat(var, "N00A");
            core = true;
          } else if (crys_num == 10) strcat(var, "N00B");
          else if (crys_num > 1 && crys_num < 10) {
            char append[10];
            sprintf(append, "P0%ix", crys_num - 1);
            strcat(var, append);
            segment = true;
          } else if (crys_num > 10) {
            char append[10];
            sprintf(append, "0%ix", crys_num - 10);
            strcat(var, append);
          }
        } else if (strcmp(colour, "G") == 0 || strcmp(colour, "W") == 0) {
          if (crys_num < 6) {
            if (det_num < 10) sprintf(var, "TIS%i%i%sN0%ix", 0, det_num, colour, crys_num);
            else sprintf(var, "TIS%i%sN0%ix", det_num, colour, crys_num);
          }
          if (crys_num > 5) {
            if (det_num < 10) sprintf(var, "TIG%i%i%s", 0, det_num, colour);
            else sprintf(var, "TIG%i%s", det_num, colour);
          }
          if (crys_num == 6) {
            strcat(var, "N00A");
            core = true;
          } else if (crys_num == 15) strcat(var, "N00B");
          else if (crys_num > 6 && crys_num < 15) {
            char append[10];
            sprintf(append, "P0%ix", crys_num - 6);
            strcat(var, append);
            segment = true;
          }
        }

        int segment_id;
        if (segment && strcmp(colour, "B") == 0)
          segment_id = (det_num - 1) * 32 + crys_num - 2;
        else if (segment && strcmp(colour, "G") == 0)
          segment_id = (det_num - 1) * 32 + crys_num + 1;
        else if (segment && strcmp(colour, "R") == 0)
          segment_id = (det_num - 1) * 32 + crys_num + 14;
        else if (segment && strcmp(colour, "W") == 0)
          segment_id = (det_num - 1) * 32 + crys_num + 17;

        if (core && strcmp(colour, "B") == 0) {
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t" << gain[(det_num - 1) * 4] << "\t" << offset[(det_num - 1) * 4] << "\t" << non_lin[(det_num - 1) * 4] << "\tTig10\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/g[%i]\" '%f'", i, gain[(det_num - 1) * 4]);
          tignames << line << "\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/o[%i]\" '%f'", i, offset[(det_num - 1) * 4]);
          tignames << line << "\n";
        } else if (core && strcmp(colour, "G") == 0) {
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t" << gain[(det_num - 1) * 4 + 1] << "\t" << offset[(det_num - 1) * 4 + 1] << "\t" << non_lin[(det_num - 1) * 4 + 1] << "\tTig10\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/g[%i]\" '%f'", i, gain[(det_num - 1) * 4 + 1]);
          tignames << line << "\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/o[%i]\" '%f'", i, offset[(det_num - 1) * 4 + 1]);
          tignames << line << "\n";
        } else if (core && strcmp(colour, "R") == 0) {
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t" << gain[(det_num - 1) * 4 + 2] << "\t" << offset[(det_num - 1) * 4 + 2] << "\t" << non_lin[(det_num - 1) * 4 + 2] << "\tTig10\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/g[%i]\" '%f'", i, gain[(det_num - 1) * 4 + 2]);
          tignames << line << "\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/o[%i]\" '%f'", i, offset[(det_num - 1) * 4 + 2]);
          tignames << line << "\n";
        } else if (core && strcmp(colour, "W") == 0) {
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t" << gain[(det_num - 1) * 4 + 3] << "\t" << offset[(det_num - 1) * 4 + 3] << "\t" << non_lin[(det_num - 1) * 4 + 3] << "\tTig10\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/g[%i]\" '%f'", i, gain[(det_num - 1) * 4 + 3]);
          tignames << line << "\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/o[%i]\" '%f'", i, offset[(det_num - 1) * 4 + 3]);
          tignames << line << "\n";
        } else if (segment) {
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t" << gains[i] << "\t" << offsets[i] << "\t" << 0 << "\tTig10\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/g[%i]\" '%f'", i, segment_gains[segment_id]);
          tignames << line << "\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/o[%i]\" '%f'", i, segment_offsets[segment_id]);
          tignames << line << "\n";
        } else {
          outfile << i << "\t" << electronicaddress << "\t" <<
            var << "\t0\t0\t0\tTig10\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/g[%i]\" '%f'", i, 0);
          tignames << line << "\n";
          sprintf(line, "set \"/Analyzer/Shared Parameters/Config/o[%i]\" '%f'", i, 0);
          tignames << line << "\n";
        }

      }

    } 
 else if(i>959 && i<1024){
	      slave = 9;
      port = (i - 960) / 10 + 1;
      channel = (i - 960) % 10;
      sprintf(electronicaddress, "0x%02x0%02x%02x", slave, port, channel);
		   if(i>959 && i<976){
	          	sprintf(var,"BAE01EP%02dx",i-944);
	          	outfile << i << "\t" << electronicaddress << "\t" << var << "\t" << bam_gains[i-904] << "\t" << bam_offsets[i-904] << "\t0" << "\tTig10\n";  		  
		    }
		    else if(i>975 && i<992){
	          	sprintf(var,"BAE01EP%02dx",i-976);
	          	outfile << i << "\t" << electronicaddress << "\t" << var << "\t" << bam_gains[i-904] << "\t" << bam_offsets[i-904] << "\t0" << "\tTig10\n";    		  
		    }
		    else if(i>999){
	          	sprintf(var,"BAE01EN%02dx",ring[i-1000]);
	          	outfile << i << "\t" << electronicaddress << "\t" << var << "\t" << bam_gains[i-912] << "\t" << bam_offsets[i-912] << "\t0" << "\tTig10\n";    		  
		    }
		    else {
			sprintf(var, "Empty");
        		outfile << i << "\t0xffffffff\t" << var << "\t0\t0\t0\tTig10\n";
		    }		
		  }


    else {
      sprintf(var, "Empty");
      outfile << i << "\t0xffffffff\t" <<
        var << "\t0\t0\t0\tTig10\n";
      sprintf(electronicaddress, "0xffffffff", slave, port, channel);
      //sprintf(line,"set \"/Analyzer/Shared Parameters/Config/Name[%i]\" '%s'",i,var);	
      //tignames << line << "\n";
      //sprintf(line,"set \"/Analyzer/Shared Parameters/Config/FSCP[%i]\" '%s'",i,electronicaddress);
      //tignames << line << "\n";
    }

    if (tig_flag) {
      sprintf(line, "set \"/Analyzer/Shared Parameters/Config/Name[%i]\" '%s'", i,
        var);
      tignames << line << "\n";
      sprintf(line, "set \"/Analyzer/Shared Parameters/Config/FSCP[%i]\" '%s'", i, electronicaddress);
      tignames << line << "\n";
    }
  }

  outfile.close();
  tignames.close();

  return true;

}

int main(int argc, char * * argv) {

  bool success = false;
  if (argc == 1) {
    success = CreateConfFile();
  }
  if (argc > 3) {
    printf("Too many inputs, max two\n");
    return 0;
  }

  if (argc == 2)
    success = CreateConfFile(argv[1]);

  if (argc == 3)
    success = CreateConfFile(argv[1], argv[2]);

  if (success) {
    printf("Config files created successfully!\n");
    return 1;
  } else {
    printf("Config creation failed!\n");
    return 0;
  }

  return 0;

}
