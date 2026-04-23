// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include "Constants.h"
#include"Utilities.h"
#include <map>


using namespace std;

void Read_parfile(const char * filename, value_type parameters[18], value_type parameter_shift_scales[18],
                  value_type& treference, int& integrator, int& interpsteps_per_period_i,
                  bool& roemer, bool& einstein, bool& shapiro, bool& aberration) {
    // old version : now use the one in parametres.cpp
    FILE *fp;
    const int sizechar = 100;
    //int returnvalue = 0;
    value_type val, dval = zero;
    int ival ;
    char line[sizechar];
    char textin[sizechar]; // Should be the same size as line otherwise sscanf fails
    string str ;

    fp = fopen(filename,"r");

    while ( fgets(line, sizechar, fp) != NULL ){
        val = zero ;
        dval = zero;
            //returnvalue = fscanf(fp, " \n" ); // Ne fonctionne pas..
      //  cout << "line " << line << endl; Test
        if (sscanf( line, "%s ", textin) == EOF )
        {
      //      cout << "skip line " << line << endl; Test
            continue ;
    }

        str = textin ;
       // cout << " str " << str << endl; Test

        if (str.find("spinfreq1") != string::npos) {
            sscanf( line, "%s %Le %Le ", textin, &val, &dval );
            parameters[1] = val ;
            parameter_shift_scales[1] = dval;
        }
        else if (str.find("spinfreq") != string::npos){ // ! Be careful with the order !
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[0] = val;
              //  cout << "lfskhglo h  " << parameters[0] << "  val " << val <<   endl; Test
                parameter_shift_scales[0] = dval;
        }
        else if (str.find("eta_p") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[2] = val;
                parameter_shift_scales[2] = dval;
        }
        else if (str.find("apsini_i") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[3] = val;
                parameter_shift_scales[3] = dval;
        }
        else if (str.find("kappa_p") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[4] = val;
                parameter_shift_scales[4] = dval;
        }
        else if (str.find("apcosi_i") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[5] = val;
                parameter_shift_scales[5] = dval;
        }
        else if (str.find("tasc_p") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[6] = val;
                parameter_shift_scales[6] = dval;
        }
        else if (str.find("period_i") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[7] = val;
                parameter_shift_scales[7] = dval;
        }
        else if (str.find("eta_b") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[8] = val;
                parameter_shift_scales[8] = dval;
        }
        else if (str.find("absini_o") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[9] = val;
                parameter_shift_scales[9] = dval;
        }
        else if (str.find("kappa_b") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[10] = val;
                parameter_shift_scales[10] = dval;
        }
        else if (str.find("abcosi_o") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[11] = val;
                parameter_shift_scales[11] = dval;
        }
        else if (str.find("tasc_b") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[12] = val;
                parameter_shift_scales[12] = dval;
        }
        else if (str.find("period_o") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[13] = val;
                parameter_shift_scales[13] = dval;
        }
        else if (str.find("masspar_p") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[14] = val;
                parameter_shift_scales[14] = dval;
        }
        else if (str.find("masspar_i") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[15] = val;
                parameter_shift_scales[15] = dval;
        }
        else if (str.find("mass_o") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[16] = val;
                parameter_shift_scales[16] = dval;
        }
        else if (str.find("deltaoman") != string::npos){
                sscanf( line, "%s %Le %Le ", textin, &val, &dval );
                parameters[17] = val;
                parameter_shift_scales[17] = dval;
        }
        else if (str.find("treference") != string::npos){
            sscanf( line, "%s %Le ", textin, &val );
            treference = val;
        }
        else if (str.find("roemer") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            roemer = ival;
        }
        else if (str.find("einstein") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            einstein = ival;
        }
        else if (str.find("shapiro") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            shapiro = ival;
        }
        else if (str.find("aberration") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            aberration = ival;
        }
        else if (str.find("interpsteps") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            interpsteps_per_period_i = ival;
        }
        else if (str.find("integrator") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            integrator = ival;
        }
        else
        {}
    }
    fclose(fp);
    return;
}



void Read_TNS_file(char * filename, value_type  * &toas, turntype * &turns, errortype * &errors, long int& ntoas ){
    long int i = 0L ;
    long int j,k =0L;
    FILE * fp ;
    const int sizechar = 150;
    char line[sizechar];
    value_type toa = zero;
    turntype turn = 0ll;
    errortype error = 0.;
    vector<long int> commentedlines;
    commentedlines.reserve(10);

    long int number_of_lines = 0L;
    int ch = 0;

    // Count the lines

    fp = fopen(filename, "r");

    ntoas = 0L;
    do
    {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            number_of_lines++;
            ntoas++;
        }
        if (ch == '#')
        {
            commentedlines.push_back(number_of_lines);
            ntoas--;
        }
    } while (ch != EOF);

    // last line doesn't end with a new line!
    // but there has to be a line at least before the last line
    if(ch != '\n' && number_of_lines != 0)
    {
        number_of_lines++;
        ntoas++;
    }

    // <----------finished counting the lines

    commentedlines.push_back(number_of_lines + 1); // TO deal with the last value of j below
    rewind(fp);



// Table allocation :
    toas = new value_type[ntoas];
    turns = new turntype[ntoas];
    errors = new errortype[ntoas];
    j =0L;
    k = 0L;
    // Is there a hidden empty line ? line number = 'Warning  : Problem in line counting in "Read_TNS_file" ! '
    for ( i = 0 ;  i < number_of_lines ; ++i){
        if (NULL == fgets(line, sizechar,fp) ) {
            cout << endl << "Warning  : Problem in line counting in Read_TNS_file ! Is there a hidden empty line ? "<< endl ;
            cout << "Line number =  " <<  i << endl << endl;
            ntoas -= 1 ;
        }
       if (i != commentedlines[j] )
       {
            if (EOF == sscanf( line, "%Le %lli %f", &toa, &turn, &error ) ) {
                cout << endl << "Warning  : Problem in line counting in Read_TNS_file ! " << i << endl << endl;
            }

        toas[k] = toa;
        turns[k] = turn ;
        errors[k] = error;
        ++k;
       }
       else
           ++j;
    }

    fclose(fp);

    return ;
}




void Write_new_tim_file(const char * datafile_in, const char * datafile_out, const value_type * SATs)
// Rewrite datafile_in .tim file replacing the sats by the value given in ''SATs'' in datafile_out. Used for creating fake tim files.
{
    FILE *fp;
    FILE *fout;
    const int sizechar = 1000;
    char line[sizechar];
    char textin[sizechar]; // Should be the same size as line otherwise sscanf fails
    char textin2[sizechar]; // Should be the same size as line otherwise sscanf fails
    

    value_type lecturenldb;
    double lecturedb;
    int isat=0;

    printf("\n\n\n !!!!!!!!!!! DOES NOT WORK IF NOT SORTED TIM FILE IN THE FIRST PLACE !!!!!!!!!!!!!!!!!!!!!!!!!! \n\n\n");


    fp = fopen(datafile_in,"r");
    fout = fopen(datafile_out,"w");

    isat = 0;
    while ( fgets(line, sizechar, fp) != NULL ){


        if (sscanf( line, "%s ", textin) == EOF )  continue ;
        if ( sscanf( line, "%s %lf %Le %[^\n]", textin, &lecturedb, &lecturenldb, textin2 ) < 4 )
        {
            fprintf(fout, "%s", line);
        }
        else
        {
            fprintf(fout, "%s %.8f %.15Lf %s\n", textin, lecturedb, SATs[isat], textin2);
            isat+=1;
        }
    }

    fclose(fp);
    fclose(fout);


    return;

}


void Write_new_tim_file(const char * datafile_in, const char * datafile_out, const value_type * SATs, const double * stddev)
// Rewrite datafile_in .tim file replacing the sats by the value given in ''SATs'' in datafile_out with error given in stddev.
// Used for creating fake tim files.
// In this versions uncertainties can be specified arbitrarily using the array stddev
{
    FILE *fp;
    FILE *fout;
    const int sizechar = 1000;
    char line[sizechar];
    char textin[sizechar]; // Should be the same size as line otherwise sscanf fails
    char textin2[sizechar]; // Should be the same size as line otherwise sscanf fails

    value_type lecturenldb;
    double lecturedb;
    double lecturestd;
    int isat=0;

    printf("\n\n\n !!!!!!!!!!! DOES NOT WORK IF NOT SORTED TIM FILE IN THE FIRST PLACE !!!!!!!!!!!!!!!!!!!!!!!!!! \n\n\n");


    fp = fopen(datafile_in,"r");
    fout = fopen(datafile_out,"w");

    isat = 0;
    while ( fgets(line, sizechar, fp) != NULL ){


        if (sscanf( line, "%s ", textin) == EOF )  continue ;

        if ( sscanf( line, "%s %lf %Le %lf %[^\n]", textin, &lecturedb, &lecturenldb, &lecturestd, textin2) < 4 )
        {
            fprintf(fout, "%s", line);
        }
        else
        {
            fprintf(fout, "%s %.8f %.15Lf %.3f %s\n", textin, lecturedb, SATs[isat], stddev[isat], textin2);
            isat+=1;
        }
    }

    fclose(fp);
    fclose(fout);


    return;
}


void Write_new_tim_file(const char * datafile_in, const char * datafile_out, const value_type * SATs, const double * stddev, 
                        const int * mask, const int warning)
// Rewrite datafile_in .tim file replacing the sats by the value given in ''SATs'' in datafile_out with error given in stddev.
// Used for creating masked tim files.
// In this versions uncertainties can be specified arbitrarily using the array stddev
// In this version, if mask[i] is True, add a comment sign "#" at the beginning of line i
// Using int for mask and warning instead of bool for easier python bindings
// warning : if <= 0 shows no warnings ; 
//           if 1 shows only warnings if datafile_out doesn't match values provided by SATs and stddev ( with 1e-14 and 1e-6 resp)
//           if >=2 shows everything, i.e. 1 and order warning.
{
    FILE *fp;
    FILE *fout;
    const int sizechar = 1000;
    char line[sizechar];
    char textin[sizechar]; // Should be the same size as line otherwise sscanf fails
    char textin2[sizechar]; // Should be the same size as line otherwise sscanf fails

    value_type lecturenldb;
    double lecturedb;
    double lecturestd;
    int scanval;
    int isat=0;
    
    if (warning>1) printf("\n\n\n !!!!!!!!!!! DOES NOT WORK IF NOT SORTED TIM FILE IN THE FIRST PLACE !!!!!!!!!!!!!!!!!!!!!!!!!! \n\n\n");


    fp = fopen(datafile_in,"r");
    fout = fopen(datafile_out,"w");

    isat = 0;
    while ( fgets(line, sizechar, fp) != NULL ){


        if (sscanf( line, "%s ", textin) == EOF )  continue ;
    
        scanval = sscanf( line, "%s %lf %Lf %lf %[^\n]", textin, &lecturedb, &lecturenldb, &lecturestd, textin2 );
        if ( scanval < 4  or textin[0]=='C' or textin[0]=='#')
        {
            fprintf(fout, "%s", line);
        }
        else
        {
            if (warning > 0) // Check for differences between internal values and values from datafile_out
            {
                if (abs(lecturenldb - SATs[isat]) > 1e-14)
                {
                    printf("\n WARNING : In Write_new_tim_file, internal SAT nb %d does not match SAT read from %s. Difference = %Le days. \n", isat, datafile_in, lecturenldb - SATs[isat]);
                }
                if (abs(lecturestd - stddev[isat]) > 1e-6)
                {
                    printf("\n WARNING : In Write_new_tim_file, internal TOA error nb %d does not match TOA error read from %s. Difference = %e microsec. \n", isat, datafile_in, lecturestd - stddev[isat]);
                }
            }
            if (mask[isat] > 0) fprintf(fout, "# ");
            fprintf(fout, "%s %.8f %.15Lf %.6f %s\n", textin, lecturedb, SATs[isat], stddev[isat], textin2);
            isat+=1;
        }
    }

    fclose(fp);
    fclose(fout);


    return;
}

void Sortoutdatafile(const char * datafile, const int maxlinesize)
// Sort out a tim file by toas and record it in datafile+"-sorted".
// maxlinesize  gives the maximum size of a line in datafile
{
    FILE *fp;
  //  const int sizechar = MAX_FILELEN;
    char line[maxlinesize];
    char textin[maxlinesize]; // Should be the same size as line otherwise sscanf fails
    char textin2[maxlinesize]; // Should be the same size as line otherwise sscanf fails
    map<value_type, string> toastosort;
    value_type lecturenldb;
    double lecturedb;
    value_type header = -100000. ;


    fp = fopen(datafile,"r");


    while ( fgets(line, maxlinesize, fp) != NULL ){

        if (sscanf( line, "%s ", textin) == EOF )  continue ;
        if ( sscanf( line, "%s %lf %Le %s", textin, &lecturedb, &lecturenldb, textin2 ) < 4 )
        {
            toastosort[header] = line;
            header += 1.;
        }
        else
        {
            toastosort[lecturenldb] = line;
        }
    }

    fclose(fp);

    strcpy(textin,datafile);
    strcat(textin,"-sorted");

    printf("\n\n *** Writing the sorted tim file in %s ***\n\n", textin);

    fp = fopen(textin, "w");
    for (map<value_type, string>::iterator it=toastosort.begin(); it!=toastosort.end(); ++it)
    {
        fprintf(fp, "%s",it->second.c_str());
    }
    fclose(fp);

    return;

}
