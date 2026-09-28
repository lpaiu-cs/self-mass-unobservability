// Request 16, GPL-3.0-or-later. Exercises the actual shared native function.
#include "Delay_brut.h"
#include "Constants.h"
#include "Spline.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
using std::vector;
int main(int argc,char**argv){
 unsetenv("TIMING16_EXTRA_SCALE");unsetenv("TIMING16_EXPORT");
 constexpr long nt=5;value_type t[nt]={-2,-1,0,1,2};
 value_type p[nt][6]={},i[nt][6]={},o[nt][6]={};value_type *pp[nt],*pi[nt],*po[nt];
 vector<vector<value_type>> states(nt,vector<value_type>(30,0));
 for(int k=0;k<nt;++k){pp[k]=p[k];pi[k]=i[k];po[k]=o[k];i[k][0]=10;o[k][0]=20;states[k][9]=3;states[k][10]=4;states[k][12]=-4;states[k][13]=3;}
 value_type n[3]={0,0,1},m[2]={.2L,.3L},out[nt],mean=0;bool bad=false;
 auto evaluate=[&](bool ein,bool shap){Delays_Brut_nogeometric(t,0,nt,2,pp,pi,po,n,1,0,0,1,ein,shap,false,out,1,mean,bad,NULL,&states,m,2,1);};
 evaluate(true,false);assert(!bad);
 const value_type u=(m[0]+m[1])/5*GMsol/clight2;
 for(int k=0;k<nt;++k)assert(fabsl(out[k]-t[k]*u)<1e-18L*(1+fabsl(t[k]*u)));
 evaluate(false,true);assert(!bad);value_type reference=out[0];
 const value_type s=-2.L*4.92521372097374e-6L/daysec*(m[0]+m[1])*logl(5/clight);
 assert(fabsl(reference-s)<1e-26L);
 // Rotate every position and the line of sight; then reorder both extra bodies.
 for(auto& row:states){std::swap(row[10],row[11]);row[10]=-row[10];std::swap(row[13],row[14]);row[13]=-row[13];}
 std::swap(n[1],n[2]);n[1]=-n[1];
 evaluate(false,true);assert(!bad&&fabsl(out[0]-reference)<1e-26L);
 for(auto& row:states)for(int k=0;k<3;++k)std::swap(row[9+k],row[12+k]);std::swap(m[0],m[1]);
 evaluate(false,true);assert(!bad&&fabsl(out[0]-reference)<1e-26L);
 m[0]=m[1]=0;for(auto&row:states)for(int k=9;k<15;++k)row[k]=0;
 evaluate(true,true);assert(!bad);for(auto x:out)assert(x==0);
 m[0]=1;evaluate(true,false);assert(bad); // positive mass at zero separation
 m[0]=-1;evaluate(false,true);assert(bad);
 n[1]=0;n[2]=1;
 m[0]=1;for(auto&row:states)row[11]=-1; // d points along n: invalid log argument
 evaluate(false,true);assert(bad);
 for(const char* invalid:{"nan","-1","1junk",""}){setenv("TIMING16_EXTRA_SCALE",invalid,1);evaluate(false,true);assert(bad);}
 unsetenv("TIMING16_EXTRA_SCALE");states[0].resize(29);evaluate(true,false);assert(bad);
 if(argc==2){
  // Nonuniform quartic data with the native natural-boundary spline.
  value_type x[5]={0,.125L,.5L,1.25L,2},y[5];for(int k=0;k<5;++k)y[k]=powl(x[k],4);
  Spline spl(x,y,5);FILE* fp=fopen(argv[1],"w");assert(fp);
  for(int k=0;k<5;++k)fprintf(fp,"%La %La %La\n",x[k],y[k],spl.y2[k]);
  for(int k=0;k<4;++k){value_type z=(x[k]+x[k+1])/2;fprintf(fp,"%La %La %La\n",z,spl(z,k,k+1),spl.Integrate(x[k],x[k+1],k,k+1,k,k+1));}
  fclose(fp);
 }
 std::cout<<"PASS: shared native Einstein/Shapiro exact controls, zero mass, sum, rotation, permutation, domain rejection\n";
}
