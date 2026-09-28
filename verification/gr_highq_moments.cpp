// Directed truncated thermal moment jets for the positive high-Q expansion.
#include "capd/intervals/MpInterval.h"
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using MI=capd::intervals::MpInterval;
using MF=capd::multiPrec::MpReal;
using capd::abs;
constexpr int MOMENTS=17,FIELDS=6,SIZE=MOMENTS*FIELDS;
using Jet=std::array<MI,SIZE>;
MI encoded(const std::string&a,const std::string&b){return MI(MF(a,MF::RoundDown),MF(b,MF::RoundUp));}
MI binary(const std::string&s){return MI(MF(std::stod(s)));}
MI power(MI x,int n){MI result(1);for(int k=0;k<n;++k)result*=x;return result;}
void exact(const MI&x){if(x.leftBound()!=x.rightBound())throw std::runtime_error("nonexact partition");}
MI symmetric(const MI&x){return MI(MF(-1),MF(1))*MI(x.rightBound());}
std::string endpoint(const MF&x,bool up){std::ostringstream s;x.put(s,up?MF::RoundUp:MF::RoundDown,128,10,60);return s.str();}
struct Panel{MI a,b;int depth;};

int main(int argc,char**argv){try{
 if(argc!=6)throw std::runtime_error("args: input rule output first count");
 MF::setDefaultPrecision(128);int first=std::stoi(argv[4]),count=std::stoi(argv[5]);if(first<0||count<1)throw std::runtime_error("invalid selection");
 std::ifstream rule(argv[2]);std::string lo,hi;std::vector<MI> nodes,weights;
 for(int i=0;i<16;++i){std::string wl,wh;if(!(rule>>lo>>hi>>wl>>wh))throw std::runtime_error("rule input");nodes.push_back(encoded(lo,hi));weights.push_back(encoded(wl,wh));}
 if(!(rule>>lo>>hi))throw std::runtime_error("rule norm");MI norm=encoded(lo,hi),error_factor=norm*power(MI(2)/3,32);
 for(int k=0;k<32;++k){MI sum(0),target=MI(1)/(k+1);for(int n=0;n<16;++n)sum+=weights[n]*power(nodes[n],k);
  if(sum.leftBound()>target.leftBound()||sum.rightBound()<target.rightBound())throw std::runtime_error("Gauss moment control");}
 std::ifstream in(argv[1]);std::ofstream out(argv[3]);if(!in||!out)throw std::runtime_error("file open");
 int index,cell,done=0;std::string bs,ss,el,eh,ps;
 while(in>>index>>cell>>bs>>ss>>el>>eh>>ps){
  if(index<first||index>=first+count)continue;
  MI beta=binary(bs),Sref=binary(ss),eta=encoded(el,eh),pcut=binary(ps);exact(beta);exact(Sref);exact(pcut);
  if(beta.leftBound()<=MF(0)||Sref.leftBound()<=MF(0)||pcut.leftBound()<=MF(0))throw std::runtime_error("invalid input");
  std::vector<Panel> stack={{MI(0),pcut,0}};Jet total{},errors{};MI covered(0);int panels=0,omitted=0,maxdepth=0;auto start=std::chrono::steady_clock::now();
  while(!stack.empty()){
   Panel panel=stack.back();stack.pop_back();MI a=panel.a,b=panel.b,c=(a+b)/2,h=(b-a)/2,R=4*h,width=b-a;exact(c);exact(h);exact(R);
   if(panel.depth>48)throw std::runtime_error("depth limit");if(panel.depth>maxdepth)maxdepth=panel.depth;
   auto split=[&](){stack.push_back({c,b,panel.depth+1});stack.push_back({a,c,panel.depth+1});};
   MI lower=c-R;if(lower.leftBound()<MF(0))lower=MI(0);MI g2=1+lower*lower-R*R;if(g2.leftBound()<=MF(0)){split();continue;}
   MI g=sqrt(g2),pmax=c+R,phase=pmax*R/(beta*g);if(phase.rightBound()>MF(1)){split();continue;}
   MI tmin=(g2-1)/(beta*(g+1)),tmax=pmax*pmax/(beta*(g+1)),exponent=MI(eta.rightBound())-tmin,qmax(1);
   if(exponent.rightBound()<MF(0)){MF e=exponent.rightBound();if(e<MF(-128))e=MF(-128);qmax=exp(MI(e));}
   std::array<MI,FIELDS> factors={MI(1),MI(1),MI(1),tmax,tmax,tmax*tmax+tmax};
   MI base=qmax/(g*Sref),ratio=pmax*pmax/(pcut*pcut),monomial=ratio,allowance=width/pcut*encoded("1e-13","1e-13");Jet envelope{},remainder{};bool skip=true,accept=true;
   for(int m=0;m<MOMENTS;++m){for(int k=0;k<FIELDS;++k){int j=m*FIELDS+k;envelope[j]=base*monomial*factors[k];remainder[j]=width*error_factor*envelope[j];
     if((width*envelope[j]).rightBound()>allowance.leftBound())skip=false;if(remainder[j].rightBound()>allowance.leftBound())accept=false;}monomial*=ratio;}
   if(!skip&&!accept){split();continue;}
   if(skip){++omitted;for(int j=0;j<SIZE;++j){MI bound=width*envelope[j];total[j]+=symmetric(bound);errors[j]+=MI(bound.rightBound());}}
   else{
    Jet quadrature{};
    for(int n=0;n<16;++n){MI p=a+width*nodes[n],gamma=sqrt(1+p*p),t=p*p/(beta*(gamma+1)),q=1/(1+exp(t-eta)),v=1-q;
     std::array<MI,FIELDS> f={q,q*v,q*v*(1-2*q),t*q*v,t*q*v*(1-2*q),(t*t*(1-2*q)-t)*q*v};
     MI ratio=p*p/(pcut*pcut),monomial=ratio,weight=width*weights[n]/(gamma*Sref);
     for(int m=0;m<MOMENTS;++m){MI pref=weight*monomial;for(int k=0;k<FIELDS;++k)quadrature[m*FIELDS+k]+=pref*f[k];monomial*=ratio;}
    }
    for(int j=0;j<SIZE;++j){total[j]+=quadrature[j]+symmetric(remainder[j]);errors[j]+=MI(remainder[j].rightBound());}
   }
   covered+=width;++panels;
  }
  exact(covered);if(covered.leftBound()!=pcut.leftBound())throw std::runtime_error("coverage");
  out<<"{\"position\":"<<index<<",\"cell\":"<<cell<<",\"panels\":"<<panels<<",\"bounded_omissions\":"<<omitted<<",\"maximum_depth\":"<<maxdepth<<",\"moments\":[";
  for(int m=0;m<MOMENTS;++m){if(m)out<<',';out<<'[';for(int k=0;k<FIELDS;++k){if(k)out<<',';int j=m*FIELDS+k;out<<"[\""<<endpoint(total[j].leftBound(),false)<<"\",\""<<endpoint(total[j].rightBound(),true)<<"\"]";}out<<']';}
  out<<"],\"maximum_remainder_upper\":\"";MF maximum(0);for(auto&x:errors)if(x.rightBound()>maximum)maximum=x.rightBound();out<<endpoint(maximum,true)<<"\",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
  if(!out)throw std::runtime_error("output failure");++done;
 }
 if(done!=count)throw std::runtime_error("incomplete coverage");std::cout<<"PASS moment states "<<done<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
