// Directed CAPD/MPFR evaluation of finite momentum intervals outside the cusp.
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
using Six=std::array<MI,6>;

std::string endpoint(const MF& x,bool up){std::ostringstream s;x.put(s,up?MF::RoundUp:MF::RoundDown,128,10,60);return s.str();}
MI encoded(const std::string& a,const std::string& b){return MI(MF(a,MF::RoundDown),MF(b,MF::RoundUp));}
MI binary(const std::string& s){return MI(MF(std::stod(s)));}
void exact(const MI& x){if(x.leftBound()!=x.rightBound())throw std::runtime_error("nonexact partition arithmetic");}
MI power(MI x,int n){MI r(1);for(int j=0;j<n;++j)r*=x;return r;}
MI symmetric(const MI& x){return MI(MF(-1),MF(1))*MI(x.rightBound());}
struct Panel{MI a,b;int depth;};

int main(int argc,char**argv){try{
 if(argc!=6)throw std::runtime_error("args: inputs rule output first count");
 MF::setDefaultPrecision(256);int first=std::stoi(argv[4]),count=std::stoi(argv[5]);
 if(first<0||count<1)throw std::runtime_error("invalid selection");
 std::ifstream rule(argv[2]);std::string lo,hi;std::vector<MI> nodes,weights;
 for(int k=0;k<16;++k){std::string wl,wh;if(!(rule>>lo>>hi>>wl>>wh))throw std::runtime_error("rule input");nodes.push_back(encoded(lo,hi));weights.push_back(encoded(wl,wh));}
 if(!(rule>>lo>>hi))throw std::runtime_error("norm input");MI norm=encoded(lo,hi);
 for(int k=0;k<32;++k){MI sum(0),truth=MI(1)/(k+1);for(int n=0;n<16;++n)sum+=weights[n]*power(nodes[n],k);
  if(sum.leftBound()>truth.leftBound()||sum.rightBound()<truth.rightBound())throw std::runtime_error("Gauss interval moment");}
 std::ifstream in(argv[1]);std::ofstream out(argv[3]);if(!in||!out)throw std::runtime_error("file open");
 int index,cell,j,done=0;std::string bs,qs,ss,els,ehs,hs,ps;
 while(in>>index>>cell>>j>>bs>>qs>>ss>>els>>ehs>>hs>>ps){
  if(index<first||index>=first+count)continue;
  MI beta=binary(bs),Q=encoded(qs,qs),Sref=binary(ss),eta=encoded(els,ehs),h_cusp=binary(hs),pcut=binary(ps);
  if(beta.leftBound()<=MF(0)||Q.leftBound()<=MF(0)||Sref.leftBound()<=MF(0)||pcut.leftBound()<=MF(0)||h_cusp.leftBound()<=MF(0)||h_cusp.rightBound()>=MF(1))throw std::runtime_error("invalid physical input");
  MI left=Q*(1-h_cusp),right=Q*(1+h_cusp);exact(left);exact(right);exact(pcut);
  std::vector<Panel> stack;MI length(0),covered(0);int panels=0,omitted=0,maxdepth=0;MI phase_max(0);
  if(pcut.leftBound()>right.rightBound()){stack.push_back({right,pcut,0});length+=pcut-right;}
  MI end=pcut.rightBound()<left.leftBound()?pcut:left;stack.push_back({MI(0),end,0});length+=end;
  exact(length);Six total{},remainders{};auto start=std::chrono::steady_clock::now();
  while(!stack.empty()){
   Panel panel=stack.back();stack.pop_back();MI a=panel.a,b=panel.b,c=(a+b)/2,h=(b-a)/2,R=4*h,width=b-a;
   exact(c);exact(h);exact(R);if(panel.depth>48)throw std::runtime_error("depth limit");
   if(panel.depth>maxdepth)maxdepth=panel.depth;
   auto split=[&](){stack.push_back({c,b,panel.depth+1});stack.push_back({a,c,panel.depth+1});};
   MI distance=abs(c-Q);if(R.rightBound()>=distance.leftBound()||R.rightBound()>=(c+Q).leftBound()){split();continue;}
   MI lower=c-R;if(lower.leftBound()<MF(0))lower=MI(0);
   MI g2=1+lower*lower-R*R;if(g2.leftBound()<=MF(0)){split();continue;}
   MI g=sqrt(g2),pmax=c+R,phase=pmax*R/(beta*g);if(phase.rightBound()>MF(1)){split();continue;}
   MI tmin=(g2-1)/(beta*(g+1)),tmax=pmax*pmax/(beta*(g+1)),exponent=MI(eta.rightBound())-tmin,qmax(1);
   if(exponent.rightBound()<MF(0)){MF e=exponent.rightBound();if(e<MF(-128))e=MF(-128);qmax=exp(MI(e));}
   MI logs=abs(log(c+Q))-log(1-R/(c+Q))+abs(log(distance))-log(1-R/distance);
   MI A=pmax*pmax/(g*Sref),B=pmax*(abs(1+c*c-Q*Q)+R*(2*c+R))/(2*Q*g*Sref);
   MI base=(A+B*logs)*qmax;Six factors={MI(1),MI(1),MI(1),tmax,tmax,tmax*tmax+tmax};
   Six envelope{},error{};MI allowance=width/pcut*encoded("1e-18","1e-18");bool skip=true,accept=true;
   for(int k=0;k<6;++k){envelope[k]=base*factors[k];error[k]=width*norm*power(MI(2)/3,32)*envelope[k];
    if((width*envelope[k]).rightBound()>allowance.leftBound())skip=false;
    if(error[k].rightBound()>allowance.leftBound())accept=false;}
   if(!skip&&!accept){split();continue;}
   if(phase.rightBound()>phase_max.rightBound())phase_max=MI(phase.rightBound());
   if(skip){++omitted;for(int k=0;k<6;++k){MI bound=width*envelope[k];total[k]+=symmetric(bound);remainders[k]+=MI(bound.rightBound());}}
   else{
    Six quadrature{};
    for(int n=0;n<16;++n){MI p=a+width*nodes[n],gamma=sqrt(1+p*p),t=p*p/(beta*(gamma+1)),q=1/(1+exp(t-eta)),v=1-q;
     MI kernel=(p*p/gamma+p*(1+p*p-Q*Q)/(2*Q*gamma)*(log(p+Q)-log(abs(p-Q))))/Sref;
     Six f={q,q*v,q*v*(1-2*q),t*q*v,t*q*v*(1-2*q),(t*t*(1-2*q)-t)*q*v};
     for(int k=0;k<6;++k)quadrature[k]+=width*weights[n]*kernel*f[k];}
    for(int k=0;k<6;++k){total[k]+=quadrature[k]+symmetric(error[k]);remainders[k]+=MI(error[k].rightBound());}
   }
   ++panels;covered+=width;
  }
  exact(covered);if(covered.leftBound()!=length.leftBound())throw std::runtime_error("partition coverage");
  out<<"{\"case\":"<<index<<",\"cell\":"<<cell<<",\"z_index\":"<<j<<",\"panels\":"<<panels<<",\"bounded_omissions\":"<<omitted<<",\"maximum_depth\":"<<maxdepth
     <<",\"phase_upper\":\""<<endpoint(phase_max.rightBound(),true)<<"\",\"covered_length\":\""<<endpoint(covered.rightBound(),true)<<"\",\"enclosures\":[";
  for(int k=0;k<6;++k){if(k)out<<',';out<<"[\""<<endpoint(total[k].leftBound(),false)<<"\",\""<<endpoint(total[k].rightBound(),true)<<"\"]";}
  out<<"],\"remainder_upper\":[";for(int k=0;k<6;++k){if(k)out<<',';out<<'"'<<endpoint(remainders[k].rightBound(),true)<<'"';}
  out<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
  if(!out)throw std::runtime_error("output failure");++done;
 }
 if(done!=count)throw std::runtime_error("incomplete case coverage");std::cout<<"PASS finite momentum cases "<<done<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
