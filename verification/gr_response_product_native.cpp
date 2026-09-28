// Reusable certified H interpolation and logarithmic product integration.
#include "capd/intervals/MpInterval.h"
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using MI=capd::intervals::MpInterval;using MF=capd::multiPrec::MpReal;using capd::abs;
constexpr int N=32;using Six=std::array<MI,6>;using Coeff=std::array<Six,N>;
MI encoded(const std::string&a,const std::string&b){return MI(MF(a,MF::RoundDown),MF(b,MF::RoundUp));}
MI fraction(const std::string&s){auto p=s.find('/');return p==std::string::npos?encoded(s,s):encoded(s.substr(0,p),s.substr(0,p))/encoded(s.substr(p+1),s.substr(p+1));}
MI binary(const std::string&s){return MI(MF(std::stod(s)));}
MI power(MI x,int n){MI y(1);for(int k=0;k<n;++k)y*=x;return y;}
void exact(const MI&x){if(x.leftBound()!=x.rightBound())throw std::runtime_error("nonexact partition or coordinate");}
MI symmetric(const MI&x){return MI(MF(-1),MF(1))*MI(x.rightBound());}
std::string endpoint(const MF&x,bool up,int digits=60){std::ostringstream s;x.put(s,up?MF::RoundUp:MF::RoundDown,128,10,digits);return s.str();}
std::string point(const MI&x){exact(x);auto s=endpoint(x.leftBound(),false,200);MI y=encoded(s,s);exact(y);if(y.leftBound()!=x.leftBound())throw std::runtime_error("partition serialization");return s;}
void json_six(std::ostream&out,const Six&v){out<<'[';for(int k=0;k<6;++k){if(k)out<<',';out<<"[\""<<endpoint(v[k].leftBound(),false)<<"\",\""<<endpoint(v[k].rightBound(),true)<<"\"]";}out<<']';}
std::vector<MI> log_moments(MI z){
 exact(z);std::vector<MI> J(N+1);MI az=abs(z);
 if(az.leftBound()==MF(1)){J[0]=2*log(MI(2))-2;for(int n=1;n<=N;++n)J[n]=MI(-2*(z.leftBound()<MF(0)&&n%2?-1:1))/(n*(n+1));return J;}
 if(az.leftBound()>=MF(1.25)){
  MI rho=1/(az+sqrt(az*az-1)),r(MF(0),rho.rightBound());std::vector<MI> ratios(N+1),Q(N+2);constexpr int padding=160;
  for(int n=N+1+padding;n>=1;--n){MI den=(2*n+1)*az-(n+1)*r;if(den.leftBound()<=MF(0))throw std::runtime_error("continued fraction denominator");r=MI(n)/den;if(n-1<=N)ratios[n-1]=r;}
  Q[0]=log((az+1)/(az-1))/2;for(int n=0;n<=N;++n)Q[n+1]=Q[n]*ratios[n];
  J[0]=(az+1)*log(az+1)-(az-1)*log(az-1)-2;
  for(int n=1;n<=N;++n)J[n]=MI(2*(z.leftBound()<MF(0)&&n%2?-1:1))*(Q[n+1]-Q[n-1])/(2*n+1);
 }else{
  MI q0=log(abs((z+1)/(z-1)))/2;std::vector<MI>P(N+2),W(N+2);P[0]=1;P[1]=z;W[0]=0;W[1]=1;
  for(int n=1;n<=N;++n){P[n+1]=((2*n+1)*z*P[n]-n*P[n-1])/(n+1);W[n+1]=((2*n+1)*z*W[n]-n*W[n-1])/(n+1);}
  J[0]=(z+1)*log(abs(z+1))-(z-1)*log(abs(z-1))-2;
  for(int n=1;n<=N;++n)J[n]=2*((P[n+1]-P[n-1])*q0-W[n+1]+W[n-1])/(2*n+1);
 }
 return J;
}
struct Panel{MI a,b;int depth;};struct Stored{MI a,b;bool skip;Six error;Coeff c;};

int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("mode and arguments");std::string mode=argv[1];MF::setDefaultPrecision(mode=="build"?128:256);
 if(mode=="moments"){
  if(argc!=3)throw std::runtime_error("moments output");std::ofstream out(argv[2]);
  for(std::string z:{"0","1","-1","7/4","2","-2","10","-10","1000000"}){auto J=log_moments(fraction(z));out<<"{\"z\":\""<<z<<"\",\"moments\":[";
   for(int n=0;n<=N;++n){if(n)out<<',';out<<"[\""<<endpoint(J[n].leftBound(),false)<<"\",\""<<endpoint(J[n].rightBound(),true)<<"\"]";}out<<"]}\n";}
  if(!out)throw std::runtime_error("moment output");return 0;
 }
 if(argc!=6&&argc!=8)throw std::runtime_error("mode inputs rule output position [polynomials Q-input]");
 int wanted=std::stoi(argv[5]),index,cell;std::string bs,ss,es,ps;std::ifstream input(argv[2]);MI beta,Sref,eta,Pcut;bool found=false;
 while(input>>index>>cell>>bs>>ss>>es>>ps){if(index!=wanted)continue;beta=binary(bs);Sref=binary(ss);eta=binary(es);Pcut=binary(ps);found=true;break;}
 if(!found||beta.leftBound()<=MF(0)||Sref.leftBound()<=MF(0)||Pcut.leftBound()<=MF(0))throw std::runtime_error("state input");
 std::ofstream out(argv[4]);if(!out)throw std::runtime_error("output open");auto started=std::chrono::steady_clock::now();
 if(mode=="build"){
  std::ifstream rule(argv[3]);std::string lo,hi,s;std::array<MI,N>nodes;std::array<std::array<MI,N>,N>connection,chebyshev;
  for(int j=0;j<N;++j){if(!(rule>>lo>>hi))throw std::runtime_error("nodes");nodes[j]=encoded(lo,hi);chebyshev[j][0]=1;chebyshev[j][1]=nodes[j];for(int k=2;k<N;++k)chebyshev[j][k]=2*nodes[j]*chebyshev[j][k-1]-chebyshev[j][k-2];}
  for(int k=0;k<N;++k)for(int n=0;n<N;++n){if(!(rule>>s))throw std::runtime_error("connection matrix");connection[k][n]=fraction(s);}
  MI factor=MI(2)/power(MI(6),N),allowance=encoded("1e-18","1e-18")/Pcut,covered(0);Six maxerror{};int count=0,omitted=0;std::vector<Panel>stack={{MI(0),Pcut,0}};
  while(!stack.empty()){
   Panel panel=stack.back();stack.pop_back();MI a=panel.a,b=panel.b,c=(a+b)/2,h=(b-a)/2,R=4*h;exact(c);exact(h);if(panel.depth>48)throw std::runtime_error("interpolation depth");
   auto split=[&](){stack.push_back({c,b,panel.depth+1});stack.push_back({a,c,panel.depth+1});};MI lower=c-R;if(lower.leftBound()<MF(0))lower=0;
   MI g2=1+lower*lower-R*R;if(g2.leftBound()<=MF(0)){split();continue;}MI g=sqrt(g2),pmax=c+R,phase=pmax*R/(beta*g);if(phase.rightBound()>MF(1)){split();continue;}
   MI tmin=(g2-1)/(beta*(g+1)),tmax=pmax*pmax/(beta*(g+1)),exponent=eta-tmin,qmax(1),L;
   if(exponent.rightBound()>=MF(0))L=1+MI(exponent.rightBound());else{MF e=exponent.rightBound();if(e<MF(-128))e=MF(-128);qmax=exp(MI(e));L=qmax;}
   Six M={qmax/g+2*beta*L,qmax*(1/g+2*beta),qmax*(1/g+2*beta),tmax*qmax/g+2*beta*(L+tmax*qmax),tmax*qmax/g+2*beta*qmax*(1+tmax),
          (tmax*tmax+tmax)*qmax/g+2*beta*(L+(tmax+tmax*tmax)*qmax)};
   Six error{};bool skip=true,accept=true;for(int k=0;k<6;++k){M[k]/=Sref;error[k]=factor*M[k];if(M[k].rightBound()>allowance.leftBound())skip=false;if(error[k].rightBound()>allowance.leftBound())accept=false;}
   if(!skip&&!accept){split();continue;}Coeff coefficients{};
   if(skip){++omitted;error=M;}else{
    Coeff values{},cheb{};
    for(int j=0;j<N;++j){MI p=c+h*nodes[j],gamma=sqrt(1+p*p),t=p*p/(beta*(gamma+1)),q=1/(1+exp(t-eta)),v=1-q,logf=log(1+exp(eta-t));
     values[j]={q/gamma+2*beta*logf,q*v/gamma+2*beta*q,q*v*(1-2*q)/gamma+2*beta*q*v,
      t*q*v/gamma+2*beta*(logf+t*q),t*q*v*(1-2*q)/gamma+2*beta*(q+t*q*v),
      (t*t*(1-2*q)-t)*q*v/gamma+2*beta*(logf+t*q+t*t*q*v)};
     for(auto&value:values[j])value/=Sref;
    }
    for(int n=0;n<N;++n)for(int j=0;j<N;++j)for(int k=0;k<6;++k)cheb[n][k]+=values[j][k]*chebyshev[j][n]*MI(n==0?1:2)/N;
    for(int n=0;n<N;++n)for(int j=n;j<N;++j){if(connection[j][n].leftBound()==MF(0)&&connection[j][n].rightBound()==MF(0))continue;for(int k=0;k<6;++k)coefficients[n][k]+=cheb[j][k]*connection[j][n];}
   }
   out<<point(a)<<' '<<point(b)<<' '<<(skip?1:0);for(int k=0;k<6;++k){out<<' '<<endpoint(error[k].rightBound(),true);if(error[k].rightBound()>maxerror[k].rightBound())maxerror[k]=MI(error[k].rightBound());}
   if(!skip)for(int n=0;n<N;++n)for(int k=0;k<6;++k)out<<' '<<endpoint(coefficients[n][k].leftBound(),false)<<' '<<endpoint(coefficients[n][k].rightBound(),true);
   out<<'\n';covered+=b-a;++count;
  }
  exact(covered);if(covered.leftBound()!=Pcut.leftBound())throw std::runtime_error("H coverage");if(!out)throw std::runtime_error("H output");
  std::cout<<"{\"position\":"<<wanted<<",\"cell\":"<<cell<<",\"panels\":"<<count<<",\"bounded_omissions\":"<<omitted<<",\"uniform_H_errors\":";json_six(std::cout,maxerror);std::cout<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 }else if(mode=="response"){
  if(argc!=8)throw std::runtime_error("response inputs");std::ifstream polynomials(argv[6]),qs(argv[7]);std::vector<Stored>panels;std::string lo,hi,s;int skip;Six maximum{};
  while(polynomials>>lo>>hi>>skip){Stored panel;panel.a=encoded(lo,lo);panel.b=encoded(hi,hi);exact(panel.a);exact(panel.b);panel.skip=skip!=0;
   for(int k=0;k<6;++k){if(!(polynomials>>s))throw std::runtime_error("H error input");panel.error[k]=encoded(s,s);if(panel.error[k].rightBound()>maximum[k].rightBound())maximum[k]=MI(panel.error[k].rightBound());}
   if(!panel.skip)for(int n=0;n<N;++n)for(int k=0;k<6;++k){if(!(polynomials>>lo>>hi))throw std::runtime_error("H coefficient input");panel.c[n][k]=encoded(lo,hi);}panels.push_back(panel);
  }
  int j;std::string qtext;while(qs>>j>>qtext){MI Q=binary(qtext);Six sum{};
   for(auto&panel:panels){if(panel.skip)continue;MI h=(panel.b-panel.a)/2,c=(panel.a+panel.b)/2;exact(h);exact(c);
    if(Q.leftBound()==MF(0)){for(int k=0;k<6;++k)sum[k]+=2*h*panel.c[0][k];continue;}
    auto plus=log_moments(-(c+Q)/h),minus=log_moments((Q-c)/h);std::array<Six,N+1>b{};
    for(int n=0;n<N;++n)for(int k=0;k<6;++k){b[n][k]+=c*panel.c[n][k];b[n+1][k]+=h*MI(n+1)/(2*n+1)*panel.c[n][k];if(n>0)b[n-1][k]+=h*MI(n)/(2*n+1)*panel.c[n][k];}
    for(int n=0;n<=N;++n)for(int k=0;k<6;++k)sum[k]+=h/(2*Q)*b[n][k]*(plus[n]-minus[n]);
   }
   for(int k=0;k<6;++k)sum[k]+=symmetric(Pcut*maximum[k]);out<<"{\"position\":"<<wanted<<",\"cell\":"<<cell<<",\"z_index\":"<<j<<",\"finite_response\":";json_six(out,sum);out<<"}\n";
  }
  if(!out)throw std::runtime_error("response output");std::cout<<"{\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 }else throw std::runtime_error("unknown mode");
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
