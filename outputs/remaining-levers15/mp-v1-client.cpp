// Request 15: reuse the frozen, audited RHS without changing Request 14.
#define main request14_main
#include "validated_variational.cpp"
#undef main
#include "capd/mpcapdlib.h"
#include "capd/dynsys/MpStepControl.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/TimeMap.hpp"
using MI=capd::MpInterval;
using MF=capd::MpFloat;
using MS=capd::dynsys::OdeSolver<capd::MpIMap,capd::dynsys::MpLastTermsStepControl>;
using MT=capd::poincare::TimeMap<MS>;
MI exact(long double x){return MI(MF(x));} // 64 significand bits, precision >= 128.
double upper(MI x){return toDouble(x.rightBound(),MPFR_RNDU);}
void mpair(std::ostream& f,MI x){f<<'['<<toDouble(x.leftBound(),MPFR_RNDD)<<','<<upper(x)<<']';}
int main(int argc,char**argv){
 try {
  if(argc!=6)throw std::runtime_error("args: ivp.hex horizon bits output.jsonl seconds");
  int bits=std::stoi(argv[3]);if(bits<128)throw std::runtime_error("precision must enclose long double exactly");
  MF::setDefaultPrecision(bits);std::cout<<std::setprecision(17);
  // Nonlinear exact rational positive control, including a nonzero initial box.
  capd::MpIMap cf("var:x;fun:x*x;");MS cs(cf,30);MT ct(cs);
  MI box=MI(1)+MI(-1,1)/1048576;capd::MpIVector cu(1);cu[0]=box;capd::MpC1Rect2Set control(cu);
  ct(MI(1)/4,control);capd::MpIVector cy=control;capd::MpIMatrix cj=control;
  MI ey=box/(1-box/4),den=1-box/4,ej=1/(den*den);
  auto contains=[](MI a,MI b){return a.leftBound()<=b.leftBound()&&a.rightBound()>=b.rightBound();};
  if(!contains(cy[0],ey)||!contains(cj[0][0],ej))throw std::runtime_error("MP nonlinear control failed");
  std::cout<<"{\"mp_nonlinear_state_jacobian_control\":true,\"bits\":"<<bits<<"}\n"<<std::flush;
  std::ifstream ivp(argv[1]);int n;ivp>>n;if(n!=N)throw std::runtime_error("body count");
  auto t0=readhex(ivp);readhex(ivp);readhex(ivp);readhex(ivp);auto timescale=readhex(ivp);
  if(t0!=0)throw std::runtime_error("t0 must be zero");
  long double p[NP];p[0]=readhex(ivp);p[1]=readhex(ivp);
  capd::MpIVector u(D);for(int i=0;i<D;++i)u[i]=exact(readhex(ivp))/(i<12?1:128);
  for(int i=2;i<NP;++i)p[i]=readhex(ivp);
  for(auto a:p)if(a<=0)throw std::runtime_error("nonpositive parameter");
  long double total=p[2];for(int i=0;i<3;++i){total+=p[3+i];mu[i]=static_cast<double>(roundl(p[3+i]/total*0x1p32L)/0x1p32L);}
  for(int i=0;i<N;++i)for(int k=0;k<N;++k){auto g=readhex(ivp),ga=readhex(ivp);if((i!=k&&g!=1)||ga!=0)throw std::runtime_error("non-GR");}
  for(int i=0;i<64;++i)if(readhex(ivp)!=0)throw std::runtime_error("non-GR beta");
  if(readhex(ivp)!=0)throw std::runtime_error("nonzero dynamic coupling");
  readhex(ivp);readhex(ivp);readhex(ivp);
  MI raw[D],out[D];for(int i=0;i<D;++i)raw[i]=u[i];to_jacobi(raw,out);for(int i=0;i<D;++i)u[i]=out[i];
  capd::MpIMatrix S(D,D),T(D,D);
  for(int k=0;k<D;++k){for(int i=0;i<D;++i)raw[i]=i==k?1:0;from_jacobi(raw,out);for(int i=0;i<D;++i)S[i][k]=out[i];to_jacobi(raw,out);for(int i=0;i<D;++i)T[i][k]=out[i];}
  double horizon=std::stod(argv[2]);direction=horizon<0?-1:1;
  capd::MpIMap flow(jacobifield,D,D,NP);for(int i=0;i<NP;++i)flow.setParameter(i,exact(p[i]));
  MS solver(flow,30);solver.setStep(MI(1)/8);MT tm(solver);tm.stopAfterStep(true);
  capd::MpC1Rect2Set s(u);int steps=0;auto start=std::chrono::steady_clock::now();
  std::ofstream log(argv[4]);log<<std::setprecision(17);double budget=std::stod(argv[5]);
  do {
   tm(exact(std::abs(horizon))*128,s);++steps;
   if(steps==1||steps%100==0||tm.completed()){
    capd::MpIVector y=S*static_cast<capd::MpIVector>(s);capd::MpIMatrix j=S*static_cast<capd::MpIMatrix>(s)*T;
    double sw=0,jw=0;MI psq(0);
    for(int i=0;i<D;++i){y[i]*=MI(i<12?1:128);sw=std::max(sw,upper(MI(y[i].rightBound())-MI(y[i].leftBound())));for(int k=0;k<D;++k){j[i][k]*=MI(i<12?1:128)/MI(k<12?1:128);jw=std::max(jw,upper(MI(j[i][k].rightBound())-MI(j[i][k].leftBound())));}}
    for(int i=0;i<3;++i){MI w=MI(y[i].rightBound())-MI(y[i].leftBound());psq+=w*w;}
    log<<"{\"step\":"<<steps<<",\"t\":";mpair(log,tm.getCurrentTime()*direction/128);
    log<<",\"max_state_width\":"<<sw<<",\"max_jacobian_width\":"<<jw<<",\"geometric_delay_width_us_upper\":"<<upper(sqrt(psq)/exact(p[0])*exact(timescale)*1000000)<<",\"state\":[";
    for(int i=0;i<D;++i){if(i)log<<',';mpair(log,y[i]);}log<<"],\"jacobian\":[";
    for(int i=0;i<D;++i)for(int k=0;k<D;++k){if(i||k)log<<',';mpair(log,j[i][k]);}log<<"]}\n"<<std::flush;
    std::cout<<"{\"step\":"<<steps<<",\"t\":";mpair(std::cout,tm.getCurrentTime()*direction/128);std::cout<<",\"max_state_width\":"<<sw<<",\"max_jacobian_width\":"<<jw<<"}\n"<<std::flush;
    if(jw>1e6)throw std::runtime_error("Jacobian width ceiling 1e6; no restart");
   }
   if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>budget)throw std::runtime_error("run budget; no restart");
  }while(!tm.completed());
  std::cout<<"{\"requested_horizon_completed\":true,\"full_timing_certificate\":false}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';std::cout<<"{\"requested_horizon_completed\":false,\"full_timing_certificate\":false}\n";return 2;}
}
