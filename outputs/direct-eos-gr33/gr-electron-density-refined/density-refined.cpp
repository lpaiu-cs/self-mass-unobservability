// Reuse the CAPD/MPFR Gauss rule and directed output from gr_plasma_interval.cpp.
#include "capd/intervals/MpInterval.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using MI=capd::intervals::MpInterval;
using MF=capd::multiPrec::MpReal;

std::string endpoint(const MF& x,bool up){
  std::ostringstream out;x.put(out,up?MF::RoundUp:MF::RoundDown,128,10,60);return out.str();
}

int main(int argc,char**argv){
 try{
  if(argc!=6)throw std::runtime_error("args: input output panels first count");
  int panels=std::stoi(argv[3]),first=std::stoi(argv[4]),count=std::stoi(argv[5]);
  if(panels!=8192||first<0||count<1)throw std::runtime_error("outside frozen rule");
  MF::setDefaultPrecision(128);
  MI outer=sqrt((3+2*sqrt(MI(6)/5))/7),inner=sqrt((3-2*sqrt(MI(6)/5))/7);
  MI nodes[4]={-outer,-inner,inner,outer};
  MI weights[4]={(18-sqrt(MI(30)))/36,(18+sqrt(MI(30)))/36,(18+sqrt(MI(30)))/36,(18-sqrt(MI(30)))/36};
  for(int k=0;k<8;++k){
   MI v(0),exact=k%2?MI(0):MI(2)/(k+1);
   for(int i=0;i<4;++i){MI power(1);for(int j=0;j<k;++j)power*=nodes[i];v+=weights[i]*power;}
   if(v.leftBound()>exact.leftBound()||v.rightBound()<exact.rightBound())throw std::runtime_error("Gauss control");
  }
  std::vector<MI> tnodes,wq;MI h=MI(10)/panels;
  for(int j=0;j<panels;++j)for(int k=0;k<4;++k){
   MI s=h*(MI(j)+MI(1)/2+nodes[k]/2);tnodes.push_back(s*s);wq.push_back(h*weights[k]);
  }
  std::ifstream in(argv[1]);std::ofstream out(argv[2]);
  if(!in||!out)throw std::runtime_error("file open failure");
  int id,cell,done=0;std::string es,bs;MI radius=MI(1)/MI(70368744177664.0);
  while(in>>id>>cell>>es>>bs){
   if(id<first||id>=first+count)continue;
   double ev=std::stod(es),bv=std::stod(bs);
   MI center{MF(ev)},beta{MF(bv)},eta=center+MI(MF(-1),MF(1))*radius;
   if(eta.leftBound()<MF(-17)||eta.rightBound()>MF(24)||!(bv>0&&bv<=.006))throw std::runtime_error("outside certificate box");
   MI sum[6]={MI(0),MI(0),MI(0),MI(0),MI(0),MI(0)},point(0);
   auto start=std::chrono::steady_clock::now();
   for(size_t i=0;i<tnodes.size();++i){
    MI t=tnodes[i],bt=beta*t,w=1+bt/2,s=sqrt(w),q=1/(1+exp(t-eta)),qc=1/(1+exp(t-center));
    MI base=wq[i]*t,H=s*(1+bt),Hb=t*s+t*(1+bt)/(4*s),Hbb=t*t/(2*s)-t*t*(1+bt)/(16*w*s);
    MI qe=q*(1-q),qee=qe*(1-2*q);
    point+=base*H*qc;
    sum[0]+=base*H*q;sum[1]+=base*H*qe;sum[2]+=base*Hb*q;
    sum[3]+=base*H*qee;sum[4]+=base*Hb*qe;sum[5]+=base*Hbb*q;
   }
   out<<"{\"position\":"<<id<<",\"cell\":"<<cell<<",\"bits\":128,\"panels\":"<<panels
      <<",\"center_D\":[\""<<endpoint(point.leftBound(),false)<<"\",\""<<endpoint(point.rightBound(),true)<<"\"],\"box_D_partials\":[";
   for(int j=0;j<6;++j){if(j)out<<',';out<<"[\""<<endpoint(sum[j].leftBound(),false)<<"\",\""<<endpoint(sum[j].rightBound(),true)<<"\"]";}
   out<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
   if(!out)throw std::runtime_error("output failure");++done;
  }
  if(done!=count)throw std::runtime_error("incomplete coverage");
  std::cout<<"PASS certified electron density rows "<<done<<"\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
