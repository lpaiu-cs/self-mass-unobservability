// Evaluate the declared electron moment quadrature with existing MPFR intervals.
#include "capd/intervals/MpInterval.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using MI=capd::intervals::MpInterval;
using MF=capd::multiPrec::MpReal;

std::string endpoint(const MF& value, bool up){
  std::ostringstream out;
  value.put(out,up?MF::RoundUp:MF::RoundDown,128,10,60);
  return out.str();
}

int main(int argc,char**argv){
 try{
  if(argc!=6)throw std::runtime_error("args: input output panels first count");
  const int panels=std::stoi(argv[3]),first=std::stoi(argv[4]),count=std::stoi(argv[5]);
  if(panels<256||panels>32768||panels%256||first<0||count<1)throw std::runtime_error("invalid bounds");
  MF::setDefaultPrecision(128);
  MI outer=sqrt((3+2*sqrt(MI(6)/5))/7),inner=sqrt((3-2*sqrt(MI(6)/5))/7);
  MI nodes[4]={-outer,-inner,inner,outer};
  MI weights[4]={(18-sqrt(MI(30)))/36,(18+sqrt(MI(30)))/36,
                 (18+sqrt(MI(30)))/36,(18-sqrt(MI(30)))/36};
  for(int k=0;k<8;++k){
    MI value(0),target=k%2?MI(0):MI(2)/(k+1);
    for(int i=0;i<4;++i){MI power(1);for(int j=0;j<k;++j)power*=nodes[i];value+=weights[i]*power;}
    if(value.leftBound()>target.leftBound()||value.rightBound()<target.rightBound())throw std::runtime_error("Gauss moment control");
  }
  MI unit=exp(MI(0)),square=sqrt(MI(4));
  if(unit.leftBound()!=MF(1)||unit.rightBound()!=MF(1)||square.leftBound()!=MF(2)||square.rightBound()!=MF(2))throw std::runtime_error("elementary control");
  std::vector<MI> tnodes,wq;MI h=MI(10)/panels;
  for(int j=0;j<panels;++j)for(int k=0;k<4;++k){
    MI s=h*(MI(j)+MI(1)/2+nodes[k]/2);tnodes.push_back(s*s);wq.push_back(h/2*weights[k]);
  }
  std::ifstream in(argv[1]);std::ofstream out(argv[2]);
  if(!in||!out)throw std::runtime_error("file open failure");
  int id,done=0;std::string es,bs;
  while(in>>id>>es>>bs){
    if(id<first||id>=first+count)continue;
    double ev=std::stod(es),bv=std::stod(bs);
    if(!(ev>=-17&&ev<=24&&bv>0&&bv<=.006))throw std::runtime_error("outside certified box");
    MI eta{MF(ev)},beta{MF(bv)};MI norm=ev<0?exp(-eta):MI(1);
    std::vector<MI> sum(26,MI(0));auto start=std::chrono::steady_clock::now();
    for(size_t i=0;i<tnodes.size();++i){
      MI t=tnodes[i],bt=beta*t,x2=bt*(2+bt),gamma=1+bt,w=x2/(gamma*gamma);
      MI value=wq[i]*2*t*sqrt(1+bt/2)*norm/(1+exp(t-eta));
      for(int j=0;j<26;++j){sum[j]+=value;value*=w;}
    }
    out<<"{\"cell\":"<<id<<",\"bits\":128,\"panels\":"<<panels<<",\"moments\":[";
    for(int j=0;j<26;++j){if(j)out<<',';out<<"[\""<<endpoint(sum[j].leftBound(),false)<<"\",\""<<endpoint(sum[j].rightBound(),true)<<"\"]";}
    out<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
    if(!out)throw std::runtime_error("output failure");++done;
  }
  if(done!=count)throw std::runtime_error("incomplete selected input coverage");
  std::cout<<"PASS MPFR interval moment rows "<<done<<"\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
