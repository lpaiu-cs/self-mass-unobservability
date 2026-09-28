// Request 14. GR specialization of Nutimo rhs_GR_nbody, GPL-3.0-or-later.
// Original equation implementation: Guillaume VOISIN (2024).
// The input reader rejects non-GR/static-SEP and dynamic-SEP coefficients.
#include "capd/capdlib.h"
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using capd::autodiff::Node;
using capd::interval;
using Vec = std::array<Node,3>;
constexpr int N=4, D=24, NP=6;
int direction=1;

Node dot(const Vec& a,const Vec& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void field(Node, Node in[], int, Node out[], int, Node p[], int) {
    // p = c, c^2, and four geometric masses. Preserve the native c^2 constant.
    Vec v[N], n[N][N]; Node vv[N], r[N][N];
    for(int i=0;i<N;++i) {
        for(int a=0;a<3;++a) v[i][a]=in[3*N+3*i+a]/p[0];
        vv[i]=dot(v[i],v[i]);
        for(int j=0;j<i;++j) {
            Vec d; for(int a=0;a<3;++a) d[a]=in[3*j+a]-in[3*i+a];
            r[i][j]=sqrt(dot(d,d)); r[j][i]=r[i][j];
            for(int a=0;a<3;++a) {n[i][j][a]=d[a]/r[i][j]; n[j][i][a]=-n[i][j][a];}
        }
    }
    for(int i=0;i<N;++i) {
        Vec acc{Node(0),Node(0),Node(0)};
        for(int j=0;j<N;++j) if(i!=j) {
            Node vn=dot(v[j],n[i][j]);
            Node s=1-4*dot(v[i],v[j])+vv[i]+2*vv[j]-1.5*vn*vn;
            Vec w; for(int a=0;a<3;++a) w[a]=4*v[i][a]-3*v[j][a];
            Node q=dot(n[i][j],w), pref=p[2+j]/(r[i][j]*r[i][j]);
            for(int a=0;a<3;++a) acc[a]=acc[a]+pref*(n[i][j][a]*s+(v[j][a]-v[i][a])*q);
            for(int k=0;k<N;++k) {
                if(k!=j) {
                    s=0.5*dot(n[i][j],n[j][k])/r[j][k]-1/r[i][j];
                    q=3.5/r[j][k];
                    pref=p[2+j]*p[2+k]/(r[i][j]*r[j][k]);
                    for(int a=0;a<3;++a) acc[a]=acc[a]+pref*(s*n[i][j][a]+q*n[j][k][a]);
                }
                if(k!=i) {
                    pref=4*p[2+j]*p[2+k]/(r[i][j]*r[i][j]*r[i][k]);
                    for(int a=0;a<3;++a) acc[a]=acc[a]-pref*n[i][j][a];
                }
            }
        }
        for(int a=0;a<3;++a) {out[3*i+a]=direction*in[3*N+3*i+a]; out[3*N+3*i+a]=direction*acc[a]*p[1];}
    }
}

long double readhex(std::istream& f) {
    std::string s; if(!(f>>s)) throw std::runtime_error("truncated input");
    char* end; long double x=strtold(s.c_str(),&end);
    if(*end || !std::isfinite(x)) throw std::runtime_error("invalid input");
    return x;
}
interval enclose(long double x) {
    double d=static_cast<double>(x),lo=d,hi=d;
    if(static_cast<long double>(d)>x) lo=std::nextafter(d,-INFINITY);
    if(static_cast<long double>(d)<x) hi=std::nextafter(d,INFINITY);
    if(!(static_cast<long double>(lo)<=x && x<=static_cast<long double>(hi))) throw std::runtime_error("conversion not enclosed");
    return interval(lo,hi);
}
double width(interval x) {return (interval(x.rightBound())-interval(x.leftBound())).rightBound();}
void pair(std::ostream& f,interval x) {f<<'['<<x.leftBound()<<','<<x.rightBound()<<']';}

void control() {
    capd::IMap f("var:x,y;fun:y,-x;"); capd::IOdeSolver solver(f,20);
    capd::ITimeMap tm(solver); capd::IVector u(2); u[0]=1;u[1]=0;
    capd::C1Rect2Set s(u); tm(interval(1),s);
    capd::IVector y=s; capd::IMatrix j=s;
    // Exact-rational alternating-series bounds are checked by the build producer.
    const interval co(0x1.14a280fb5068bp-1,0x1.14a280fb5068dp-1);
    const interval si(0x1.aed548f090cedp-1,0x1.aed548f090cefp-1);
    auto contains=[](interval a,interval b){return a.leftBound()<=b.leftBound()&&a.rightBound()>=b.rightBound();};
    if(!contains(y[0],co)||!contains(y[1],-si)||!contains(j[0][0],co)||!contains(j[0][1],si)||!contains(j[1][0],-si)||!contains(j[1][1],co)) throw std::runtime_error("oscillator enclosure control failed");
    std::cout<<"{\"oscillator_state_and_jacobian_control\":true}\n";
}

int main(int argc,char** argv) {
    std::cout<<std::setprecision(17);
    try {
        control();
        if(argc==1) return 0;
        if(argc!=5) throw std::runtime_error("args: ivp.hex native-rhs.hex horizon output.jsonl");
        std::ifstream ivp(argv[1]); int bodies;ivp>>bodies;
        if(bodies!=N) throw std::runtime_error("not four bodies");
        auto t0=readhex(ivp),tmin=readhex(ivp),tmax=readhex(ivp),length=readhex(ivp),timescale=readhex(ivp);
        if(t0!=0) throw std::runtime_error("only t0=0 exported IVP supported");
        long double par[NP];par[0]=readhex(ivp);par[1]=readhex(ivp);
        capd::IVector u(D);for(int i=0;i<D;++i) u[i]=enclose(readhex(ivp));
        for(int i=2;i<NP;++i)par[i]=readhex(ivp);
        for(int i=0;i<N;++i)for(int j=0;j<N;++j) {
            auto g=readhex(ivp),gamma=readhex(ivp);
            if((i!=j && g!=1)||gamma!=0) throw std::runtime_error("non-GR coefficients");
        }
        for(int i=0;i<N*N*N;++i)if(readhex(ivp)!=0)throw std::runtime_error("nonzero beta");
        auto A=readhex(ivp);readhex(ivp);readhex(ivp);readhex(ivp);
        if(A!=0)throw std::runtime_error("nonzero dynamic coupling");
        capd::IMap f(field,D,D,NP);capd::LDMap fl(field,D,D,NP);
        for(int i=0;i<NP;++i){f.setParameter(i,enclose(par[i]));fl.setParameter(i,par[i]);}
        std::ifstream samples(argv[2]); int count=0,contained=0;long double maxscaled=0;
        while(samples>>std::ws && samples.peek()!=EOF) {
            readhex(samples);capd::LDVector x(D);capd::IVector xi(D);
            for(int i=0;i<D;++i){x[i]=readhex(samples);xi[i]=enclose(x[i]);}
            auto y=fl(x);auto yi=f(xi);
            for(int i=0;i<D;++i) {
                auto z=readhex(samples);maxscaled=std::max(maxscaled,fabsl(y[i]-z)/(1+fabsl(z)));
                if(yi[i].leftBound()<=z&&z<=yi[i].rightBound())++contained;
            }++count;
        }
        // This is a transcription regression, not a proof by finite samples.
        if(count!=64||maxscaled>1e-16L||contained!=count*D)throw std::runtime_error("native RHS regression failed");
        std::cout<<"{\"native_samples\":"<<count<<",\"contained_components\":"<<contained<<",\"max_scaled_long_double_difference\":"<<maxscaled<<",\"native_tmin\":"<<tmin<<",\"native_tmax\":"<<tmax<<",\"length_m\":"<<length<<",\"timescale_s\":"<<timescale<<"}\n";
        double horizon=std::stod(argv[3]); direction=horizon<0?-1:1;
        capd::IMap flow(field,D,D,NP);for(int i=0;i<NP;++i)flow.setParameter(i,enclose(par[i]));
        capd::IOdeSolver solver(flow,20);solver.setAbsoluteTolerance(1e-18);solver.setRelativeTolerance(1e-18);
        capd::ITimeMap tm(solver);tm.stopAfterStep(true);capd::C1Rect2Set s(u);
        std::ofstream log(argv[4]);log<<std::setprecision(17);
        auto start=std::chrono::steady_clock::now(); int steps=0;
        try {
            do {
                tm(interval(std::abs(horizon)),s);++steps;
                capd::IVector y=s;capd::IMatrix j=s;
                double sw=0,jw=0;for(int i=0;i<D;++i){sw=std::max(sw,width(y[i]));for(int k=0;k<D;++k)jw=std::max(jw,width(j[i][k]));}
                interval psq(0);for(int i=0;i<3;++i){interval w(width(y[i]));psq+=w*w;}
                interval roemer_width=sqrt(psq)/enclose(par[0])*enclose(timescale)*interval(1000000);
                if(steps==1||steps%100==0||tm.completed()||jw>1e6) {
                    log<<"{\"step\":"<<steps<<",\"t\":";pair(log,tm.getCurrentTime()*interval(direction));
                    log<<",\"max_state_width\":"<<sw<<",\"max_jacobian_width\":"<<jw<<",\"geometric_delay_width_us_upper\":"<<roemer_width.rightBound()<<",\"state\":[";
                    for(int i=0;i<D;++i){if(i)log<<',';pair(log,y[i]);}log<<"],\"jacobian\":[";
                    for(int i=0;i<D;++i)for(int k=0;k<D;++k){if(i||k)log<<',';pair(log,j[i][k]);}log<<"]}\n"<<std::flush;
                    std::cout<<"{\"step\":"<<steps<<",\"t\":";pair(std::cout,tm.getCurrentTime()*interval(direction));std::cout<<",\"max_state_width\":"<<sw<<",\"max_jacobian_width\":"<<jw<<"}\n"<<std::flush;
                }
                if(jw>1e6)throw std::runtime_error("declared Jacobian-width ceiling 1e6 exceeded; no restart");
                if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>1800)throw std::runtime_error("1800 second run budget; no restart");
            }while(!tm.completed());
            std::cout<<"{\"requested_horizon_completed\":true,\"D2_full_timing_certificate\":false}\n";
        }catch(const std::exception& e){std::cout<<"{\"requested_horizon_completed\":false,\"D2_full_timing_certificate\":false}\n";std::cerr<<e.what()<<'\n';return 2;}
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
