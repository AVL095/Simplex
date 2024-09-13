#include <string.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector> 
#include <functional>

using namespace std;

double normal_dist_pdf(double,double,double);
double normal_dist_cdf(double,double,double);
double normal_dist_inv_cdf(double,double,double);
int neldermead(vector<double>&x0, double eps,double(*func)(vector<double>));
vector<vector<double>>CovMatrixMleN(int n,vector<double>x,vector<int>r,double a, double s);
vector<vector<double>>CovMatrixMleW(int n,vector<double>x,vector<int>r,double c, double b);
double NormalMinFunction(vector<double>);
double WeibullMinFunction(vector<double>);
vector<vector<double>>InverseMatrix(vector<vector<double>>,int);
void MLE_Normal(string);
void MLE_Weibull(string);

struct ne_simp {
 int n;
 vector <double>x;
 vector<int>r;
};
ne_simp nesm;

//######################################################

double  NormalMinFunction(vector<double>xsimpl) {
    double s1,s2,s3,s4,z,psi,p,d,c1,c2;
    int i,kx;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
    if (xsimpl[0]<=0) return 10000;
    if (xsimpl[1]<=0) return 10000;

    for (i=0;i<nesm.n;i++) {
            z = (nesm.x[i]-xsimpl[0])/xsimpl[1];
            d = normal_dist_pdf(z,0,1);
            p = normal_dist_cdf(z,0,1);
            psi = d / (1. - p);
            s1 +=(1.-nesm.r[i])*(nesm.x[i]-xsimpl[0]);
            s2 += (1.-nesm.r[i])*pow(nesm.x[i]-xsimpl[0],2);
            s3 += nesm.r[i]*psi;
            s4 += nesm.r[i]*psi*z;
            kx+=1-nesm.r[i];
    }
    c1=s1+xsimpl[1]*s3;
    c2=s2+pow(xsimpl[1],2)*(s4-kx);
    z=c1*c1+c2*c2;
    return z;
}


//#################MLE Weibull Minimized Function#######################

double WeibullMinFunction(vector<double>xsimpl) {
 double s1,s2,s3,z,b,c;
 int i,k;
 if (xsimpl[0]<=0) return(10000000.);
   s1=0;s2=0;s3=0;k=0;
   b=xsimpl[0];
 for(i=0;i<nesm.n;i++) {
   k+=(1-nesm.r[i]);
   s1+=pow(nesm.x[i],b);
 }
   c=s1/k;

for(i=0;i<nesm.n;i++) {
  z=(pow(nesm.x[i],b))/c;
  s3+=z*log(z);
  s2+=(1-nesm.r[i])*log(z);
}
 c=s3-s2-k;
 return c*c;
}


//##########################################################

vector<vector<double>>CovMatrixMleW(int n,vector<double> x,vector<int>r,double c, double b) {

    int i, k;
    double s1, s2, z, cpw, ckow;
    cpw = log(c); ckow = 1 / b; s1 = 0; s2 = 0; k = 0;
    vector<vector<double>>v(n,vector<double>(n));
    for (i = 0; i < n; i++) {
        z = (log(x[i]) - cpw) / ckow;
        s1 += (1 -r[i]) * z;
        s2 += z * z * exp(z);
        k += (1 - r[i]);
    }
    v[0][0]=double(k)/double(n);v[0][1]=(k + s1)/n;
    v[1][0]=(k+s1)/n;v[1][1]=(k+s2)/n;
    v=InverseMatrix(v,2);
    return v;
}

//############################################

vector<vector<double>>CovMatrixMleN(int n,vector<double> x,vector<int>r,double a,double s) {

    double z, p, d, s1, s2, s3, psi;
    int j, k;
    s1 = 0; s2 = 0; s3 = 0; k = 0;
    vector<vector<double>>v(n,vector<double>(n));
    for (j = 0; j < n; j++) {
        z = (x[j] - a) / s;
        p = normal_dist_cdf(z,0,1);
        d =normal_dist_pdf(-0.5*z*z,0,1); 
        psi = d / (1 - p);
        s1 += r[j] * psi * (psi - z);
        s2 += r[j] * psi * z * (z * (psi - z) - 1);
        s3 += r[j] * psi * (z * (psi - z) - 1);
        k += (1 - r[j]);
    }

    v[0][0]=(k+s1)/n;v[0][1]=s3/n;
    v[1][0]=s3/n;v[1][1]=(2*k+s2)/n;
    v=InverseMatrix(v,2);
    return v;
   
}

//######################################################################

 vector<vector<double>>InverseMatrix(vector<vector<double>>a,int n) {
        double temp;
        int i,j,k;

        vector<vector<double>>e(n,vector<double>(n));

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                e[i][j] = 0;
                 if (i == j)
                    e[i][j] = 1;
            }

        for (k = 0; k < n; k++) {
            temp = a[k][k];
            for (j = 0; j < n; j++) {
                a[k][j] /= temp;
                e[k][j] /= temp;
            }
 
            for (i = k + 1; i < n; i++) {
                temp = a[i][k];
                for (j = 0; j < n; j++) {
                    a[i][j] -= a[k][j] * temp;
                    e[i][j] -= e[k][j] * temp;
                }
            }
        }
       for (k = n - 1; k > 0; k--) {
            for (i = k - 1; i >= 0; i--)  {
                temp = a[i][k];
                 for (j = 0; j < n; j++)  {
                    a[i][j] -= a[k][j] * temp;
                    e[i][j] -= e[k][j] * temp;
                }
            }
        }
    return e; 
    }

//##########################################################

double normal_dist_pdf(double x,double mu,double sigma) {
    double diff,pi;
    pi=3.14159265358979;
    diff=(x-mu)/sigma;
    return(exp(-0.5*diff * diff) /(sigma*sqrt(2.*pi)));
}
//############################################

double normal_dist_cdf(double x,double mu,double sigma) {
    return(0.5*(1.0+erf((x-mu)/(sigma *sqrt(2.0)))));
}

//###################################################


double normal_dist_inv_cdf(double p, double mu, double sigma) {

    double q,r,num,x,den;
  
    q = p - 0.5;
    if(fabs(q) <= 0.425) {
        r = 0.180625 - q * q;
        num = (((((((2.5090809287301226727e+3 * r +3.3430575583588128105e+4) * r + 6.7265770927008700853e+4) * r + 4.5921953931549871457e+4) * r +
                     1.3731693765509461125e+4) * r +1.9715909503065514427e+3) * r +1.3314166789178437745e+2) * r +3.3871328727963666080e+0) * q;
        den = (((((((5.2264952788528545610e+3 * r +2.8729085735721942674e+4) * r + 3.9307895800092710610e+4) * r + 2.1213794301586595867e+4) * r +
                     5.3941960214247511077e+3) * r + 6.8718700749205790830e+2) * r + 4.2313330701600911252e+1) * r +1.0);
        x = num / den;
        return(mu + (x * sigma));
     }

    if(q <= 0.0)  {
      r=p;
   }
    else {
      r=1.0-p;
    }
    r =sqrt(-log(r));
    if(r <= 5.0) {
        r = r - 1.6;
        num = (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r +2.41780725177450611770e-1) * r +1.27045825245236838258e+0) * r +
                     3.64784832476320460504e+0) * r + 5.76949722146069140550e+0) * r +4.63033784615654529590e+0) * r +1.42343711074968357734e+0);
        den = (((((((1.05075007164441684324e-9 * r +5.47593808499534494600e-4) * r +1.51986665636164571966e-2) * r + 1.48103976427480074590e-1) * r +
                     6.89767334985100004550e-1) * r +1.67638483018380384940e+0) * r +2.05319162663775882187e+0) * r +
                     1.0);
    } else {
        r = r - 5.0;
        num = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +1.24266094738807843860e-3) * r + 2.65321895265761230930e-2) * r +
                     2.96560571828504891230e-1) * r +1.78482653991729133580e+0) * r + 5.46378491116411436990e+0) * r + 6.65790464350110377720e+0);
        den = (((((((2.04426310338993978564e-15 * r +1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r +7.86869131145613259100e-4) * r +
                     1.48753612908506148525e-2) * r +1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) * r +1.0);
    }
    x = num / den;
    if(q < 0.0)  x = -x;
    return(mu + (x * sigma));

}

//###################################################################

int neldermead(vector<double>&x0, double eps,double(*func)(vector<double>)) {

    double rho,chi, psi, sigma, nonzdelt, zdelt;
    double maxfun,fval,fxr,fxe,fxc,fxcc;
    int i,j,k,N,maxiter,iterations,doshrink;

    rho = 1;chi = 2;psi = 0.5;sigma = 0.5;nonzdelt = 0.05;zdelt = 0.00025;
    N = x0.size();maxiter = N * 200; maxfun = 1.e20;

    vector<vector<double>>sim(N+1,vector<double>(N,0.0));
    vector<double> fsim(N + 1);
    sim[0]= x0;
    vector<double>y = x0;

    for (k = 0; k < N; k++) {
        if (y[k] != 0) {
            y[k] = (1 + nonzdelt) * y[k];
        } else {
            y[k] = zdelt;
        }
        sim[k + 1] = y;
    }
  
    for (i = 0; i < N + 1; i++)   fsim[i] = func(sim[i]);

    iterations = 0;
    fval = maxfun;

    while (true) {
        if (fval <= eps || iterations >= maxiter) break;
        vector<double> xbar(N,0.0);
        for (i = 0; i < N; i++) {
            for ( j = 0; j < N; j++) xbar[j] += sim[i][j];
        }
        for (i = 0; i < N; i++) xbar[i] /= N;
        
        vector<double> xr(N);
        for (i = 0; i < N; i++) xr[i] = (1 + rho) * xbar[i] - rho * sim[N][i];
        fxr = func(xr);

        doshrink = 0;
        if (fxr < fsim[0]) {
            vector<double> xe(N);
            for (i = 0; i < N; i++)  xe[i] = (1 + rho * chi) * xbar[i] - rho * chi * sim[N][i];
            fxe = func(xe);
            if (fxe < fxr) {
                sim[N] = xe;fsim[N] = fxe;
            } else {
                sim[N] = xr;fsim[N] = fxr;
            }
        } else {
            if (fxr < fsim[N - 1]) {
                sim[N] = xr;fsim[N] = fxr;
            } else {
                if (fxr < fsim[N]) {
                    vector<double> xc(N);
                    for (i = 0; i < N; i++) xc[i] = (1 + psi * rho) * xbar[i] - psi * rho * sim[N][i];
                    fxc = func(xc);
                    if (fxc <= fxr) {
                        sim[N] = xc; fsim[N] = fxc;
                    } else {
                        doshrink = 1;
                    }
                } else {
                    vector<double> xcc(N);
                    for (i = 0; i < N; i++)  xcc[i] = (1 - psi) * xbar[i] + psi * sim[N][i];
                    fxcc = func(xcc);
                    if (fxcc < fsim[N]) {
                        sim[N] = xcc;fsim[N] = fxcc;
                    } else {
                        doshrink = 1;
                    }
                }
                if (doshrink) {
                    for (j = 1; j < N + 1; j++) {
                        for (i = 0; i < N; i++)     sim[j][i] = sim[0][i] + sigma * (sim[j][i] - sim[0][i]);
                        fsim[j] = func(sim[j]);
                    }
                }
            }
        }

        iterations++;
        fval = func(sim[0]);
        vector<int> ind(N + 1);
        for (i = 0; i < N + 1; i++) ind[i] = i;
        sort(ind.begin(), ind.end(), [&](int i, int j) { return fsim[i] < fsim[j]; });
        vector<vector<double>> sim_sorted(N + 1, vector<double>(N));
        vector<double> fsim_sorted(N + 1);
        for (i = 0; i < N + 1; i++) {
            sim_sorted[i] = sim[ind[i]];
            fsim_sorted[i] = fsim[ind[i]];
        }
        sim = sim_sorted; fsim = fsim_sorted;

    }  //end while


    x0=sim[0];
    return iterations;

}

//###############################################################

void MLE_Normal(string ff) {
    
    int i,j,k,icount;
    string s1;
    double cp,cko,q,eps,z;

    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");
    inp >> s1; //text
    inp >>nesm.n; //sample size
    inp >> s1;
    for (i = 0; i <nesm.n; i++) {inp>>z;nesm.x.push_back(z);}
    inp >> s1; //text
    for (i=0;i<nesm.n; i++) {inp>>j;nesm.r.push_back(j);}
    inp.close();
//###########################################################
 
     vector<double>xsimpl; 
     vector<vector<double>>v(2,vector<double>(2));

    cp=0;cko=0;k=0;
    for(i=0;i<nesm.n;i++) {
        k+=(1-nesm.r[i]); // количество наблюдений
        cp+=(1-nesm.r[i])*nesm.x[i];
        cko+=(1-nesm.r[i])*nesm.x[i]*nesm.x[i];
    }
    cp/=k; //выборочное среднее по наблюдениям
    cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям
    xsimpl.push_back(cp);xsimpl.push_back(cko);
    v[0][0]=1.;v[1][1]=0.5;v[0][1]=0.;v[1][0]=0.;q=0;eps=1.e-15;icount=0;
    if(k!=nesm.n) {
     icount=neldermead(xsimpl,eps,NormalMinFunction);
     q=NormalMinFunction(xsimpl);
     v=CovMatrixMleN(nesm.n,nesm.x,nesm.r,xsimpl[0],xsimpl[1]);
    }

 out << "Method:" << ff << "\n";
 out << "n=" <<nesm.n << "\n";
 out << "X" << "\n";
 for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
 out << "\n";
 out << "R" << "\n";
 for (i = 0; i<nesm.n; i++)  out <<nesm.r[i] << " , ";
 out << "\n";
 out << "cp*=" <<setprecision(12) << fixed <<cp << "\n";
 out << "cko*="<<setprecision(12) << fixed <<cko << "\n";
 out << "Q="<<setprecision(12) << q << "\n";
 out << "icount="<<icount<<endl;
 out << "cp="<<setprecision(12) << fixed << xsimpl[0]<<"\n";
 out << "cko="<<setprecision(12) << fixed <<xsimpl[1]<<"\n";
 out << "v11=" << v[0][0] << "\n";
 out << "v12=" << v[0][1] << "\n";
 out << "v21=" << v[1][0] << "\n";
 out << "v22=" << v[1][1] << "\n";
 out.close();

 v.clear();nesm.r.clear();nesm.x.clear();xsimpl.clear();

}

//#########################################################################

void MLE_Weibull(string ff) {
 int i,j,k,icount;
 string s1;
 double cp, cko,eps,q,s,b,c,aw,sw,z;

 vector<double>logx;
 vector<double>xsimpl;
//##################################################################		   
 ifstream inp("Inp/" + ff + ".inp");
 inp >> s1;
 inp >>nesm.n;
 inp >> s1;
 for(i=0;i<nesm.n;i++) {inp>>z;nesm.x.push_back(z);}
 inp >> s1;
 for(i= 0;i<nesm.n;i++) {inp>>j;nesm.r.push_back(j);}
 inp.close();
//###############################################################
 cp = 0;cko=0;k=0;
 for(i=0;i<nesm.n;i++) logx.push_back(log(nesm.x[i]));
 for(i=0;i<nesm.n;i++) {
   k +=(1-nesm.r[i]); // количество наблюдений
   cp+=(1-nesm.r[i])*logx[i];
   cko+=(1-nesm.r[i])*logx[i]*logx[i];
   }
   cp/=k; //выборочное среднее по наблюдениям
   cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям

   vector<vector<double>>v(2,vector<double>(2,0.0));

 q=0;xsimpl.push_back(1./cko);eps=1.e-15;icount=0;
 icount=neldermead(xsimpl,eps,WeibullMinFunction);
 q=WeibullMinFunction(xsimpl);

 b=xsimpl[0];
 s = 0.;
 for(i=0;i<nesm.n;i++) s+=pow(nesm.x[i],b);
 c=pow(s/k,1/b);
 v=CovMatrixMleW(nesm.n,nesm.x,nesm.r,c,b);
 sw=1./b;aw=log(c);

 ofstream out("Out/" + ff + ".out");
  out << "Method:" << ff << "\n";
  out << "n=" <<nesm.n << "\n";
  out << "X" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
  out << "\n";
  out << "log(X)" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<logx[i] << " , ";
  out << "\n";
  out << "R" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<nesm.r[i] << " , ";
  out << "\n";
  out << "cp*="<<setprecision(12)<< fixed<<cp << "\n";
  out << "cko*="<<setprecision(12)<< fixed<<cko << "\n";
  out << "Q=" << q << "\n";
  out <<"icount="<<icount<<endl;
  out << "b="<<setprecision(12)<< fixed <<b<<"\n";
  out << "c="<<setprecision(12) << fixed <<c<< "\n";
  out << "aw="<<setprecision(12)<< fixed <<aw<<"\n";
  out << "sw="<<setprecision(12) << fixed <<sw<< "\n";
  out << "\n";
  out << "v11=" << v[0][0] << "\n";
  out << "v12=" << v[0][1] << "\n";
  out << "v21=" << v[1][0] << "\n";
  out << "v22=" << v[1][1] << "\n";
  out.close();
 
  logx.clear();nesm.r.clear();xsimpl.clear();nesm.x.clear();v.clear();

}

//########################################################


int main() {

 string ff;
 ifstream inp("main.inp");
 inp>>ff;
 inp.close();
 if(ff=="MLE_Normal") MLE_Normal(ff);
 if(ff=="MLE_Weibull") MLE_Weibull(ff);

return 0;
}


