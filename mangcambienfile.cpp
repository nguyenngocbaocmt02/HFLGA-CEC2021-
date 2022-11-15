#include<bits/stdc++.h>
#include<dirent.h>
using namespace std;
#define MUTATION_RATE 0.1
#define CROSSOVER_RATE 0.75
#define POP_SIZE 250
#define GENERATIONS 400
#define nm 2.0
#define nc 2.0
#define GENERATIONSt 1000
#define CROSSOVER_RATE2 0.75
#define MUTATION_RATE2 0.1
#define ENCODE_RATE 0.5
const float Emax =10800,Emin=540,Pm=1,v=5;
double T=20000,U=4,Emc=108000;
string tenfile[109];
double dis[502][502];
struct Sensor {
	double x,y,p,energy,hesosac,w;
};
int n;
vector<Sensor> sensors;
struct Individual{
	vector<double> mhpath;
	vector<int> path;
	vector<double> chargetime;
	double fitness;
	double tsac;
	double tdichuyen;
	double fitnesspath;
void check() {
cout<<"STT"<<"|"<<"nangluong_bandau"<<"|"<<"congsuat"<<"|"<<"thoigianchosac"<<"|"<<"nangluong_batdausac"<<"|"<<"thoigiansac"<<"|"<<"nangluongsaukhisac"<<"|"<<"nangluongT"<<endl;
	  double dichuyen=0,sac=0;
	  vector<bool> charged(n+2,0);
	  cout<<endl<<setprecision(7)<<fixed;
	  for(int i=1;i<path.size()-1;i++) {
	  	dichuyen+=dis[path[i]][path[i-1]]/v;
	  	if(dichuyen+sac>T) break;
		  charged[path[i]]=1;
cout<<path[i]<<"|"<<sensors[path[i]].energy<<"|"<<sensors[path[i]].p<<"|"<<dichuyen+sac<<"|"<<sensors[path[i]].energy-sensors[path[i]].p*(dichuyen+sac)<<"|"<<chargetime[i]<<"|"<<sensors[path[i]].energy-sensors[path[i]].p*(dichuyen+sac)+chargetime[i]*(U-sensors[path[i]].p)<<"|"<<sensors[path[i]].energy-sensors[path[i]].p*(T)+chargetime[i]*(U);
	    sac+=chargetime[i];
	    if(sensors[path[i]].energy-sensors[path[i]].p*(dichuyen+sac)<=Emin||sensors[path[i]].energy-sensors[path[i]].p*(T)+chargetime[i]*(U)<=Emin) cout<<"        CHET";
	cout<<endl;
	}
	for(int i=1;i<=n;i++) {
		if(charged[i]==1) continue;
		cout<<i<<"|"<<sensors[i].energy<<"|"<<sensors[i].p<<"|"<<sensors[i].energy-sensors[i].p*T;
		if(sensors[i].energy-sensors[i].p*T<Emin) cout<<"            CHET";
		cout<<endl;
	}
}
void check_path() {
	for(int i=0;i<path.size();i++) {
		cout<<path[i]<<" ";
	}
	cout<<endl;
}
void check_chargetime() {
	for(int i=0;i<chargetime.size();i++) {
		cout<<chargetime[i]<<" ";
	}
	cout<<endl;
};
};
Individual solution;
struct Population{
	vector<Individual> pop;
	Individual best;
};
Population collect;
int random(int minN, int maxN){
    return minN + rand() % (maxN + 1 - minN);
}
double double_rand( double min, double max ) {
    double scale = abs(rand()) / (double) RAND_MAX; /* [0, 1.0] */
    return min + scale * ( max - min );      /* [min, max] */
}
void readfile(int qwer,int seed,double vol) {
	srand(seed);
	string s=tenfile[qwer];
	cout<<s<<": "<<seed<<endl;
	ifstream itt(s.c_str());
	if (itt.fail())
	cout << "Failed to open this file!!" <<endl;
    int sodong=0;string q;char c;
	while (itt) {
	    int dem=0;
        while (itt && c != '\n') {
        itt.get(c);
        dem++;
        }
    if(dem>=2) sodong += 1;
    itt.get(c);
    }
    n=sodong-1;
	collect.pop.clear();
	sensors.clear();
	sensors.resize(n+1);
	ifstream it(s.c_str());
	it>>sensors[0].x>>sensors[0].y;
	for(int i=1;i<=n;i++) {
		it>>sensors[i].x>>sensors[i].y>>sensors[i].p>>sensors[i].energy;
		sensors[i].w=sensors[i].energy/sensors[i].p;
	}
	for(int i=0;i<n+1;i++) {
		for(int j=0;j<n+1;j++) {
			dis[i][j]=sqrt(pow(sensors[i].x-sensors[j].x,2)+pow(sensors[i].y-sensors[j].y,2));
			if(i==j) {dis[i][j]=999999.0;continue;}
		}
	}
	it.close();
}
double tinhfitness1(Individual ob) {
	double tongkhoangcach=0;
	double f1=0,f2=0;
	for(int i=1;i<ob.path.size()-1;i++) {
		tongkhoangcach+=(dis[ob.path[i]][ob.path[i-1]])/v;
		f1+=tongkhoangcach/(sensors[ob.path[i]].w);
	}
	tongkhoangcach=0;
	for(int i=1;i<ob.path.size()-1;i++) {
		tongkhoangcach+=(dis[ob.path[i]][ob.path[i-1]])/v;
        f2+=abs(tongkhoangcach/(sensors[ob.path[i]].w)-f1/n);
	}
	return f1/2+f2/2;
}
vector<Individual> pmx(int a, int b) {
	int j=random(1,n-1);
	int k=random(j+1,n);
	Individual child1,child2;
	child1.path.resize(n+2);
	child2.path.resize(n+2);
	vector<int> cotronga(n+2,-1),cotrongb(n+2,-1),thea(n+2,0),theb(n+2,0);
	for(int i=0;i<=n;i++) {
		thea[i]=i;
		theb[i]=i;
	}
	for(int i=j;i<=k;i++) {
		cotronga[collect.pop[a].path[i]]=i;
		cotrongb[collect.pop[b].path[i]]=i;
	}
	for(int i=j;i<=k;i++) {
		if(cotrongb[collect.pop[a].path[i]]>=0) continue;
		int vitri=i;
		while(cotronga[collect.pop[b].path[vitri]]>=0) {
			vitri=cotronga[collect.pop[b].path[vitri]];
		} 
		thea[collect.pop[b].path[vitri]]=collect.pop[a].path[i];
    }
    for(int i=j;i<=k;i++) {
		if(cotronga[collect.pop[b].path[i]]>=0) continue;
		int vitri=i;
		while(cotrongb[collect.pop[a].path[vitri]]>0) {
			vitri=cotrongb[collect.pop[a].path[vitri]];
		} 
		theb[collect.pop[a].path[vitri]]=collect.pop[b].path[i];
    }
    child1=collect.pop[a],child2=collect.pop[b];
    for(int i=j;i<=k;i++) {
    	child1.path[i]=collect.pop[b].path[i];
    	child2.path[i]=collect.pop[a].path[i];
	}
	for(int i=1;i<=n;i++) {
		if(i>=j&&i<=k) continue;
    	child1.path[i]=thea[child1.path[i]];
    	child2.path[i]=theb[child2.path[i]];
	}
	child1.fitnesspath=tinhfitness1(child1);
	child2.fitnesspath=tinhfitness1(child2);
	vector<Individual> kqq;
	kqq.push_back(child1);
	kqq.push_back(child2);
	return kqq;
}
Individual CIM(Individual a) {
	int p=random(2,n-1);
	Individual child;
	child.path.push_back(0);
	for(int i=p;i>=1;i--) {
		child.path.push_back(a.path[i]);
	}
	for(int i=n;i>p;i--) {
		child.path.push_back(a.path[i]);
	}
	child.path.push_back(0);
	child.fitnesspath=tinhfitness1(child);
	return child;
}
Individual fswap(Individual a) {
	int p=random(1,n);
	int q=random(1,n);
	while(q==p) q=random(1,n);
	int t;
	t=a.path[p];
    a.path[p]=a.path[q];
	a.path[q]=t;
	a.fitnesspath=tinhfitness1(a);
	return a;
}
vector<Individual> singlepoint(int a, int b) {
	int p=random(2,n-1);
	Individual child1,child2;
	vector<bool> ks(n+1,1);
	child1.path.push_back(0);
	for(int i=1;i<=p;i++) {
		child1.path.push_back(collect.pop[a].path[i]);
		ks[collect.pop[a].path[i]]=0;
	}
	for(int i=1;i<=n;i++) {
		if(ks[collect.pop[b].path[i]]==1) child1.path.push_back(collect.pop[b].path[i]);
	}
	child1.path.push_back(0);
	vector<bool> ks2(n+1,1);
	child2.path.push_back(0);
	for(int i=1;i<=p;i++) {
		child2.path.push_back(collect.pop[b].path[i]);
		ks2[collect.pop[b].path[i]]=0;
	}
	for(int i=1;i<=n;i++) {
		if(ks2[collect.pop[a].path[i]]==1) child2.path.push_back(collect.pop[a].path[i]);
	}
	child2.path.push_back(0);
	child1.fitnesspath=tinhfitness1(child1);
	child2.fitnesspath=tinhfitness1(child2);
	vector<Individual> kqq;
	kqq.push_back(child1);
	kqq.push_back(child2);
	return kqq;
} 
Individual output(Individual ob) {
	double dichuyen=0,sac=0,cc=0,sonutchet=0;
	vector<int> t1(n+1,0);
	vector<int> t2(n+1,0);
	vector<bool> charged(n+2,0);
	for(int i=1;i<ob.path.size()-1;i++) {
		dichuyen+=dis[ob.path[i]][ob.path[i-1]]/v;
		if(dichuyen+sac>T) {dichuyen=T-sac;break;};
		charged[ob.path[i]]=1;
		if(sensors[ob.path[i]].energy-(dichuyen+sac)*sensors[ob.path[i]].p<Emin) t1[ob.path[i]]=1;
		if(sensors[ob.path[i]].energy-(dichuyen+sac)*sensors[ob.path[i]].p+ob.chargetime[i]*(U-sensors[ob.path[i]].p)>Emax) {
			double ttt=-Emax+sensors[ob.path[i]].energy-(dichuyen+sac)*sensors[ob.path[i]].p+ob.chargetime[i]*(U-sensors[ob.path[i]].p);
			ttt/=abs(U-sensors[ob.path[i]].p);
			ob.tsac-=(ttt);
			ob.chargetime[i]-=ttt;
		}
		if(t1[ob.path[i]]==1) {ob.chargetime[i]=0;continue;} else {
		sac+=ob.chargetime[i];
		if(dichuyen+sac>T) {
		double hh=dichuyen+sac-T;
		ob.chargetime[i]-=hh;
		sac-=hh;
		break;
	    }
	    }
	    if(i==ob.path.size()-2) dichuyen+=dis[ob.path[i]][0]/v;
	}
	for(int i=1;i<ob.path.size()-1&&charged[ob.path[i]]==1;i++) {
		if(sensors[ob.path[i]].energy-(T)*sensors[ob.path[i]].p+ob.chargetime[i]*U<Emin) t2[ob.path[i]]=1; 
	}
	for(int i=1;i<=n;i++) {
		if(charged[i]==1) continue;
		if(sensors[i].energy-(T)*sensors[i].p<Emin) t2[i]=1;
	}
	for(int i=1;i<=n;i++) {
		if(t1[i]==1||t2[i]==1) sonutchet++;
    }
    sonutchet+=abs(dichuyen*Pm)/Emc;
    ob.fitness=sonutchet;
    ob.tdichuyen=dichuyen;
    ob.tsac=sac;
    return ob;
}
Individual init_fitness(Individual ob) {
	double dichuyen=0,sac=0,cc=0,sonutchet=0;
	Individual xx=ob;
	vector<int> t1(n+1,0);
	vector<int> t2(n+1,0);
	vector<bool> charged(n+2,0);
	for(int i=1;i<ob.path.size()-1;i++) {
		dichuyen+=dis[ob.path[i]][ob.path[i-1]]/v;
		if(dichuyen+sac>T) {dichuyen=T-sac;break;};
		charged[ob.path[i]]=1;
		if(sensors[ob.path[i]].energy-(dichuyen+sac)*sensors[ob.path[i]].p<Emin) t1[ob.path[i]]=1;
		if(sensors[ob.path[i]].energy-(dichuyen+sac)*sensors[ob.path[i]].p+ob.chargetime[i]*(U-sensors[ob.path[i]].p)>Emax) {
			double ttt=-Emax+sensors[ob.path[i]].energy-(dichuyen+sac)*sensors[ob.path[i]].p+ob.chargetime[i]*(U-sensors[ob.path[i]].p);
			ttt/=abs(U-sensors[ob.path[i]].p);
			ob.tsac-=(ttt);
			ob.chargetime[i]-=ttt;
		}
		if(t1[ob.path[i]]==1) {ob.chargetime[i]=0;continue;} else {
		sac+=ob.chargetime[i];
		if(dichuyen+sac>T) {
		double hh=dichuyen+sac-T;
		ob.chargetime[i]-=hh;
		sac-=hh;
		break;
	    }
	    }
	    if(i==ob.path.size()-2) dichuyen+=dis[ob.path[i]][0]/v;
	}
	for(int i=1;i<ob.path.size()-1&&charged[ob.path[i]]==1;i++) {
		if(sensors[ob.path[i]].energy-(T)*sensors[ob.path[i]].p+ob.chargetime[i]*U<Emin) t2[ob.path[i]]=1; 
	}
	for(int i=1;i<=n;i++) {
		if(charged[i]==1) continue;
		if(sensors[i].energy-(T)*sensors[i].p<Emin) t2[i]=1;
	}
	for(int i=1;i<=n;i++) {
		if(t1[i]==1||t2[i]==1) sonutchet++;
    }
    sonutchet+=abs(dichuyen*Pm)/Emc;
    ob.fitness=sonutchet;
    ob.chargetime=xx.chargetime;
    ob.tsac=xx.tsac;
    return ob;
}
void Create_Pop () {
collect.pop.clear();
	Individual individual;
	vector<int> base;
	base.push_back(0);
	for(int i=1;i<=n;i++) {
		base.push_back(i);
	}
	base.push_back(0);
	for(int i=0;i<POP_SIZE;i++) {
		random_shuffle(base.begin()+1,base.end()-1);
		individual.path=base;
		individual.fitnesspath=tinhfitness1(individual);
		collect.pop.push_back(individual);
	}
}
void reproduction(){
	Individual parent1, parent2;
	int m1,m2,p1,p2;
	while (collect.pop.size() < POP_SIZE*2){
		m1=random(0,POP_SIZE-1);
		do{
		p1=random(0,POP_SIZE-1);
		} while(p1==m1);
		parent1=collect.pop[p1];
		parent2=collect.pop[m1];
		if (double_rand(0, 1) <= 0.5){
		vector<Individual> child ;
		if(double_rand(0,1)<0.6) child=pmx(p1,m1);
		else child=singlepoint(m1,p1);
		if(double_rand(0,1)<0.1) {
			if(double_rand(0,1)<0.6) child[0]=CIM(child[0]);
			else child[0]=fswap(child[0]);
		}
		if(double_rand(0,1)<0.1) {
			if(double_rand(0,1)<0.6) child[1]=CIM(child[1]);
			else child[1]=fswap(child[1]);
		}
		collect.pop.push_back(child[0]);
		collect.pop.push_back(child[1]);
		}
	}
}
Individual solve() {
	Create_Pop();
	int gen = 0;	
	while (gen++ < GENERATIONS){
		for(int i=0;i<collect.pop.size();i++) {
			for(int j=i+1;j<collect.pop.size();j++) {
				if(collect.pop[j].fitnesspath<collect.pop[i].fitnesspath) swap(collect.pop[i],collect.pop[j]);
			}
		}
		collect.pop.resize(POP_SIZE);
		collect.best=collect.pop[0];
		reproduction();
	}
	collect.pop.resize(POP_SIZE);
	return collect.best;
}
void init2() {
	double tongthoigiansac=solution.tsac;
	double conlai;
	for(int i=0;i<collect.pop.size();i++) {
		collect.pop[i]=solution;
		collect.pop[i].chargetime.resize(n+2);
		conlai=tongthoigiansac;
		for(int j=0;j<collect.pop[i].chargetime.size()-2;j++) {
           double cantren=conlai;
           if(cantren>(Emax-Emin)/(U-sensors[collect.pop[i].path[j]].p)) cantren=(Emax-Emin)/(U-sensors[collect.pop[i].path[j]].p);
           if(cantren<0) cantren=0;
		   double t=double_rand(0,cantren);
           collect.pop[i].chargetime[j]=t;
           conlai-=(t);
		}
		collect.pop[i]=init_fitness(collect.pop[i]);
	}
	return ;
}
Individual mut2(Individual a) {
	int i=random(1,a.chargetime.size()-2);
	int j;
	do {
		j=random(1,a.chargetime.size()-2);
	} while (i==j);
	double dieuchinh1=(double_rand(0,abs((Emax-Emin)/(U-sensors[a.path[i]].p)-a.chargetime[i]))),dieuchinh2=abs(double_rand(0,a.chargetime[j]));
	if(dieuchinh1>dieuchinh2) dieuchinh1=dieuchinh2;
	a.chargetime[i]+=dieuchinh1;
	a.chargetime[j]-=dieuchinh1;
	a=init_fitness(a);
	return a;
}
vector<Individual> SPAH (Individual a, Individual b) {
	Individual child1=a,child2=b;
	int point=random(1,a.path.size()-2);
	double gr=a.tsac;
	double beta=double_rand(-0.5,0.5);
	for(int i=point+1;i<=a.path.size()-2;i++) {
		child1.chargetime[i]=b.chargetime[i];
		child2.chargetime[i]=a.chargetime[i];
	}
	child1.chargetime[point]=abs((1-beta)*a.chargetime[point]+beta*(b.chargetime[point]));
	child2.chargetime[point]=abs((beta)*a.chargetime[point]+(1-beta)*(b.chargetime[point]));
	double t1=0,t2=0;
	for(int i=1;i<=a.chargetime.size()-2;i++) {
		t1+=child1.chargetime[i];
		t2+=child2.chargetime[i];
	}
	if(t1>gr) {
		double heso=t1/gr;
		for(int i=1;i<=child1.chargetime.size()-2;i++) child1.chargetime[i]/=heso;
	} 
	else {
		double plus=gr-t1;
		for(int i=1;i<=child1.chargetime.size()-2&&plus>0;i++) {
			double t=double_rand(0,abs((Emax-Emin)/(U-sensors[child1.path[i]].p)-child1.chargetime[i]));
			if(t>plus) {
				t=plus;
			}
			child1.chargetime[i]+=t;
			plus-=t;
		}
	}
	if(t2>gr) {
		double heso=t2/gr;
		for(int i=1;i<=child2.chargetime.size()-2;i++) child2.chargetime[i]/=heso;
	} 
	else {
		double plus=gr-t2;
		for(int i=1;i<=child2.chargetime.size()-2&&plus>0;i++) {
			double t=double_rand(0,abs((Emax-Emin)/(U-sensors[child2.path[i]].p)-child2.chargetime[i]));
			if(t>plus) {
				t=plus;
			}
			child2.chargetime[i]+=t;
			plus-=t;
		}
	}
	child1=init_fitness(child1);
	child2=init_fitness(child2);
	vector<Individual> q;
	q.push_back(child1);
	q.push_back(child2);
	return q;
}
void reproduction2() {
	Individual parent1, parent2;
	int m1,p1;
	for(int i=0;i<POP_SIZE/2;i++) {
		m1=random(0,POP_SIZE-1);
		do{
		p1=random(0,POP_SIZE-1);
		} while(p1==m1);
		parent1=collect.pop[p1];
		parent2=collect.pop[m1];
		if (double_rand(0, 1) <= 1){
		vector<Individual> child = SPAH(parent1, parent2);
		if(double_rand(0,1)<MUTATION_RATE2) {
			child[0]=mut2(child[0]);
		}
		if(double_rand(0,1)<MUTATION_RATE2) {
			child[1]=mut2(child[1]);
		}
		collect.pop.push_back(child[0]);
		collect.pop.push_back(child[1]);
	}
    }
}
Individual solve2() {
	int gen = 0;	
	while (gen++<GENERATIONSt){
		for(int i=0;i<collect.pop.size();i++) {
			for(int j=i+1;j<collect.pop.size();j++) {
				if(collect.pop[j].fitness<collect.pop[i].fitness) swap(collect.pop[i],collect.pop[j]);
			}
		}
		collect.pop.resize(POP_SIZE);
		collect.best=collect.pop[0];
		reproduction2();
	}
	return collect.best;
}
void init_T(double uuu) {
	double pp=(double)n/(uuu*uuu);
	double tmove=(n*sqrt(abs(2/pp)))/v;
	T=(Emc-Pm*tmove)/5+tmove;
}
int main() { 
    fstream fout,fout2,fout3;                      
    DIR *dpdf;
    struct dirent *epdf;
    dpdf=opendir("dulieu");
    int ks=0,i=1;
    if (dpdf !=NULL){
         while (epdf=readdir(dpdf)){
         	if(ks++<2) continue;
         	tenfile[i]="dulieu/"+(string) epdf->d_name;
         	i++;
        }
    }
    closedir(dpdf);
    T=20000;
    string fileghi="GACS/seed=1-10.txt";
    string fileghi2="GACS/thoigiandichuyenSeed=1-10.txt";
    fout.open(fileghi.c_str(),ios::out);
	fout2.open(fileghi2.c_str(),ios::out);
	double E[5]={108000,162000,216000,270000,324000};
	for(double vol=1;vol<=3;vol=vol+0.5) {
     	U=4*(vol);
        for(int i=1;i<=1;i++) {
            for(int seed=1;seed<=10;seed++) {
            readfile(i,seed,8.0);
            solution = solve();
            solution.tdichuyen=0;
            for(int i=1;i<solution.path.size()-1;i++) {
		    solution.tdichuyen+=dis[solution.path[i]][solution.path[i-1]]/v;
	        }
	        solution.tdichuyen+=dis[solution.path[solution.path.size()-2]][0]/v;
	        solution.tsac=min((Emc-solution.tdichuyen*Pm)/U,T-solution.tdichuyen);
	        if(solution.tsac<0) solution.tsac=0;
            init2();
		    collect.best=collect.pop[0];
            solution=solve2();
            solution=output(solution);
            solution.check();
            fout<<(int)solution.fitness<<'\t';
            fout2<<solution.tdichuyen<<'\t';
            }
            fout<<'\n';
            fout2<<'\n';
        }
    }
    fout.close();
    fout2.close();
    return 0;
}
