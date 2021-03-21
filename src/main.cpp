#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include<cstring>
#include <ctime>
#include <cstdlib>
#include <sys/socket.h> 
#include <unistd.h>
// #include<windows.h>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

#define DEBUG
// #define CHROM
#define DRAW

Mat srcImage;
vector<Point> drawPoint;
 
const int CUS_MAXN = 100;
const double pi = acos(-1);
const int VEH_MAX = 10;
const int INF = 1 << 30;
 
int generation = 100000;	//代数 
double mutate_p_swap = 600;   //1---1000 
double mutate_p_rev = 20; //1---1000
double mutate_p_shift = 800;
const int PRO = 10000; 	//变异概率的分母 
// const string NAME="A-n33-k6.vrp";
 
class Customer {
public:
	int x, y;
	int demand;
	double polar;
	double distance;
	Customer();
	Customer(int, int, int);
 
	bool operator<(const Customer&) const;
 
	// double operator-(const Customer&);
};
 
class Vehicle
{
public:
	vector<int> cus_vec;
	double x;
	double y;
	double polar;
	double distance;
	double length;
	int cap;
	int cap_remain;

	void clear();
	bool push(int c); // 输入的是顾客的index
 
	int get_num(); // 返回顾客的数列
	void get_coor(); // 得到这辆货车的新的坐标
 
	bool operator<(const Vehicle& b) const;
	void optimate();
	Vehicle();
	~Vehicle();
};
 
class Individual
{
public:
	int chromosome[CUS_MAXN]; // 染色体序列
	Vehicle veh[VEH_MAX];
	
	double fitness;		//适应度
	double unfitness;
 
	void update();
	bool is_satisfy();
	double get_ufit();
	double get_fit();
	
	Individual& operator = (const Individual &b); 
	Individual();
	~Individual();
};
 
 
 
int cus_num; // 顾客数
int veh_num; // 货车数
int opt_dis; // 最优解
int best_ind; // 当前最优解
int CAP; // 每辆货车的容量
 
double tightness; // 宽松度
double Rc, Rd; // ~~~
 
Individual best;//当前得到的最优解 
const int popul_size = 30; //每一代人口数量
double dis[CUS_MAXN][CUS_MAXN]; //距离矩阵

Customer cus[CUS_MAXN];
Individual popul[popul_size];

	
//初始解 
void ini_popul();
//binary tournament 
int selection(int tournament_size);
//交叉变异 
Individual gene_two_point_crossover(int curr1,int curr2,Individual &ans1, Individual &ans2); 
//进化换代 
void evolution(); 
 
//突变 
void gene_swap(Individual &ans);
void gene_shift(Individual &ans);
void gene_reverse(Individual &ans);
 
double get_polar(double, double); // 角度

double get_dis_of_cus(const Customer&, const Customer&);
double get_dis(double, double);

void get_d();
 
bool ini_data();

// void draw_picture();
 
 
//输入坐标返回极坐标
double get_polar(double x, double y) {
	double polar;
	if (x > 0 && y >= 0) {
		return polar = atan(1.0 * y / x);
	}
	if (x == 0) {
		if (y > 0) return polar = pi / 2;
		else return polar = 3 * pi / 2;
	}
	if (x < 0 && y >= 0) {
		return polar = pi - atan(-1.0 * y / x);
	}
	if (x < 0 && y < 0) {
		return polar = pi + atan(1.0 * y / x);
	}
	if (x > 0 && y < 0) {
		return polar = 3 * pi / 2 + atan(-1.0 * x / y);
	}
}

double get_dis(double x, double y) {
	return sqrt(x * x + y * y);
}
 
double get_dis_of_cus(const Customer& a, const Customer& b) {
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
 
void get_d() {
	for (int i = 0; i <=cus_num; i++) {
		for (int j = 0; j <= cus_num; j++) {
			if (i==j){
				dis[i][j]=0;
			}
			else if (i==0){
				dis[i][j]=get_dis(cus[j-1].x,cus[j-1].y);
			}
			else if (j==0){
				dis[i][j]=get_dis(cus[i-1].x,cus[i-1].y);
			}
			dis[i][j] = get_dis_of_cus(cus[i-1], cus[j-1]);
		}
	}
}
 
bool ini_data() {
	ifstream in("../data/A-n33-k5.vrp");
	//ifstream in(NAME);
	if (!in.is_open()) return false;
	
	string s;
	if (getline(in, s)) {
		cout << s << endl;
		getline(in, s);
		istringstream is(s);
		is >> veh_num >> opt_dis;
		// cout << "veh_num = " << veh_num << endl;
		// cout << "opt_dis = " << opt_dis << endl;
		getline(in, s);
		is.str(s); is.clear();
		cout << is.str() << endl;
		getline(in, s);
		is.str(s); is.clear();
		is >> cus_num;
		cout << "cus_num = " << cus_num << endl;

		getline(in, s);
		is.str(s); is.clear();
		is >> CAP;
		// cout << "CAP = " << CAP << endl;
		// cout << "popul_size = " << popul_size << endl;
		for (int i=0 ; i<popul_size ; i++){
			for (int j=0 ; j<veh_num ; j++){
				popul[i].veh[j].clear();	
			}
		}
		getline(in, s);
		int t;
		//获取相对坐标
		int x,y;
		getline(in, s);
		is.str(s); is.clear();
		is >> t;
		is >> x >> y;
		drawPoint.push_back(Point(x*7, y*7));
		for (int i = 0; i < cus_num-1; i++) {
			getline(in, s);
			is.str(s); is.clear();
			is >> t;
			is >> cus[i].x >> cus[i].y; //求出各个点的相对坐标
			drawPoint.push_back(Point(cus[i].x*7, cus[i].y*7));
			cus[i].x -= x;
			cus[i].y -= y;
			cus[i].polar = get_polar(cus[i].x, cus[i].y);
			#ifdef DEBUG
			cout << t << "\t" << cus[i].x << "\t" << cus[i].y << endl;
			#endif
			cus[i].distance = get_dis(cus[i].x, cus[i].y);
			
		}
 
		getline(in, s);
		double sum_all = 0;
		getline(in, s);
		is.str(s); is.clear();
		is >> t;
		is >> t;
		for (int i = 0; i < cus_num-1; i++) {
			getline(in, s);
			is.str(s); is.clear();
			is >> t;
			is >> cus[i].demand;
			sum_all += cus[i].demand;
			#ifdef DEBUG
			cout << i+1 << "\t" << cus[i].demand << endl;
			#endif
		}
		tightness = 1.0 * sum_all / (veh_num * CAP);
		if (tightness >= 0.97) {
			Rc = 0.6;
		}
		else if (tightness >= 0.94){
			Rc=0.75;
		}
		else Rc = 0.9;
	}
	in.close();
	cus_num-=1;
	sort(cus , cus + cus_num);
	get_d();
	cout << "data initial successful!" << endl;
	//Sleep(5000);
	return 1;
}
 
 
 
void ini_popul(){
	//srand((unsigned)time(NULL)); 
	
	int num;
	int p = 0;
	//cout<<popul_size<<endl;
	
	//for(int i=0;i<popul_size;){
	for(int i=0;i<popul_size/2;)
	{
		num = rand()%cus_num;
		//cout<<num<<endl;
		p = 0;
		bool f = false;
		for (int j=0 ; j<veh_num ; j++)
		{
			popul[i].veh[j].clear();
		}
		for (int j=0 ; j<cus_num;)
		{
			int k = (num+j)%cus_num;
			// cout<<k<<veh_num<<endl;
			// cout<<"v[i].cap"<<popul[i].veh[p].cap_remain<<endl;
			if (p >= veh_num)
			{
				//cout<<"veh_num"<<p;
				//p=0;
				f=false;
				break;
			}
			if (popul[i].veh[p].push(k)) 
			{
				//cout<<"remain: "<<popul[i].veh[p].cap_remain<<endl; 
				popul[i].chromosome[k] = p;
				if (p == veh_num-1){
					f = true;
				}
				j++;
			}
			else p++;
		}
		//cout<<p<<endl;
		if (f){
			popul[i].update();
			// cout<<"kkkkkkkkkk"<<endl;
			i++;
			// cout<<i<<endl;
		}
		
	//	cout<<"kk"<<endl;
	}
	
	
	/*
	double ave_demand = 0;
	for (int i = 0; i < cus_num; i++)
		ave_demand += cus[i].demand;
	ave_demand /= veh_num;
	int seed[VEH_MAX] = {0};
	int index[VEH_MAX];
	int demands = 0;
	int k = 0;
	// 分类，index[] 是种子对应的顾客编号。
	for (int i = 0; i < cus_num; i++) {
		demands += cus[i].demand;
		if (demands > ave_demand && k < veh_num) {
			demands = cus[i].demand;
			++k;
		}
		// 每个种子为最大的距离
		if (cus[i].distance > seed[k]) {
			seed[k] = cus[i].distance;
			index[k] = i;
		}
	}
	int seed_index[2][CUS_MAXN] = {0}; // 每一个顾客的两个最优种子的编号
	int belong_index[2][CUS_MAXN]; // 种子对应车辆的编号
	for (int i = 0; i < cus_num; i++) {
		double dif0 = 1000000, dif1 = 1000000;
		for (int j = 0; j < veh_num; j++) {
			double dif = get_dis_of_cus(cus[i], cus[index[j]]);
			if (dif < dif0) { // 距离小于第一个种子
				dif1 = dif0;
				seed_index[1][i] = seed_index[0][i];
				seed_index[0][i] = index[j];
				belong_index[1][i] = belong_index[0][i];
				belong_index[0][i] = j;
				dif0 = dif;
			} else if (dif < dif1) { // 距离小于第二个种子
				seed_index[1][i] = index[j];
				belong_index[1][i] = j;
				dif1 = dif;
			}
		}
	}
	for (int i = 0; i < cus_num; i++) {
		if (seed_index[1][i] == 0) {
			seed_index[1][i] = seed_index[0][i];
			belong_index[1][i] = belong_index[0][i];
		}
	}
	// 获得两个种子的概率
	int pro[2][CUS_MAXN];
	for (int i = 0; i < cus_num; i++) {
		for (int j = 0; j < 2; j++) {
			pro[j][i] = dis[0][i + 1] + dis[i + 1][seed_index[j][i] + 1]
			            - dis[0][seed_index[j][i] + 1];
		}
		for (int j = 0; j < 2; j++) {
			pro[j][i] = (int)(100.0 * pro[j][i] / (pro[0][i] + pro[1][i]));
		}
	}
	srand((unsigned)time(0));
	// 根据概率获取15个初始解
	// int init_num = 15;
	for (int i = 0 ; i < popul_size; i++) {
	//for (int i = popul_size/2 ; i < popul_size; i++) {
		bool isV[VEH_MAX] = {0};
		int cnt = 0;
		for (int j = 0; j < cus_num; j++) {
			int random = rand() % 100;
			if (random > pro[0][j]) {
				// popul[i].veh[belong_index[0][j]].push(j - 1);
				popul[i].chromosome[j] = belong_index[0][j];
			} else {
				// popul[i].veh[belong_index[1][j]].push(j - 1);
				popul[i].chromosome[j] = belong_index[1][j];
			}
			if (!isV[popul[i].chromosome[j]]) {
				++cnt;
				isV[popul[i].chromosome[j]] = true;
			}
		}
		if (cnt == veh_num)
			popul[i].update();
		// ++init_num;
	}
	
	*/
}
 
 
int selection(int tournament_size){
	int selected = rand() % popul_size;
	for (int i = 1; i < tournament_size; i++){
		int r = rand() % popul_size;
		if (popul[r].fitness > popul[selected].fitness){
			selected = r;
		}
	}
	return selected;
}
 
 
void evolution(){
	for(int T=0;T<generation;T++){
		int father1=selection(2);
		int father2=selection(2);
		Individual son[2];
		Individual tmp;
	//	cout<<father1<<" "<<father2<<endl; 
		tmp=gene_two_point_crossover(father1,father2,son[0],son[1]);
	//	cout<<"xxx"<<endl;
		bool found=0;
		for(int j=0;j<popul_size;j++){
			if(popul[j].fitness >= tmp.fitness && popul[j].unfitness >= tmp.unfitness){
				popul[j]=tmp;found=1;break;
			}
		}
		if(!found) for(int j=0;j<popul_size;j++){ 
			if(popul[j].fitness >= tmp.fitness && popul[j].unfitness <= tmp.unfitness){
				popul[j]=tmp;found=1;break;
			}
		}
		if(!found) for(int j=0;j<popul_size;j++){
			if(popul[j].fitness <= tmp.fitness && popul[j].unfitness >= tmp.unfitness){
				popul[j]=tmp;found=1;break;
			}
		} 
		if(!found) for(int j=0;j<popul_size;j++){	
			if(popul[j].fitness <= tmp.fitness && popul[j].unfitness <= tmp.unfitness){
				popul[j]=tmp;found=1;break;
			}
		}
	}
}
 
 
Individual gene_two_point_crossover(int curr1,int curr2,Individual &ans1, Individual &ans2){
	//int T=3; 
	//while(T--){
		ans1=popul[curr1];
		ans2=popul[curr2];
		
	//while(1){
		//cout<<curr1<<" cur"<<curr2<<endl;
		//curr1=selection(2);
		//curr2=selection(2);
	
	
		/*
		for (int i=0;i<cus_num;i++){
			cout<<ans1.chromosome[i]<<' ';
		}
		cout<<endl;
		for (int i=0;i<cus_num;i++){
			cout<<ans2.chromosome[i]<<' ';
		}
		cout<<endl;*/
	//交叉 
		int pt1 = rand() % cus_num;
		int pt2 = rand() % cus_num;
		
		if (pt1 > pt2) swap(pt1, pt2);
		//cout<<pt1<<" "<<pt2<<endl;
		
		
		for(int i=pt1;i<=pt2;i++)
			swap(ans1.chromosome[i], ans2.chromosome[i]);
		/*
		for (int i=0;i<cus_num;i++){
			cout<<ans1.chromosome[i]<<' ';
		}
		cout<<endl;
		for (int i=0;i<cus_num;i++){
			cout<<ans2.chromosome[i]<<' ';
		}
		
		cout<<endl;
		cout<<endl;
		Sleep(5000);*/
		//变异 
		//cout<<pt1<<" x "<<pt2<<endl;
		
		int p_mutate = (rand()%PRO);//变异概率范围：0-0.01
		if (p_mutate < mutate_p_swap){
			gene_swap(ans1);
			gene_swap(ans2);
		}	
	
		double p_mutate2 = (rand()%PRO);
	//	cout<<"rand"<<p_mutate2<<endl;
		if (p_mutate2 < mutate_p_shift){
			gene_shift(ans1);
			gene_shift(ans2);
		}
		ans1.update();
		ans2.update();
		if (ans1.fitness>=ans2.fitness&&ans1.unfitness>=ans2.unfitness){
			return ans1;
		}
		else return ans2;
		/*
		 
		for (int i=0;i<cus_num;i++){
			cout<<ans1.chromosome[i]<<' ';
		}
		cout<<endl;
		for (int i=0;i<cus_num;i++){
			cout<<ans2.chromosome[i]<<' ';
		}
		cout<<endl;
		cout<<endl;
		*/
		
		/*
		if (ans1.is_satisfy()) {
			ans1.update();	
			if (ans2.is_satisfy()) {
				ans2.update();
				if (ans1.fitness>ans2.fitness){
					return ans1;
				}
				else return ans2;			
			}
			else return ans1;
		}
		if (ans2.is_satisfy()) {
			ans2.update();
			return ans2;	
		}*/
		//Sleep(10000);
	//}
	//return ans1;
}
 
 
void gene_swap(Individual &ans){
	//int len = cus_num;
	int pt1 = rand() % cus_num;
	int pt2 = rand() % cus_num;
	//if (pt1 > pt2) swap(pt1, pt2);
	//cout<<pt1<<" swap "<<pt2<<endl;
 
	swap(ans.chromosome[pt1], ans.chromosome[pt2]);
//	ans.fitness = ans.compute_fitness(ans);//gai
}
void gene_shift(Individual &ans){
	//int len = cus_num;
	int pt1 = rand() % cus_num;
	//int pt2 = rand() % cus_num;
	//if (pt1 > pt2) swap(pt1, pt2);
	if (pt1%2==0&&pt1+1<cus_num){
		swap(ans.chromosome[pt1], ans.chromosome[pt1+1]);	
	}
	if (pt1%2==1&&pt1-1>=0){
		swap(ans.chromosome[pt1], ans.chromosome[pt1-1]);	
	}
//	ans.fitness = ans.compute_fitness(ans);//gai
}
 
 
 
void gene_reverse(Individual &ans){
	int len = cus_num;
	int pt1 = rand() % len; //pt1: [0, cus_num-1]
	int pt2 = rand() % len; //pt2: [0, cus_num-1]
	if (pt1 > pt2) swap(pt1, pt2);
	
	for (int i = pt1, j = pt2; i <= j; i++, j--){
		swap(ans.chromosome[i], ans.chromosome[j]);
	}
	//ans.fitness = ans.compute_fitness(ans);// gai
}
 
 

 
Customer::Customer() {
	x = y = demand = 0;
}
 
Customer::Customer(int x, int y, int d) {
	this->x = x;
	this->y = y;
	this->demand = d;
	polar = get_polar(x, y);
	distance = get_dis(x, y);
}
 
bool Customer::operator<(const Customer & c) const {
	if (polar == c.polar) return distance < c.distance;
	return polar < c.polar;
}
 
 
 
 
Vehicle::Vehicle() {
	cap = CAP;
	cap_remain = CAP;
}
 
 
Vehicle::~Vehicle() {
}
 
void Vehicle::clear() {
	cus_vec.clear();
	cap_remain = CAP;
}
 
bool Vehicle::push(int c) {
	if (cap_remain >= cus[c].demand) {
		cus_vec.push_back(c);
	//	cout<<c<<" c num and demand "<<cus[c].demand<<endl;;
		cap_remain -= cus[c].demand;
	//	cout<<"cap "<<cap_remain<<endl;
	//	cout<<"custom "<<c<<endl;
	//	cout<<"push true"<<endl;
		return 1;
	}
	else {
		// 根据Rc判断是否接受解
		if (cap_remain * 1.0 / cus[c].demand > Rc) {
			cus_vec.push_back(c);
			cap_remain -= cus[c].demand;
			return 1;
		}
		return 0;
	}
}
 
int Vehicle::get_num() {
	return cus_vec.size();
}
 
void Vehicle::get_coor() {
	double sumx = 0, sumy = 0;
	int l = get_num();
	for (int i = 0; i < l; i++) {
		sumx += cus[cus_vec[i]].x;
		sumy += cus[cus_vec[i]].y;
	}
	x = sumx / l;
	y = sumy / l;
	polar = get_polar(x, y);
	distance = get_dis(x, y);
}
 
bool Vehicle::operator<(const Vehicle& v)const {
	if (polar != v.polar) {
		return polar < v.polar;
	}
	else return distance < v.distance;
}
 
void Vehicle::optimate() {
	bool isv[cus_num];
	int l=cus_vec.size();
	memset(isv,0,sizeof(isv));
	vector<int> v;
	double len=0;
	int now=0;
	double len_min=INF;
	int t;
	int min_p=0;
//	cout<<"veh"<<endl;
	for (int k=0;k<l;k++){
		len_min=INF;
		for (int i=0;i<l;i++){
			t=cus_vec[i];
			if (!isv[t]&&dis[now][t+1]<len_min){
				min_p=t;
				len_min=dis[now][t+1];
			}
		}
//		cout<<len_min<<endl;
		len+=len_min;
		isv[min_p]=1;
		now=min_p+1;
		v.push_back(min_p);	
	}
//	cout<<"end"<<endl;
	len+=dis[min_p+1][0];
	cus_vec=v;
	length=len;
	return ;
	/*
	for (int i = 0; i < cus_vec.size() - 1; i++)
	{
		double i_x = cus[cus_vec[i]].x;
		double i_y = cus[cus_vec[i]].y;
		double i1_x = cus[cus_vec[i + 1]].y;
		double i1_y = cus[cus_vec[i + 1]].y;
		if (get_dis(i1_x, i1_y) < get_dis(i_x, i_y) / 2)
		{
			swap(cus_vec[i], cus_vec[i + 1]);
		}
	}*/
}
 
Individual::Individual() {
	fitness=0;
	unfitness=INF;
	//cout<<"ini_indi"<<endl;
	for (int i=0;i<veh_num;i++){
		veh[i].clear();
	}
}
Individual::~Individual() {}
 
Individual& Individual::operator = (const Individual &b){
	//memcpy(chromosome,b.chromosome,cus_num);
	for (int j=0;j<cus_num;j++){
		chromosome[j]=b.chromosome[j];
	}
	for (int i=0;i<veh_num;i++){
		veh[i]=b.veh[i];
	}
	fitness=b.fitness;
	unfitness=b.unfitness;
}
 
bool Individual::is_satisfy(){
	for (int i=0;i<veh_num;i++){
			veh[i].clear();
	}
	bool isv[veh_num];
	memset(isv,0,veh_num);
	for (int i=0;i<cus_num;i++){
		//if (chromosome[i]<0||chromosome[i]>=veh_num) cout<<"no"<<endl;
		//cout<<i<<" remain "<<veh[chromosome[i]].cap_remain<<endl;
		isv[chromosome[i]]=1;
		if (veh[chromosome[i]].cap_remain>0)
		{
			if (!veh[chromosome[i]].push(i)) return 0;	
		}
		else return 0;
	}
	for (int i=0;i<veh_num;i++){
		if (!isv[i]) return 0;
	}
	return 1;
}
void Individual::update() {
		for (int i=0;i<veh_num;i++){
			veh[i].clear();
		}
		/*
		for (int i=0;i<cus_num;i++){
			if (veh[chromosome[i]].cap_remain>0)
				if (!veh[chromosome[i]].push(i)) return 0;
			else return 0;
			//veh[i].get_coor();
		}*/
		for (int i=0;i<cus_num;i++){
			//if (veh[chromosome[i]].cap_remain>0)
			//每次生成儿子不用判断是否合理 
			//veh[chromosome[i]].push(i);
			veh[chromosome[i]].cus_vec.push_back(i);
			veh[chromosome[i]].cap_remain-=cus[i].demand;
			//veh[i].get_coor();
		}
		for (int i=0;i<veh_num;i++){
			veh[i].get_coor();
		}
		get_ufit();
		get_fit();
		if (fitness>best.fitness&&unfitness==0){
			best=*this;
		}
		sort(veh,veh+veh_num);
		for (int i=0;i<veh_num;i++){
			for (int j=0;j<veh[i].cus_vec.size();j++){
				chromosome[veh[i].cus_vec[j]]=i;
			}
		}
		return ; 
}

double Individual::get_ufit() {
	double sum = 0;
	for (int i = 0; i < veh_num; i++) {
		//计算总超过的容量
		if (veh[i].cap_remain < 0) {
			sum -= veh[i].cap_remain;
		}
	}
	return unfitness = sum / CAP;
}
 
 
double Individual::get_fit() {
	double len = 0;
	for (int i = 0; i < veh_num; i++) {
		//double sum = 0;
		//int l = veh[i].cus_vec.size();
		veh[i].optimate();
		len+=veh[i].length;
	}
	return fitness = 1 / len;
}

// void draw_picture(){
// 	srcImage = cv::Mat::zeros(cv::Size(700, 700), CV_8UC1);
// 	cout << "drawPoint = " << drawPoint.size() << endl;
// 	for(int i=0 ; i<drawPoint.size() ; i++)
// 	{
// 		cv::circle(srcImage, drawPoint[i], 3, Scalar(255,255,255), -1, LINE_AA);
// 	}
// 	int i=0, j=0;
// 	// for (i=0;i<veh_num;i++)
// 	// {
// 		line(srcImage, drawPoint[0], drawPoint[best.veh[i].cus_vec[0]+1], Scalar(255,255,255), 3);
// 		for (j=0;j<best.veh[i].cus_vec.size()-1;j++)
// 		{
// 			// cout<<" "<<best.veh[i].cus_vec[j]<<" ";
// 			line(srcImage, drawPoint[best.veh[i].cus_vec[j]+1], drawPoint[best.veh[i].cus_vec[j]+2], Scalar(255,255,255), 7);
// 		};
// 		line(srcImage, drawPoint[best.veh[i].cus_vec[j]+2], drawPoint[0], Scalar(255,255,255), 7);
// 	// }
// 	imshow("点位图", srcImage);
// 	waitKey(0);
// }

/// 
int main(){
	if(!ini_data()){
		cout<<"have no file"<<endl;
		return 0;	
	}
	srand((unsigned)time(0)); 
	#ifdef DEBUG
	cout << "veh_num = " << veh_num << endl;
	cout << "cus_num = " << cus_num << endl;
	cout << "CAP = " << CAP << endl;
	cout << "optimal = " << opt_dis << endl;
	cout << "tightness = " << tightness << endl;
	#endif

	//计算出初始解
	ini_popul();

	//显示染色体序列
	#ifdef CHROM
		Individual tmp;
		tmp=popul[0];
		for (int j=0;j<cus_num;j++){
				cout<<tmp.chromosome[j]<<" ";
			}
			cout<<endl;
		
		for (int i=0;i<popul_size;i++){
			for (int j=0;j<cus_num;j++){
				cout<<popul[i].chromosome[j]<<" ";
			}
			cout<<endl;
		}
	#endif


	for (int i=0;i<veh_num;i++){
		for (int j=0;j<best.veh[i].cus_vec.size();j++){
			cout<<best.veh[i].cus_vec[j]<<" ";
		}
		cout<<endl;
	}
	if (best.fitness==0) cout<<"000"<<endl;
	else cout<<1/(best.fitness)<<endl;
	cout<<endl;
	cout<<endl;
	//cout<<Rc<<" "<<endl;
	evolution();
	
	//bool isv[33]={0};
	for (int i=0;i<veh_num;i++){
		for (int j=0;j<best.veh[i].cus_vec.size();j++){
			cout<<" "<<best.veh[i].cus_vec[j]<<" ";
		//	isv[best.veh[i].cus_vec[j]]=1;
		}
		cout<<endl;
		//cout<<"veh_num "<<i<<endl;
	}
	/*
	for (int i=0;i<cus_num;i++) {
		
		if (isv[i]==0) cout<<"false"<<endl; 
	}*/
	if (best.fitness==0) cout<<"000"<<endl;
	else cout<<1/(best.fitness)<<endl;
	// draw_picture();
} 
