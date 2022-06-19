#include <iostream>
#include <string>
#include <fstream>
//#include <cmath>
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

/////////////////////////
////	Impedance	////
///////////////////////
struct Impedance{
	double rel, img;
	std::complex<double> z;
	std::complex<double> y;
	Impedance(){
		rel = 0;
		img = 0;
		z.real(rel);
		z.imag(img);
	}

	std::complex<double> getAdmittance(){
		Impedance conj;
		double d = (rel * rel) + (img * img);
		conj.rel = rel/d;
		conj.img = -img/d;
		
		y.real(conj.rel);
		y.imag(conj.img);
		return y;

	}

	std::complex<double> getImpedance(){
		z.real(rel);
		z.imag(img);
		return z;
	}

	Impedance operator + (Impedance z){
		Impedance s;
		s.rel = rel + z.rel;
		s.img = img + z.img;
		return s;
	}

	Impedance operator * (Impedance z){
		Impedance s;
		s.rel = (this->rel * z.rel - this->img * z.img);
		s.img = (this->rel * z.img + this->img * z.rel);
		return s;
	}
	
	Impedance operator ||(Impedance z){
		Impedance p;
		Impedance s = *this + z;
		Impedance m = *this * z;
		Impedance conj;
		conj.rel = s.rel;
		conj.img = -s.img;
		double d = (s.rel * s.rel) + (s.img * s.img);
		p.rel = m.rel/d;
		p.img = m.img/d;
		p = p * conj;
		return p;
	}
}nul;

/////////////////////////
////	Voltage		////
///////////////////////
struct Voltage{
	double rel, img;
	std::complex<double> V;
	
	Voltage(){
		rel = 0;
		img = 0;
		V.real(rel);
		V.imag(img);
	}

	Voltage operator + (Voltage z){
		Voltage s;
		s.rel = rel + z.rel;
		s.img = img + z.img;
		return s;
	}

	Voltage operator * (Voltage z){
		Voltage s;
		s.rel = (this->rel * z.rel - this->img * z.img);
		s.img = (this->rel * z.img + this->img * z.rel);
		return s;
	}

	std::complex<double> getVoltage(){
		V.real(rel);
		V.imag(img);
		return V;

	}

	void setVoltage(std::complex<double> c)
	{
		V = c;
		rel = V.real();
		img = V.imag();
	}
};

/////////////////////////
////	Current		////
///////////////////////
struct Current{
	double rel, img;
	std::complex<double> I;
	
	Current(){
		rel = 0;
		img = 0;
		I.real(rel);
		I.imag(img);
	}

	Current operator + (Current z){
		Current s;
		s.rel = rel + z.rel;
		s.img = img + z.img;
		return s;
	}

	Current operator * (Current z){
		Current s;
		s.rel = (this->rel * z.rel - this->img * z.img);
		s.img = (this->rel * z.img + this->img * z.rel);
		return s;
	}

	std::complex<double> getCurrent(){
		I.real(rel);
		I.imag(img);
		return I;

	}
	void setcurrnt(std::complex<double> c)
	{
		I=c;
		rel=I.real();
		img=I.imag();
	}



};

/////////////////////////
////	Component	////
///////////////////////
class Component{
protected:
	int node1, node2;

public:
	Component(int n1, int n2){ setNodes(n1, n2); }
	void setNodes(int n1, int n2){ node1 = n1;		node2 = n2; }
	int getNode1(){ return node1; }
	int getNode2(){ return node2; }
	virtual Impedance getZ() = 0;
};

/////////////////////////////////////////
////	Independent Voltage Source	////
///////////////////////////////////////
class IndepVolSrc : public Component{
protected:
	double Vmax;
	double phi;
	double angularFreq;
	Voltage value;
public:
	IndepVolSrc(double v, double w, double a, int n1, int n2) : Component(n1, n2), value() 
	{ Vmax = v;	angularFreq = w;		phi = a; 
	value.rel +=getVmax() * cos(getPhi() * (22.0 / 7) / 180);
	value.img += getVmax() * sin(getPhi() * (22.0 / 7) / 180);
	}
	IndepVolSrc(int n1, int n2): Component(n1, n2){ }
	double getVmax(){ return Vmax; }
	double getPhi(){ return phi; }
	double getOmiga(){ return angularFreq; }
	Impedance getZ(){ return nul; }
	Voltage getValue() { return value; }
};

/////////////////////////////////////////
////	Dependent Voltage Source	////
///////////////////////////////////////
class DepVolSrc : public IndepVolSrc{
	double coeff;
	int depN1, depN2;				//	Node1 & Node2 of The Source We Depend On
public:
	DepVolSrc(double c, int compN1, int compN2, int depN1, int depN2) : IndepVolSrc(depN1, depN2){
		coeff = c;
		this->depN1 = depN1;
		this->depN2 = depN2;
	}

	int getDepN1(){ return depN1; }
	int getDepN2(){ return depN2; }
	void setSrc(double Vm, double angle, double w){
		Vmax = (Vm > 0)? coeff*Vm : -coeff*Vm;
		phi = angle;
		angularFreq = w;
	}
};

/////////////////////////////////////////
////	Independent Current Source	////
///////////////////////////////////////
class IndepCrntSrc : public Component{
protected:
	double Imax;
	double phi;
	double angularFreq;
	Current value;
public:
	IndepCrntSrc(double i, double w, double a, int n1, int n2) : Component(n1, n2),value()
	{ Imax = i;	angularFreq = w;	phi = a;
	value.rel += getImax() * cos(getPhi() * (22.0 / 7.0) / 180.0);
	value.img += getImax() * sin(getPhi() * (22.0 / 7.0) / 180.0);
	}
	IndepCrntSrc(int n1, int n2): Component(n1, n2){ }
	double getImax(){ return Imax; }
	double getPhi(){ return phi; }
	double getOmiga(){ return angularFreq; }
	Impedance getZ(){ return nul; }
	Current getValue() { return value; }
};

/////////////////////////////////////////
////	Dependent Current Source	////
///////////////////////////////////////
class DepCrntSrc : public IndepCrntSrc{
	double coeff;
	int depN1, depN2;				//	Node1 & Node2 of The Source We Depend On
public:
	DepCrntSrc(double c, int compN1, int compN2, int depN1, int depN2) : IndepCrntSrc(depN1, depN2){
		coeff = c;
		this->depN1 = depN1;
		this->depN2 = depN2;
	}

	int getDepN1(){ return depN1; }
	int getDepN2(){ return depN2; }
	void setSrc(double Im, double angle, double w){
		Imax = (Im > 0)? coeff*Im : -coeff*Im;
		phi = angle;
		angularFreq = w;
	}
};

/////////////////////////
////	Resistor	////
///////////////////////
class Resistor : public Component{
	Impedance z;

public:
	Resistor(double r, int n1, int n2) : Component(n1, n2){ z.rel = r;	z.img = 0; }
	Impedance getZ(){ return z; }
};

/////////////////////////
////	Capacitor	////
///////////////////////
class Capacitor : public Component{
	Impedance z;

public:
	Capacitor(double c, double w, int n1, int n2) : Component(n1, n2){ z.rel = 0;	z.img = -1/(w*c); }
	Impedance getZ(){ return z; }
};

/////////////////////////
////	Inductor	////
///////////////////////
class Inductor : public Component{
	Impedance z;

public:
	Inductor(double l, double w, int n1, int n2) : Component(n1, n2){ z.rel = 0;	z.img = w*l; }
	Impedance getZ(){ return z; }
};


class Node
{
public:
	int num;
	Voltage V;
	int i;
	Node() { i = 0; }
	void setV(std::complex<double> v) { V.setVoltage(v); }
	void setN(int n) { num = n; }
	int getN() { return num; }
	std::complex<double> getV() { return V.getVoltage(); }
};

class branch
{
public:
	Impedance impBranch;
	Voltage volt;
	Current amb;
	int node1;
	int node2;
	Current BranchCurrent;
	branch():impBranch(),volt(),amb()
	{
	}
	void printinfo()
	{
		cout<<"Branch ( "<<node1<<" , "<<node2<<" )"<<endl;
		cout<<"Impedance = ";
		cout<<impBranch.getImpedance();
		cout<<endl;
		cout<<"Voltage = ";
		cout<<volt.getVoltage();
		cout<<endl;
		cout<<"Current = ";
		cout<<amb.getCurrent();
		cout<<endl;
	}
	bool isv()
	{
		if(volt.rel==0 && volt.img==0)
			return false;
		else return true;
	}
	bool isimp()
	{
		if(impBranch.rel==0 && impBranch.img==0)
			return false;
		else return true;
	}


};


class realBranch
{
public:
	Node* N1, * N2;
	Impedance impBranch;
	Current BranchCurrent;
	Voltage volt;
	Component* BranchComponents[10];
	int BranchComponentCount = 0;
	realBranch() :impBranch(), volt(), BranchCurrent()
	{
		impBranch.rel = 0;
		impBranch.img = 0;
		volt.setVoltage(0);
	}
	void setN1(Node* n)
	{
		N1 = n;
	}

	void setN2(Node* n)
	{
		N2 = n;
	}

	Node* getN1()
	{
		return N1;
	}

	Node* getN2()
	{
		return N2;
	}

	void AddComponent(Component* c)
	{
		BranchComponents[BranchComponentCount++] = c;
	}

	void calculateCurrent()
	{
		bool flag = false;
		for(int i = 0; i< BranchComponentCount; i++)
		{
			if (dynamic_cast<IndepCrntSrc*>(BranchComponents[i]) != NULL)
			{
				flag = true;
				/////work here
				BranchCurrent = ((IndepCrntSrc*)BranchComponents[i])->getValue();
			}
			else
			{
				impBranch = impBranch + BranchComponents[i]->getZ();
			}
		}
		if (!flag)
		{
			////work here 
			if(abs(N1->getV()) > abs(N2->getV()))
				BranchCurrent.setcurrnt((N1->getV() - N2->getV()) / impBranch.getImpedance());
			else
				BranchCurrent.setcurrnt((N2->getV() - N1->getV()) / impBranch.getImpedance());
		}
	}


	Current getCurrent()
	{
		return BranchCurrent;
	}

	bool IsInside(Component* c)
	{
		for (int i = 0; i < BranchComponentCount; i++)
		{
			if (BranchComponents[i] == c)
				return true;
		}
		return false;
	}
	void printInfo()
	{
		cout << "I(";
		cout << getN1()->getN();
		for (int i = 0; i < BranchComponentCount; i++)
		{
			cout << "," << BranchComponents[i]->getNode2() ;
		}
		cout << ") = "<<BranchCurrent.getCurrent()<<endl;
	}
};

struct nonsimpleknown
{
	int node;
	Voltage vb;
	int n;
	nonsimpleknown()
	{
		node=-1;
		vb.img=0;
		vb.rel=0;
		n=0;
	}

};

struct nonsimplnode
{
	int node;
	int index;
	nonsimplnode()
	{
		node=0;
		index=0;
	}

};

void deletenode(int);

void voltknown();

bool supernode(int& , int&, complex<double>&);

void CalculateNodes(Node arrN[], int n_N, realBranch arrB[], int n_B, Component* arrC[], int n_C);

bool isNonSimpleNode(int *arr, int n, int node);

void analyseComp(Component* pComp, int n1, int n2);

void nodeanalysis();

void VoltToCurrent();

void LoadInputFile(Component* [], int&, string);

class FILE_NOT_FOUND {};

nonsimpleknown nonSknown[5];
int num_nonSknown=0;

nonsimplnode nsnode[6];
int n_nsnode=0;

branch B[20];
Node N[10];
realBranch RB[20];
int n_B = 0;
int n_RB = 0;
Impedance branchImp[10];				//	Contains The Impedance Of Each Branch
Voltage branchVol[10];					//	Contains The Value Of The Voltage Source In Each Branch -->> if exist
Current branchCrnt[10];					//	Contains The Value Of The Current Source In Each Branch -->> if exist
int x = 0;								//	Index Of The Branch
int numOfConnectedBranchesToNode[5];	//	����� ���� ���� ���
int nonSimpleNodes[6], numOfNonSimpNodes = 0;
int numOfNodes = 2;
int main(){

	Component* arr[10];
	int num = 0;

	string InputFileName;
	while (true)
	{
		cout << "Please, enter the name of the file to " << endl;;
		cin >> InputFileName;
		try
		{
			LoadInputFile(arr, num, InputFileName);
			break;
		}
		catch (FILE_NOT_FOUND)
		{
			cout << "file not found" << endl;
		}
	}
	/*IndepCrntSrc i1(2,0,0,0,1);
	arr[num++]= &i1;

	Resistor r1(2,0,1);
	arr[num++]= &r1;

	IndepVolSrc v1(2,0,0,1,2);
	arr[num++]= &v1;

	Resistor r2(10,1,2);
	arr[num++]= &r2;

	Resistor r3(4,2,0);
	arr[num++]= &r3;

	IndepCrntSrc i2(7,0,0,2,0);
	arr[num++]= &i2;*/
	//	Getting The Non-Simple Nodes
	
	N[0].setN(arr[0]->getNode1());
	N[1].setN(arr[0]->getNode2());
	bool flag1, flag2;

	for (int i = 1; i < num; i++)
	{
		flag1 = flag2 = false;
		for (int k = 0; k < numOfNodes; k++)
		{
			if (N[k].getN() == arr[i]->getNode1())
			{
				flag1 = true;
			}

			if (N[k].getN() == arr[i]->getNode2())
			{
				flag2 = true;
			}
		}
		if (!flag1)
			N[numOfNodes++].setN(arr[i]->getNode1());
		else if (!flag2)
			N[numOfNodes++].setN(arr[i]->getNode2());
	}

	for(int i = 0; i < num; i++){
		int n1, n2, n3, n4, repeatN1 = 0, repeatN2 = 0;
		n1 = arr[i]->getNode1();
		n2 = arr[i]->getNode2();
		for(int j = i; j < num; j++){
			n3 = arr[j]->getNode1();
			n4 = arr[j]->getNode2();
			if(n1 == n3 || n1 == n4)
				repeatN1++;
			if(n2 == n3 || n2 == n4)
				repeatN2++;
		}
		if(repeatN1 > 2){
			if(!isNonSimpleNode(nonSimpleNodes, numOfNonSimpNodes, n1))
				nonSimpleNodes[numOfNonSimpNodes++] = n1;
		}
		if(repeatN2 > 2){
			if(!isNonSimpleNode(nonSimpleNodes, numOfNonSimpNodes, n2))
				nonSimpleNodes[numOfNonSimpNodes++] = n2;
		}

		//	Setting The Dependant Sources
		if(dynamic_cast<DepVolSrc*>(arr[i]) != NULL){
			DepVolSrc* v = dynamic_cast<DepVolSrc*>(arr[i]);
			int depN1 = v->getDepN1(),
				depN2 = v->getDepN2();

			for(int k = 0; k < num; k++){
				if(((arr[k]->getNode1() == depN1 && arr[k]->getNode2() == depN2) || (arr[k]->getNode1() == depN2 && arr[k]->getNode2() == depN1)) && dynamic_cast<IndepVolSrc*>(arr[k]) != NULL){
					IndepVolSrc* v2 = dynamic_cast<IndepVolSrc*>(arr[k]);
					v->setSrc(v2->getVmax(), v2->getPhi(), v2->getOmiga());
					break;
				}
			}
		}
		if(dynamic_cast<DepCrntSrc*>(arr[i]) != NULL){
			DepCrntSrc* crnt = dynamic_cast<DepCrntSrc*>(arr[i]);
			int depN1 = crnt->getDepN1(),
				depN2 = crnt->getDepN2();

			for(int k = 0; k < num; k++){
				if(((arr[k]->getNode1() == depN1 && arr[k]->getNode2() == depN2) || (arr[k]->getNode1() == depN2 && arr[k]->getNode2() == depN1)) && dynamic_cast<IndepCrntSrc*>(arr[k]) != NULL){
					IndepCrntSrc* crnt2 = dynamic_cast<IndepCrntSrc*>(arr[k]);
					crnt->setSrc(crnt2->getImax(), crnt2->getPhi(), crnt2->getOmiga());
					break;
				}
			}
		}
	}

	int refNode = nonSimpleNodes[0];
	
	for(int i = 1; i < numOfNonSimpNodes; i++)
	{

		for(int j = 0; j < num; j++)
		{
			int n1 = arr[j]->getNode1();
			int n2 = arr[j]->getNode2();
			int startNode;

			if(n2 == nonSimpleNodes[i])
			{
				B[n_B].node1=n2;
				numOfConnectedBranchesToNode[i]++;
				startNode = n1;
				B[n_B].node2=startNode;
				analyseComp(arr[j], n1, n2);

				while(!isNonSimpleNode(nonSimpleNodes, numOfNonSimpNodes, startNode))
				{
					for(int k = 0; k < num; k++){
						if(k == j)
							continue;
						n1 = arr[k]->getNode1();
						n2 = arr[k]->getNode2();
						if(startNode == n2)
						{
							startNode = n1;
							B[n_B].node2=startNode;
							analyseComp(arr[k], n1, n2);
							break;
						}
						else if(startNode == n1)
						{
							startNode = n2;
							B[n_B].node2=startNode;
							analyseComp(arr[k], n1, n2);
							break;
						}
					}
				}
				B[n_B++].node2=startNode;
				x++;
			}
			else if(n1 == nonSimpleNodes[i])
			{
				B[n_B].node1=n1;
				numOfConnectedBranchesToNode[i]++;
				startNode = n2;
				B[n_B].node2=startNode;
				analyseComp(arr[j], n1, n2);

				while(!isNonSimpleNode(nonSimpleNodes, numOfNonSimpNodes, startNode))
				{
					for(int k = 0; k < num; k++){
						if(k == j)
							continue;
						n1 = arr[k]->getNode1();
						n2 = arr[k]->getNode2();
						if(startNode == n2){
							startNode = n1;
							B[n_B].node2=startNode;
							analyseComp(arr[k], n1, n2);
							break;
						}else if(startNode == n1){
							startNode = n2;
							B[n_B].node2=startNode;
							analyseComp(arr[k], n1, n2);
							break;
						}
					}
				}
				B[n_B++].node2=startNode;
				x++;
			}
		}
	}
	VoltToCurrent();
	//voltknown();

	for(int i = 0; i < numOfNonSimpNodes; i++)
		cout << nonSimpleNodes[i] << endl;


	RB[n_RB++].setN1(&N[0]);
	RB[0].setN2(&N[1]);
	RB[0].AddComponent(arr[0]);
	int PreviousNode;
	flag1 = flag2 = false;
	
	for (int i = 1; i < num; i++)
	{
		for (int k = 0; k < n_RB; k++)
		{
			flag1 = false;
			PreviousNode = RB[k].getN2()->getN();
			if (PreviousNode == arr[i]->getNode1() && !isNonSimpleNode(nonSimpleNodes, numOfNonSimpNodes, PreviousNode))
			{
				RB[k].AddComponent(arr[i]);
				RB[k].setN2(&N[arr[i]->getNode2()]);
			}
			else 
				flag1 = true;
		}
			if(flag1)
			{
			RB[n_RB++].setN1(&N[arr[i]->getNode1()]);
			RB[n_RB - 1].setN2(&N[arr[i]->getNode2()]);
			RB[n_RB - 1].AddComponent(arr[i]);
			}
	}

	cout << "\nImpedence\n";
	for(int i = 0; i < x; i++){
		cout << branchImp[i].getImpedance() << endl;
	}
	cout << "\nAdmittance\n";
	for(int i = 0; i < x; i++){
		cout << branchImp[i].getAdmittance() << endl;
	}
	cout << "\nCurrentSource\n";
	for(int i = 0; i < x; i++){
		cout << branchCrnt[i].getCurrent() << endl;
	}
	cout << "VoltageSource\n";
	for(int i = 0; i < x; i++){
		cout << branchVol[i].getVoltage() << endl;
	}

	for(int i = 1; i < numOfNonSimpNodes; i++){
		cout << numOfConnectedBranchesToNode[i] << endl;
	}

	/*for(int i=0 ;i<numOfConnectedBranchesToNode[2];i++)
		cout << branchImp[i+3].getImpedance() << endl;*/
	for(int i=0;i<n_B;i++)
	{
		B[i].printinfo();
		cout<<endl;
	}

	for(int i=0; i<numOfNonSimpNodes; i++)
	{
		nsnode[n_nsnode].node=nonSimpleNodes[i];
		nsnode[n_nsnode].index=i;
	}

	nodeanalysis();

	for (int i = 0; i < n_RB; i++)
	{
		RB[i].calculateCurrent();
		RB[i].printInfo();
	}

	CalculateNodes(N, numOfNodes, RB, n_RB, arr, num);
	cout << "V0 = 0" << endl;
	for (int i = 1; i < numOfNodes; i++)
	{
		cout << "V" << i << " = " << N[i].getV()<<endl;
	}



	system("pause");
	return 0;
}

bool isNonSimpleNode(int *arr, int n, int node){
	for(int i = 0; i < n; i++){
		if(arr[i] == node)
			return true;
	}
	return false;
}

void analyseComp(Component* pComp, int n1, int n2){
	if(dynamic_cast<Resistor*>(pComp) != NULL){
		Resistor* r = dynamic_cast<Resistor*>(pComp);
		branchImp[x] = branchImp[x] + r->getZ();
		B[n_B].impBranch=branchImp[x];
		
	}
	else if(dynamic_cast<Capacitor*>(pComp) != NULL){
		Capacitor* c = dynamic_cast<Capacitor*>(pComp);
		branchImp[x] = branchImp[x] + c->getZ();
		B[n_B].impBranch=branchImp[x];
	}
	else if(dynamic_cast<Inductor*>(pComp) != NULL){
		Inductor* l = dynamic_cast<Inductor*>(pComp);
		branchImp[x] = branchImp[x] + l->getZ();
		B[n_B].impBranch=branchImp[x];
	}
	else if(dynamic_cast<IndepVolSrc*>(pComp) != NULL){
		IndepVolSrc* v = dynamic_cast<IndepVolSrc*>(pComp);
		branchVol[x].rel += v->getVmax()*cos(v->getPhi() * (22.0/7)/180);
		branchVol[x].img += v->getVmax()*sin(v->getPhi() * (22.0/7)/180);
		
		if(B[n_B].node2==n1)
		{
			B[n_B].volt.rel=branchVol[x].rel;
			B[n_B].volt.img=branchVol[x].img;
		}
		else
		{
			B[n_B].volt.rel=-branchVol[x].rel;
			B[n_B].volt.img=-branchVol[x].img;
		}
	
		if(n1 > n2){
			branchVol[x].rel *= -1;
			branchVol[x].img *= -1;
			
		}

		
	}else{
		IndepCrntSrc* v = dynamic_cast<IndepCrntSrc*>(pComp);
		branchCrnt[x].rel += v->getImax()*cos(v->getPhi() * (22.0/7.0)/180.0);
		branchCrnt[x].img += v->getImax()*sin(v->getPhi() * (22.0/7.0)/180.0);
		if(B[n_B].node2==n2)
		{
			B[n_B].amb.rel=-branchCrnt[x].rel;
			B[n_B].amb.img=-branchCrnt[x].img;
		}
		else
		{
			B[n_B].amb.rel=branchCrnt[x].rel;
			B[n_B].amb.img=branchCrnt[x].img;
		}
		if(n1 > n2)
		{
			branchCrnt[x].rel *= -1;
			branchCrnt[x].img *= -1;
			
		}
	}
}


void nodeanalysis()
{
	if (numOfNonSimpNodes != 0)
	{
		int n = numOfNonSimpNodes - 1;
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> y(n, n);
		y.setZero();
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cu(n, 1);
		cu.setZero();
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> R(n, 1);
		R.setZero();
		for (int i = 1; i < numOfNonSimpNodes; i++)
		{
			for (int j = 0; j < n_B; j++)
			{
				if (B[j].node1 == nonSimpleNodes[i] && B[j].node2 == 0)
				{
					if (B[j].isimp())
						y(i - 1, i - 1) = y(i - 1, i - 1) + 1.0 / B[j].impBranch.getImpedance();
					cu(i - 1, 0) = cu(i - 1, 0) + B[j].amb.getCurrent();
				}
				else if (B[j].node1 == nonSimpleNodes[i] && B[j].node2 != 0)
				{


					for (int k = 1; k < numOfNonSimpNodes; k++)
					{
						if (B[j].node2 == nonSimpleNodes[k] && i != k)
						{
							if (B[j].isimp())
							{
								y(i - 1, i - 1) = y(i - 1, i - 1) + 1.0 / B[j].impBranch.getImpedance();
								y(i - 1, k - 1) = y(i - 1, k - 1) - 1.0 / B[j].impBranch.getImpedance();
							}
							cu(i - 1, 0) = cu(i - 1, 0) + B[j].amb.getCurrent();
						}

					}

				}

			}

		}

		voltknown();
		for (int t = 0; t < num_nonSknown; t++)
		{
			int m = nonSknown[t].n;
			cu.col(0).array() -= y.col(m - 1).array();
			cu.col(0).array() *= nonSknown[t].vb.getVoltage();
			y.col(m - 1).setZero();
			y.row(m - 1).setZero();
			y(m - 1, m - 1) = 1;
			cu(m - 1, 0) = nonSknown[t].vb.getVoltage();


		}

		{
			int i = 0; int k = 0; complex<double> e;
			if (supernode(i, k, e))
			{
				y.row(i - 1) = y.row(i - 1) + y.row(k - 1);
				y.row(k - 1).setZero();
				cu.row(i - 1) = cu.row(i - 1) + cu.row(k - 1);
				cu.row(k - 1).setZero();
				y(k - 1, k - 1) = 1;
				if (k < i)
				{
					y(k - 1, i - 1) = -1;
					cu(k - 1, 0) = e;
				}
				else
				{
					y(k - 1, i - 1) = -1;
					cu(k - 1, 0) = -e;
				}

			}
		}





		R = y.inverse() * cu;

		for (int i = 1; i < numOfNonSimpNodes; i++)
		{
			for (int k = 0; k < numOfNodes; k++)
				if (N[k].getN() == nonSimpleNodes[i])
				{
					N[k].setV(R(i - 1, 0));
				}
		}

		cout << y << endl << endl;
		cout << cu << endl << endl;
		cout << R << endl << endl;

	}
}

void VoltToCurrent()
{
	for(int i=0; i<n_B ;i++)
	{
		if(B[i].isv() && B[i].isimp())
		{
			B[n_B].node1=B[i].node1;
			B[n_B].node2=B[i].node2;
			
			B[n_B].amb.setcurrnt(B[i].volt.getVoltage()/B[i].impBranch.getImpedance());
			
			n_B++;
			B[i].volt.rel=0;
			B[i].volt.img=0;
		}

	}

}


void voltknown()
{
	for(int i=1; i<6 ; i++)
	{
		for (int j=0 ; j<n_B ;j++)
		{
			if( B[j].node1==nonSimpleNodes[i]  && B[j].node2==0 && !B[j].isimp() && B[j].isv())
			{
				nonSknown[num_nonSknown].node=B[j].node1;
				nonSknown[num_nonSknown].n=i;
				nonSknown[num_nonSknown].vb=B[j].volt;
				num_nonSknown++;
				/*nonSimpleNodes[i]=0;
				for(int k=i; k<5 ;k++)
				{
					if(nonSimpleNodes[k+1]!=0)
					{
						nonSimpleNodes[k]=nonSimpleNodes[k+1];
						nonSimpleNodes[k+1]=0;
					}

				}
				numOfNonSimpNodes--;*/

			}

		}

	}

}


void deletenode( int i)
{
	nonSimpleNodes[i]=0;
	for(int k=i; k<5 ;k++)
	{
		if(nonSimpleNodes[k+1]!=0)
		{
			nonSimpleNodes[k]=nonSimpleNodes[k+1];
			nonSimpleNodes[k+1]=0;
		}

	}
	numOfNonSimpNodes--;
}

bool supernode( int& p , int& q , complex<double>& c)
{
	for(int i=1; i<numOfNonSimpNodes; i++)
	{
		for(int j=0;j<n_B;j++)
		{
			if(B[j].node1==nonSimpleNodes[i] && B[j].node2!=0)
			{
				for(int k=1; k<numOfNonSimpNodes; k++)
				{
					if(B[j].node2==nonSimpleNodes[k] && i!=k)
					{
						if(!B[j].isimp() && B[j].isv())
						{
							p=i;
							q=k;
							c=B[j].volt.getVoltage();
							/*deletenode(i);
							deletenode(k);*/
							return true;
						}

					}

				}

			}
		}
	}
	return false;
}

void LoadInputFile(Component* arr[], int& N, string FileName)
{
	ifstream InputFile("Input\\" + FileName + ".txt");
	if (InputFile.is_open())
	{
		string ComponentType;
		string ComponentName;
		int  n1, n2, depN1, depN2;
		double value, phi;
		int Omega = 0;
		while (!InputFile.eof())
		{
			InputFile >> ComponentType;

			if (ComponentType == "w")
			{
				InputFile >> value;
				Omega = value;
			}
			else if (ComponentType == "res")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value;
				arr[N++] = new Resistor(value, n1, n2);
			}
			else if (ComponentType == "vsrc")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value >> phi;
				arr[N++] = new IndepVolSrc(value, Omega, phi, n1, n2);
			}
			else if (ComponentType == "vcvs")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value >> depN1 >> depN2;		//	Value here refer to the Coefficient
				arr[N++] = new DepVolSrc(value, n1, n2, depN1, depN2);
			}
			else if (ComponentType == "isrc")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value >> phi;
				arr[N++] = new IndepCrntSrc(value, Omega, phi, n1, n2);
			}
			else if (ComponentType == "cccs")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value >> depN1 >> depN2;		//	Value here refer to the Coefficient
				arr[N++] = new DepCrntSrc(value, n1, n2, depN1, depN2);
			}
			else if (ComponentType == "cap")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value;
				arr[N++] = new Capacitor(value, Omega, n1, n2);
			}
			else if (ComponentType == "ind")
			{
				InputFile >> ComponentName >> n1 >> n2 >> value;
				arr[N++] = new Inductor(value, Omega, n1, n2);
			}
			else if (ComponentType == "-1")
			{
				InputFile.close();
				break;
			}
		}
	}
	else
		throw FILE_NOT_FOUND();
}


void CalculateNodes(Node arrN[], int n_N, realBranch arrB[], int n_B, Component* arrC[], int n_C)
{
	if (numOfNonSimpNodes != 0) {
		for (int i = 1; i < n_N; i++)
		{
			for (int k = 0; k < n_C; k++)
			{
				if (arrC[k]->getNode2() == arrN[i].getN())
				{
					for (int j = 0; j < n_B; j++)
					{
						if (arrB[j].IsInside(arrC[k]))
						{
							if (dynamic_cast<IndepVolSrc*>(arrC[k]) == NULL && dynamic_cast<IndepCrntSrc*>(arrC[k]) == NULL)
							{
								Voltage PreviousV;
								PreviousV.setVoltage(0);
								for (int l = 0; l < numOfNodes; l++)
								{
									if (arrN[i].getN() == arrN[l].getN() + 1)
									{
										PreviousV.setVoltage(arrN[l].getV());
									}
								}
								arrN[i].setV(PreviousV.getVoltage() - (arrB[j].getCurrent().getCurrent() * arrC[k]->getZ().getImpedance()));
								break;
							}
							else if (dynamic_cast<IndepVolSrc*>(arrC[k]) != NULL)
							{
								arrN[i].setV(((IndepVolSrc*)arrC[k])->getValue().getVoltage());
								break;
							}
						}
					}

				}
			}
		}
	}
	else
	{
		double R = 0;
		double v = 0;
		for (int k = 0; k < n_C; k++)
		{
			if (dynamic_cast<IndepVolSrc*>(arrC[k]) != NULL)
				v += ((IndepVolSrc*)arrC[k])->getValue().getVoltage().real();
		}
		for (int k = 0; k < n_C; k++)
		{
			if (dynamic_cast<IndepVolSrc*>(arrC[k]) == NULL)
				R += arrC[k]->getZ().rel;
		}

		double C = (v / R);
		for (int i = 1; i < n_C; i++)
		{
			for (int k = 0; k < n_C; k++)
			{
				N[k].setV(arrC[k]->getZ().rel * C);
			}
		}

	}
}