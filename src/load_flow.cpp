//============================================================================
// Name        : load_flow2.cpp
// Author      : Yeap Yew Ming
// Version     : ver2
// Copyright   : Your copyright notice
// Description : Fast Decouple Load Flow Analysis
//============================================================================

//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
#include <iomanip>

using namespace std;
using std::getline;

bool InArray (float Haystack[], int size, float needle) // function to check if an element exists in array
{
  for (int i = 0; i < size; i++)
  {
    if (Haystack[i] == needle)
      return true;
  }

  return false;
}

int main (void)
{ int i=0, j=0, counter=0;
string line;
string a;
string array[20][20];
float array_int[20][20];
complex<float> Y[10][10],B[10][10],D[10][10],X[10][10],c[10];
float base=100;
float vmag[10],ang[10],pload[10],qload[10],pgen[10],qgen[10];
float Pmismatch[10]={0},Qmismatch[10]={0},newvoltage[10]={0},newangle[10]={0};
int busno[10],bustype[10];
float k[10][10]={0},bb[10][10]={0},dd[10][10]={0},cc[10][10]={0};
float maxerror=1;
int counter1,slack;
float pv[5];
int m1=1,m2=1,p=0;
float P=0,Psigma=0;
float Qq=0, Qsigma=0;
int l,n,e,q;
int ll=1;
float conv=0.01;
int iteration=0;
int no_of_bus=0;
int no_of_pv=1;


//open line data file
ifstream datafile;
datafile.open("ieee9_linedata.txt");

if (!datafile)
{
	cerr<< "Line data file is not found" <<endl;
	system("PAUSE");	
	exit (1);
}

string token;
stringstream iss;

while (getline(datafile,line))	//extract data from file into array
{	
  iss<<line;
  j=0;
  while (getline(iss,token,'\t'))
  {
	  array[i][j] = token;
	  a=array[i][j];
	  istringstream convert(a); //convert string to double

	  if (!(convert>>array_int[i][j]))
	  {
		  array_int[0][0]=0;
	  }

	  ++j;
  }
  iss.clear();
  ++i;
  counter++;
}

datafile.close();

//open bus data file
ifstream busdata;
busdata.open("ieee9_busdata.txt");

if (!busdata)
{
	cerr<< "Bus data file is not found" <<endl;
	system("PAUSE");	
	exit (1);
}

int b=1;

while(! busdata.eof())
{
	busdata>> busno[b]>>bustype[b] >> vmag[b]>> ang[b]>> pgen[b]>>qgen[b]>>pload[b]>>qload[b];
	b++;
}  


counter1=b-1; //how many bus?

cout<<"\n_________________________________________________________________________";
cout<<"\n                      **START OF POWER FLOW ANALYSIS**             "<<endl;
cout<<"\n_________________________________________________________________________"<<endl;

cout << "Please key in the stopping criterion (ENTER for default 0.01): ";

if (cin.get() == '\n')
{
	conv=0.01;
}
else
{
	cin>>conv;
}
cout<<endl;

//find no of bus
float max1=0, max2=0, bus_no=0;

for(int i=0; i<counter; i++)
{
	if (array_int[i][0]>max1)
		max1=array_int[i][0];
	
	if (array_int[i][1]>max2)
		max2=array_int[i][1];
}

if (max1>max2)
	bus_no=max1;
else
	bus_no=max2;

cout<<"No of bus: "<<bus_no<<endl<<endl;

//non-diagonal element of Y matrix
for (int k=0; k<counter; k++)
{
	for (int q=1; q<=bus_no; q++)
	{
		for (int u=1; u<=bus_no; u++)
		{
			if ((array_int[k][0]==q || array_int[k][0]==u) && (array_int[k][1]==u || array_int[k][1]==q))
			{
				X[q][u]=complex<float>(1,0)/-(complex<float>(array_int[k][2],array_int[k][3]));
				B[q][u]=complex<float>(0,array_int[k][4]);
				
			}

		}
	}
}

//diagonal element of Y matrix
for (int q=1; q<=bus_no; q++)
{
	for (int u=1; u<=bus_no; u++)
	{
		if (q==u)
		{
			for (int t=1; t<=bus_no; t++)
			{
				Y[q][q]=-X[q][t]+B[q][t]+Y[q][q];
			}
		}
		else
		{	
			Y[q][u]=X[q][u];
		}

	}	
}

//output Y matrix
cout<<"Y matrix: "<<endl;
for(int i=1;i<=bus_no;i++)
{
	for(int j=1;j<=bus_no;j++)
		cout<<Y[i][j]<<"\t";
	cout<<endl;
}
cout<<endl;

//find slack bus and pv bus
for (int k=1;k<=counter1;k++)
{
	if (bustype[k]==0)
	{
		slack=busno[k];
		no_of_bus++;
	}

	if (bustype[k]==2)
	{
		pv[no_of_pv]=busno[k];
		no_of_bus++;
		no_of_pv++;

	}
}

	//form FD matrix
for (int k=0; k<counter; k++)
{
	for (int q=1; q<=bus_no; q++)
	{
		for (int u=1; u<=bus_no; u++)
		{
			if ((array_int[k][0]==q || array_int[k][0]==u) && (array_int[k][1]==u || array_int[k][1]==q))
			{
				bb[q][u]=1/array_int[k][3];
									
			}

		}
	}
}

//diagonal element of fast decoupled matrix
for (int q=1; q<=bus_no; q++)
{
	for (int u=1; u<=bus_no; u++)
	{
		if (q==u)
		{
			for (int t=1; t<=bus_no; t++)
			{
				k[q][q]=k[q][q]+(-bb[q][t]);

			}
		}
		else
		{	
			k[q][u]=bb[q][u];

		}

	}
	
}


//determine B' matrix
for (int i=1;i<=bus_no-1;i++)
{
	for (int j=1;j<=bus_no-1;j++)
	{
		if (j>=slack && i<=slack)
			cc[i][j]=k[i][j+1];			//cc is B'
		if (i>=slack && j<=slack)
			cc[i][j]=k[i+1][j];
		if (i>=slack && j>=slack)
			cc[i][j]=k[i+1][j+1];
		if (i<slack && j<slack)
			cc[i][j]=k[i][j];

		}
	}

//output B' matrix
cout<<"B' matrix: "<<endl;
for(int i=1;i<=bus_no-1;i++)
{
	for(int j=1;j<=bus_no-1;j++)
	{
		cout<<cc[i][j]<<"  ";
	}
	cout<<endl;
}
cout<<endl;

//inverse B' matrix
n=bus_no-1;
for (int i=1;i<=n;i++)
{
	cc[i][i]=1.0/cc[i][i];
	for (int j=1;j<=n;j++)
		{
			if (j!=i)
				{
					cc[j][i] = cc[j][i] * cc[i][i];
					for (int k=1;k<=n;k++)
						{
							if (k!=i)
								{
									cc[j][k] = cc[j][k] - cc[j][i]*cc[i][k];
									if (j==n)
										{
											cc[i][k] = -cc[i][i]* cc[i][k];
										}
								}
						}
				}
		}
}

e = n-1;
for (int l = 1;l<=e;l++)
	{
		cc[n][l] = -cc[n][n]*cc[n][l];
	}	

cout<<endl;

//output inverse of B' matrix
cout<<"Inverse of B' matrix: "<<endl;
for(int i=1;i<=n;i++)
{
	for(int j=1;j<=n;j++)
	{
		cout<<cc[i][j]<<"  ";
	}
	cout<<endl;
}
cout<<endl;

//form B'' matrix
for (int k=0; k<counter; k++)
{
	for (int q=1; q<=bus_no; q++)
	{
		for (int u=1; u<=bus_no; u++)
		{
			if ((array_int[k][0]==q || array_int[k][0]==u) && (array_int[k][1]==u || array_int[k][1]==q))
			{
				bb[q][u]=array_int[k][3]/(pow(array_int[k][3],2)+pow(array_int[k][2],2));
		
			}
	
		}
	}
}

//clearing matrix k
for (int i=1;i<=bus_no;i++)
	k[i][i]=0;

//diagonal element of B'' matrix
for (int q=1; q<=bus_no; q++)
{
	for (int u=1; u<=bus_no; u++)
	{
		if (q==u)
		{
			for (int t=1; t<=bus_no; t++)
			{
				k[q][q]=k[q][q]+(-bb[q][t]);
			}
		}
		else
		{	
			k[q][u]=bb[q][u];
		}	
	}

}
//determine B'' matrix
for(int i=1;i<=bus_no;i++)
{

	for(int j=1;j<=bus_no;j++)
	{
		if (i!=slack && j!=slack && !InArray(pv,no_of_pv,i) && !InArray(pv,no_of_pv,j))
		{	dd[m1][m2]=k[i][j];				//dd is B''
			m2++;
	
			if (m2==((bus_no-no_of_bus)+1))
			{	m2=1;
				m1++;
			}
		}
	
	}
}

//outputof B'' matrix
cout<<"B'' matrix: "<<endl;
for(int i=1;i<=bus_no-no_of_bus;i++)
{
	for(int j=1;j<=bus_no-no_of_bus;j++)
		{cout<<dd[i][j]<<"  ";}
	cout<<endl;
}
cout<<endl;
	
//inverse B'' matrix
n=bus_no-no_of_bus;
for (int i=1;i<=n;i++)
{
	dd[i][i]=1.0/dd[i][i];
	for (int j=1;j<=n;j++)
		{
			if (j!=i)
				{
					dd[j][i] = dd[j][i] * dd[i][i];
					for (int k=1;k<=n;k++)
						{
							if (k!=i)
								{
									dd[j][k] = dd[j][k] - dd[j][i]*dd[i][k];
									if (j==n)
										{
											dd[i][k] = -dd[i][i]* dd[i][k];
										}
								}
						}
				}
		}
}

e = n-1;
for (int l = 1;l<=e;l++)
	{
		dd[n][l] = -dd[n][n]*dd[n][l];
	}	

cout<<endl;
	
//output inverse of B'' matrix
cout<<"Inverse of B'' matrix: "<<endl;
for(int i=1;i<=n;i++)
{
	for(int j=1;j<=n;j++)
	{
		cout<<dd[i][j]<<"  ";
	}
	cout<<endl;

}
cout<<endl;

//iteration starts here
while (maxerror>=conv)
{	
	//calculate power mismatch
   for(int i=1;i<=bus_no;i++)
   {
	   P=0;
		Psigma=0;
		
	   if(bustype[i]!=0)
	   {
		   for(int j=0;j<counter;j++)
		   {
			   if(array_int[j][0]==i || array_int[j][1]==i)
			   {
				   if(array_int[j][0]==i)
				   {p=array_int[j][1];}
				   if(array_int[j][1]==i)
				   {p=array_int[j][0];}

				   P=P+vmag[i]*vmag[p]*abs(Y[i][p])*cos(arg(Y[i][p])-ang[i]+ang[p]);

			   }
			   
		   }
		  
		   Psigma=pow(vmag[i],2)*abs(Y[i][i])*cos(arg(Y[i][i]))+P;
	
		
			if(bustype[i]==1)
			{	l=busno[i];
				Pmismatch[l]=pgen[i]/base-pload[i]/base-Psigma;
			}
		
			if(bustype[i]==2)
			{	l=busno[i];
				Pmismatch[l]=pgen[i]/base-pload[i]/base-Psigma;
			}
	
	   }

   }

	
	int lll=1;
	float change=0;
		   
	//compute the change in angle
	for (int i=1;i<=bus_no-1;i++)
	{	for (int j=1;j<=bus_no;j++)
		{
			if (bustype[j]!=0)
			{
				int pp=busno[j];
				change=change+cc[i][ll]*Pmismatch[pp]/vmag[j];
				ll++;

			}								
		}

		if (bustype[lll]==0)
			lll++;
		newangle[lll]=change;
		change=0;
		ll=1;
		lll++;
	}

   //update angle
   for (int i=1;i<=bus_no;i++)
   {	
	   if (bustype[i]!=0)
		ang[i]=ang[i]-newangle[i];
		//cout<<"angle at bus "<<i<<" is "<<ang[i]<<endl;
   }

 //calculate Q mismatch
	for(int i=1;i<=bus_no;i++)
	{
	   Qq=0;
	   Qsigma=0;
		
	   if(bustype[i]!=0 || bustype[i]!=2)
	   {
		   for(int j=0;j<counter;j++)
		   {
			   if(array_int[j][0]==i || array_int[j][1]==i)
			   {
				   if(array_int[j][0]==i)
				   {q=array_int[j][1];}
				   if(array_int[j][1]==i)
				   {q=array_int[j][0];}
	
				   Qq =Qq+vmag[i]*vmag[q]*abs(Y[i][q])*sin(arg(Y[i][q])-ang[i]+ang[q]);
	
			   }
			   
		   }
		  
		   Qsigma=-pow(vmag[i],2)*abs(Y[i][i])*sin(arg(Y[i][i]))-Qq;
	
			if(bustype[i]==1)
			{	
				l=busno[i];
				Qmismatch[l]=qgen[i]/base-qload[i]/base-Qsigma;
			}
	
	   }				   
	}

	//compute the change in voltage
	change=0,lll=1;
	for (int i=1;i<=bus_no-no_of_bus;i++)
	{	for (int j=1;j<=bus_no;j++)
		{
			if (bustype[j]==1)
			{
				int pp=busno[j];
				change=change+dd[i][ll]*Qmismatch[pp]/vmag[j];
				ll++;
	
			}								
		}
	
		if (bustype[lll]==0 || bustype[lll]==2)
			lll++;
		newvoltage[lll]=change;
		change=0;
		ll=1;
		lll++;
	 }
			

	//update new voltage
	for (int i=1;i<=bus_no;i++)
   {	
	   if (bustype[i]!=0 || bustype[i]!=2)
		vmag[i]=vmag[i]-newvoltage[i];
   }
 

	//find the max error
	max1=0;
   for(int i=0; i<counter; i++)
	 {
		if (abs(Pmismatch[i])>max1)
		max1=abs(Pmismatch[i]);
		
	}

	maxerror=max1;
	iteration=iteration+1;
	cout<<"maxerror at iteration "<< iteration<<" is: "<<maxerror<<endl;//
}
//end of iteration

//find power generation at slack bus
P=0,Qq=0;
for(int i=1;i<=bus_no;i++)
{
	if (bustype[i]==0)
	{
		for(int j=0;j<counter;j++)
			   {
				   if(array_int[j][0]==i || array_int[j][1]==i)
				   {
					   if(array_int[j][0]==i)
					   {p=array_int[j][1];}
					   if(array_int[j][1]==i)
					   {p=array_int[j][0];}

						P=P+(vmag[i]*vmag[p]*abs(Y[i][p])*cos(arg(Y[i][p])-ang[i]+ang[p]))*base;
						Qq=Qq+(vmag[i]*vmag[p]*abs(Y[i][p])*sin(arg(Y[i][p])-ang[i]+ang[p]))*base;

				   }
					pgen[i]=(pow(vmag[i],2)*abs(Y[i][i])*cos(arg(Y[i][i])))*base+P;
					qgen[i]=-(pow(vmag[i],2)*abs(Y[i][i])*sin(arg(Y[i][i])))*base-Qq;
			   }
	}
}

//find power generation at PV bus


for(int i=1;i<=bus_no;i++)
{	Qq=0;
	if (bustype[i]==2)
	{
		for(int j=0;j<counter;j++)
			   {
				   if(array_int[j][0]==i || array_int[j][1]==i)
				   {
					   if(array_int[j][0]==i)
					   {p=array_int[j][1];}
					   if(array_int[j][1]==i)
					   {p=array_int[j][0];}
	
						Qq=Qq+(vmag[i]*vmag[p]*abs(Y[i][p])*sin(arg(Y[i][p])-ang[i]+ang[p]))*base;
				   }
					qgen[i]=-(pow(vmag[i],2)*abs(Y[i][i])*sin(arg(Y[i][i])))*base-Qq+qload[i];
			   }
	}
}

cout<<endl;
cout<<"Program stops at iteration "<<iteration<<" with max error "<<maxerror<<endl<<endl;
cout<<"\n_________________________________________________________________________";
cout<<"\n                            Updated Bus Data             "<<endl;
cout<<"\n_________________________________________________________________________"<<endl;
cout<< setw(4) <<"  BUS  "<<"VOLTAGE  "<<"   ANGLE  "<<"     <<<<<<<GENERATION>>>>>>"<<"      <<<<<<<LOAD>>>>>>>"<<endl;
cout<< setw(4) <<"  NO. "<<"   MAG. "<<"     DEGREE "<<"        MW   "<<"    Mvar   "<<"       MW    "<<"     Mvar   "<<endl<<endl;
for (int i=1;i<=bus_no;i++)
cout<<setw(4)<<setprecision(4)<<setiosflags(ios::fixed)<<busno[i]<<"  "<<setprecision(4)<<setiosflags(ios::fixed)<<vmag[i]<<"     "<<setprecision(4)<<setiosflags(ios::fixed)<<ang[i]<<"        "<<setprecision(4)<<setiosflags(ios::fixed)<<pgen[i]<<"    "<<setprecision(4)<<setiosflags(ios::fixed)<<qgen[i]<<"       "<<setprecision(4)<<setiosflags(ios::fixed)<<pload[i]<<"       "<<setprecision(4)<<qload[i]<<endl;

//clear B,X,k,bb matrix
for (int i=1;i<=bus_no+1;i++)
{
	for (int i=1;i<=bus_no+1;i++)
	{
		D[i][j]=0;
		X[i][j]=0;
		k[i][j]=0;
		bb[i][j]=0;
	}
}

//calculate the load flow
for(int i=1;i<=bus_no;i++)
{
	c[i]=polar(vmag[i],ang[i]);
}

cout<<"______________________________________________________________________"<<endl;
cout<< setw(4) <<"    Bus  "<<"    To Bus  "<<"   Active power (MW)  "<<"  Reactive power (MVAR)"<<endl;
cout<<"______________________________________________________________________"<<endl;

for(int i=1;i<=bus_no;i++)
{
	for(int j=1;j<=bus_no;j++)
	{
		D[i][j]=(-Y[i][j])*(c[i]-c[j])+B[i][j]*c[i];
		X[i][j]=c[i]*conj(D[i][j]);
		k[i][j]=real(X[i][j])*base;
		bb[i][j]=imag(X[i][j])*base;
	
	}
}

for(int i=1;i<=bus_no;i++)
{
	for(int j=0;j<counter;j++)
		   {
			   if(array_int[j][0]==i || array_int[j][1]==i)
			   {
				   if(array_int[j][0]==i)
				   {p=array_int[j][1];}
				   if(array_int[j][1]==i)
				   {p=array_int[j][0];}

					cout<< setprecision(4)<< setiosflags(ios::fixed)  <<"     "<<i<<"\t       "<< p; 
					cout<< setprecision(4)<< setiosflags(ios::fixed)  <<"\t"<<k[i][p]<<"  \t      "<<bb[i][p]<<"\t   "<<endl;

			   }

		   }
		
}
cout<<"\n_________________________________________________________________________";
cout<<"\n                      **END OF ANALYSIS**             "<<endl;
cout<<"\n_________________________________________________________________________"<<endl;
    
    
system("PAUSE");
return 0;
}

